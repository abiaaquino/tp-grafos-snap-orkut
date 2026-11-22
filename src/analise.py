"""Orquestração das Partes I–III. Cada etapa grava figuras (figuras/), tabelas
(resultados/*.csv e *.tex) e um JSON (resultados/json/<etapa>.json); a etapa
``consolidar`` junta tudo em resultados/resultados.json, resultados.md e
valores.tex (macros para o relatório). Usado pelo notebook e por scripts/run_all.py.
"""
from __future__ import annotations

import json
import math
import time

import numpy as np
import pandas as pd
from scipy import stats

from . import algoritmos as A
from . import leipotencia as LP
from . import metricas as M
from . import modelos as Mo
from . import robustez as R
from .centralidade import brandes, pagerank
from .distancias import caminho_medio_amostrado, para_igraph, todos_para_todos, validar_com_bfs_propria
from .figuras import CORES, plt, salvar
from .medicao import (ajuste_linear_nm, ajuste_proporcional, cronometrar, info_maquina,
                      intervalo_confianca)
from .relatorio import RES, escrever_macros, ler_json, salvar_json, salvar_tabela
from .tratamento import PROC, carregar_lcc, com_pesos, pesos_embeddedness, subgrafo_bfs

SEED = 42


class Contexto:
    """Grafo carregado + objetos derivados reaproveitados entre etapas."""

    def __init__(self, rapido: bool = False):
        self.rapido = rapido
        self.G, self.nos = carregar_lcc()
        self.w = pesos_embeddedness(self.G)
        self.Gw = com_pesos(self.G, self.w)
        self.ig = para_igraph(self.G)
        self.maquina = info_maquina()
        self.cache = {}
        RES.mkdir(exist_ok=True)


def preparar(rapido: bool = False) -> Contexto:
    return Contexto(rapido)


def _cronometrado(fn):
    t0 = time.perf_counter()
    r = fn()
    return r, time.perf_counter() - t0


# ============================================================== PARTE I
def parte1(ctx: Contexto) -> dict:
    G, nos = ctx.G, ctx.nos
    graus = G.graus()
    rg = M.resumo_graus(graus)
    comp = M.componentes(G.adj)
    dist, t_dist = _cronometrado(lambda: todos_para_todos(G, g=ctx.ig))
    ctx.cache["dist"] = dist
    valida = validar_com_bfs_propria(G, np.random.default_rng(SEED).choice(G.n, 20, replace=False))
    cl, t_cl = _cronometrado(lambda: M.clusterizacao(G.adj))
    ctx.cache["clus"] = cl
    # conferência independente (igraph, em C)
    trans_ig = ctx.ig.transitivity_undirected()
    cmed_ig = ctx.ig.transitivity_avglocal_undirected(mode="zero")
    dens = M.densidade(G.n, G.m)

    res = {"vertices": G.n, "arestas": G.m, **rg, "densidade": dens,
           "n_componentes": comp["n_componentes"], "tamanhos_componentes": comp["tamanhos"][:10],
           "diametro": dist["diametro"], "raio": dist["raio"], "caminho_medio": dist["caminho_medio"],
           "pares_ordenados": dist["pares_ordenados"], "distancias_exatas": dist["exato"],
           "hist_distancias": dist["hist_distancias"].tolist(),
           "segundos_distancias": t_dist, "validacao_bfs_propria_20_fontes": valida,
           "triangulos": cl["triangulos"], "clusterizacao_media": cl["clusterizacao_media"],
           "transitividade": cl["transitividade"], "segundos_triangulos": t_cl,
           "conferencia_igraph": {"transitividade": trans_ig, "clusterizacao_media": cmed_ig},
           "centro_tamanho": int((dist["excentricidades"] == dist["raio"]).sum()),
           "periferia_tamanho": int((dist["excentricidades"] == dist["diametro"]).sum())}

    # ---- tabelas
    t = pd.DataFrame([
        ("Vértices |V|", G.n), ("Arestas |E|", G.m), ("Grau mínimo", rg["grau_min"]),
        ("Grau máximo", rg["grau_max"]), ("Grau médio", rg["grau_medio"]), ("Grau mediano", rg["grau_mediano"]),
        ("Densidade", dens), ("Componentes conexas", comp["n_componentes"]),
        ("Tamanho da maior componente", comp["maior"]), ("Diâmetro", dist["diametro"]),
        ("Raio", dist["raio"]), ("Comprimento médio dos caminhos", dist["caminho_medio"]),
        ("Triângulos", cl["triangulos"]), ("Clusterização média", cl["clusterizacao_media"]),
        ("Transitividade (global)", cl["transitividade"]),
    ], columns=["Medida", "Valor"])
    t["Valor"] = t["Valor"].astype(object)
    salvar_tabela(t, "parte1_medidas", "Medidas estruturais do grafo G27 (com-Orkut, LCC tratada).",
                  "tab:parte1", casas=4, alinhamento="lr")
    tc = pd.DataFrame({"Componente": np.arange(1, min(10, comp["n_componentes"]) + 1),
                       "Vértices": comp["tamanhos"][:10]})
    salvar_tabela(tc, "parte1_componentes", "Tamanho das componentes conexas.", "tab:componentes", casas=0)
    th = pd.DataFrame({"Distância": np.arange(1, len(dist["hist_distancias"])),
                       "Pares ordenados": dist["hist_distancias"][1:],
                       "Fração (%)": 100 * dist["hist_distancias"][1:] / dist["pares_ordenados"]})
    salvar_tabela(th, "parte1_distancias", "Distribuição das distâncias entre pares de vértices.",
                  "tab:distancias", casas={"Distância": 0, "Pares ordenados": 0, "Fração (%)": 2})

    calc = pd.DataFrame([
        ("Vértices e arestas", "sim", "contagem direta da lista de arestas"),
        ("Grau mín./máx./médio/mediano", "sim", "vetor de graus (CSR)"),
        ("Distribuição de graus", "sim", "histograma log-log e CCDF (k = 0 excluído)"),
        ("Densidade", "sim", "2m / (n(n-1))"),
        ("Componentes conexas e tamanhos", "sim", "BFS própria"),
        ("Diâmetro, raio, caminho médio", "sim", "exato: BFS de todos os vértices (igraph, C), em blocos; conferido com BFS própria"),
        ("Triângulos", "sim", "implementação própria, orientação por grau O(m·√m); conferida com igraph"),
        ("Clusterização média e transitividade", "sim", "a partir dos triângulos por vértice"),
        ("Visualização reduzida", "sim", "150 vértices de maior grau, layout Kamada-Kawai"),
        ("Componentes fortemente conexas / Euler dirigido", "não", "não se aplica: com-Orkut é não direcionado"),
        ("Floyd-Warshall no grafo inteiro", "não", "O(n³) ≈ 9,8e11 passos e matriz de 787 MB; rodado só em subgrafos"),
        ("Pesos reais nas arestas", "não", "o SNAP não fornece pesos; usamos peso unitário e um peso derivado (vizinhos em comum)"),
    ], columns=["Medida", "Calculada", "Método / justificativa"])
    salvar_tabela(calc, "parte1_calculadas", "Medidas calculadas e não calculadas, com o método.",
                  "tab:calculadas", alinhamento="lcp{8cm}")

    # ---- figuras
    k, pk, cc = M.distribuicao_graus(graus)
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    ax[0].hist(graus, bins=60, color=CORES["azul"])
    ax[0].set(xlabel="grau k", ylabel="nº de vértices", title="Histograma (escala linear)")
    ax[1].loglog(k, pk, "o", ms=3, color=CORES["azul"])
    ax[1].set(xlabel="grau k", ylabel="P(k)", title="Distribuição de graus (log-log)")
    ax[2].loglog(k, cc, "-", color=CORES["laranja"])
    ax[2].set(xlabel="grau k", ylabel="P(K ≥ k)", title="CCDF (log-log)")
    salvar(fig, "p1_distribuicao_graus")

    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    h = dist["hist_distancias"]
    ax.bar(np.arange(1, len(h)), 100 * h[1:] / h[1:].sum(), color=CORES["azul"])
    ax.axvline(dist["caminho_medio"], color=CORES["vermelho"], ls="--", label=f"média = {dist['caminho_medio']:.2f}")
    ax.set(xlabel="distância (saltos)", ylabel="% dos pares", title="Distribuição das distâncias")
    ax.legend()
    salvar(fig, "p1_distribuicao_distancias")

    _visualizacao(ctx)
    salvar_json(res, "parte1")
    return res


def _visualizacao(ctx, top: int = 150):
    import networkx as nx
    G, nos = ctx.G, ctx.nos
    graus = G.graus()
    sel = np.argsort(-graus, kind="stable")[:top]
    mapa = {int(x): i for i, x in enumerate(sel)}
    H = nx.Graph()
    H.add_nodes_from(range(len(sel)))
    for x in sel:
        for y in G.adj[int(x)]:
            if y in mapa and mapa[y] > mapa[int(x)]:
                H.add_edge(mapa[int(x)], mapa[y])
    pos = nx.kamada_kawai_layout(H)
    com = nos["comunidade_principal"].to_numpy()[sel]
    cats = {c: i for i, c in enumerate(pd.Series(com).value_counts().index)}
    cmap = plt.get_cmap("tab10")
    cor = [cmap(cats[c] % 10) if cats[c] < 9 else (0.6, 0.6, 0.6, 1) for c in com]
    fig, ax = plt.subplots(figsize=(8, 8))
    nx.draw_networkx_edges(H, pos, ax=ax, alpha=0.12, width=0.5)
    nx.draw_networkx_nodes(H, pos, ax=ax, node_size=20 + 200 * graus[sel] / graus.max(), node_color=cor,
                           linewidths=0.3, edgecolors="k")
    ax.set_title(f"{top} vértices de maior grau ({H.number_of_edges()} arestas entre eles);\n"
                 "cor = comunidade principal (9 mais frequentes; cinza = demais)")
    ax.axis("off")
    salvar(fig, "p1_visualizacao_top150")


# ============================================================== PARTE II
COMPLEXIDADE = {"BFS": "O(n+m)", "DFS": "O(n+m)", "Euler (verificação)": "O(n+m)",
                "Dijkstra": "O((n+m) log n)", "Bellman-Ford": "O(n·m) [efetivo O(k·m)]",
                "Tarjan (pontes/articulações)": "O(n+m)", "Prim": "O(m log n)", "Kruskal": "O(m log m)",
                "Floyd-Warshall (numpy)": "O(n³)", "Floyd-Warshall (Python)": "O(n³)"}


def _funcoes(Gw):
    wadj = Gw.adj_ponderada()
    arcos = A.arcos_nao_direcionados(Gw.u.tolist(), Gw.v.tolist(), Gw.w.tolist())
    ul, vl, wl = Gw.u.tolist(), Gw.v.tolist(), Gw.w.tolist()
    s = 0
    return {
        "BFS": lambda: A.bfs(Gw.adj, s),
        "DFS": lambda: A.dfs(Gw.adj, s),
        "Euler (verificação)": lambda: A.verificar_euler(Gw.adj),
        "Dijkstra": lambda: A.dijkstra(wadj, s),
        "Bellman-Ford": lambda: A.bellman_ford(Gw.n, arcos, s),
        "Tarjan (pontes/articulações)": lambda: A.tarjan_pontes_articulacoes(Gw.adj),
        "Prim": lambda: A.prim(wadj, s),
        "Kruskal": lambda: A.kruskal(Gw.n, ul, vl, wl),
    }


def _f_teorica(nome, n, m, rodadas=None):
    if nome in ("BFS", "DFS", "Euler (verificação)", "Tarjan (pontes/articulações)"):
        return n + m
    if nome in ("Dijkstra",):
        return (n + m) * math.log2(n)
    if nome == "Prim":
        return m * math.log2(n)
    if nome == "Kruskal":
        return m * math.log2(m)
    if nome == "Bellman-Ford":
        return rodadas * 2 * m
    if nome.startswith("Floyd"):
        return n ** 3
    raise KeyError(nome)


def parte2(ctx: Contexto) -> dict:
    G, Gw, w = ctx.G, ctx.Gw, ctx.w
    reps = 5 if ctx.rapido else 30
    reps_esc = 3 if ctx.rapido else 10
    res = {"maquina": ctx.maquina, "pesos": {
        "definicao": "w(u,v) = 1/(1 + |N(u) ∩ N(v)|) (embeddedness); peso unitário usado para validar contra BFS",
        "w_min": float(w.min()), "w_mediana": float(np.median(w)), "w_max": float(w.max()),
        "frac_w_igual_1": float((w == 1).mean())}}

    # ---------- execuções e validações no grafo inteiro
    s = 0
    dist_bfs, pai = A.bfs(G.adj, s)
    ref_ig = ctx.ig.distances(source=[s])[0]
    ordem_dfs, pai_dfs = A.dfs(G.adj, s)
    prof = [0] * G.n
    for x in ordem_dfs[1:]:
        prof[x] = prof[pai_dfs[x]] + 1
    d_unit, _ = A.dijkstra(com_pesos(G, np.ones(G.m)).adj_ponderada(), s)
    d_w, pai_w = A.dijkstra(Gw.adj_ponderada(), s)
    bf, rodadas, neg = A.bellman_ford(G.n, A.arcos_nao_direcionados(G.u.tolist(), G.v.tolist(), w.tolist()), s)
    hops = []
    for x in range(G.n):
        h, y = 0, x
        while pai_w[y] >= 0:
            y, h = pai_w[y], h + 1
        hops.append(h)
    euler = A.verificar_euler(G.adj)
    pontes, artic = A.tarjan_pontes_articulacoes(G.adj)
    import networkx as nx
    Hnx = nx.Graph(list(zip(G.u.tolist(), G.v.tolist())))
    ok_tarjan = (set(pontes) == {tuple(sorted(e)) for e in nx.bridges(Hnx)}
                 and artic == set(nx.articulation_points(Hnx)))
    arv_p, tot_p = A.prim(Gw.adj_ponderada(), 0)
    arv_k, tot_k = A.kruskal(G.n, G.u.tolist(), G.v.tolist(), w.tolist())
    g2 = para_igraph(Gw)
    mst_ig = g2.spanning_tree(weights="weight")
    tot_ig = float(sum(mst_ig.es["weight"]))
    grau_pontes = [len(G.adj[a]) + len(G.adj[b]) for a, b in pontes]
    res["execucao"] = {
        "fonte": s, "fonte_id_original": int(G.ids[s]),
        "bfs": {"alcancados": sum(d >= 0 for d in dist_bfs), "excentricidade_fonte": max(dist_bfs),
                "vertices_por_nivel": np.bincount(dist_bfs).tolist(),
                "confere_igraph": list(map(int, ref_ig)) == dist_bfs},
        "dfs": {"alcancados": len(ordem_dfs), "profundidade_max_arvore": max(prof)},
        "euler": {k: v for k, v in euler.items() if k != "impares_exemplo"},
        "dijkstra_unitario_igual_bfs": all(abs(d_unit[x] - dist_bfs[x]) < 1e-9 for x in range(G.n)),
        "dijkstra_ponderado": {"dist_max": max(d_w), "dist_media": float(np.mean(d_w)),
                               "saltos_medios_no_caminho_ponderado": float(np.mean(hops)),
                               "saltos_max": int(max(hops))},
        "bellman_ford": {"rodadas_ate_convergir": rodadas, "limite_teorico_n_menos_1": G.n - 1,
                         "ciclo_negativo": neg,
                         "igual_dijkstra": all(abs(bf[x] - d_w[x]) < 1e-9 for x in range(G.n))},
        "tarjan": {"pontes": len(pontes), "pontos_articulacao": len(artic), "confere_networkx": ok_tarjan,
                   "pontes_com_folha": sum(1 for a, b in pontes if min(len(G.adj[a]), len(G.adj[b])) == 1)},
        "agm": {"peso_prim": tot_p, "peso_kruskal": tot_k, "peso_igraph": tot_ig, "arestas": len(arv_k),
                "prim_igual_kruskal": abs(tot_p - tot_k) < 1e-6, "confere_igraph": abs(tot_k - tot_ig) < 1e-6,
                "frac_arestas_peso_1_na_agm": float(np.mean([e[2] == 1.0 for e in arv_k]))},
    }

    # ciclo negativo: basta UMA aresta negativa num grafo não direcionado (u-v-u)
    sub = subgrafo_bfs(Gw, 1000, seed=SEED)
    wneg = sub.w.copy()
    wneg[0] = -0.5
    _, rod_neg, tem_neg = A.bellman_ford(sub.n, A.arcos_nao_direcionados(sub.u.tolist(), sub.v.tolist(), wneg.tolist()), 0)
    res["execucao"]["bellman_ford_demo_ciclo_negativo"] = {"n": sub.n, "m": sub.m, "aresta_negativa": -0.5,
                                                            "detectou": tem_neg, "rodadas": rod_neg}

    # ---------- tempos no grafo inteiro (n >= 30 -> IC Normal)
    funcs = _funcoes(Gw)
    linhas = []
    tempos_brutos = {}
    for nome, f in funcs.items():
        t = cronometrar(f, reps=reps)
        tempos_brutos[nome] = t.tolist()
        ic = intervalo_confianca(t)
        linhas.append({"Algoritmo": nome, "Complexidade": COMPLEXIDADE[nome], "Repetições": ic["n"],
                       "Média (ms)": 1e3 * ic["media"], "DP (ms)": 1e3 * ic["dp"],
                       "IC95 inf. (ms)": 1e3 * ic["ic_inf"], "IC95 sup. (ms)": 1e3 * ic["ic_sup"],
                       "Método do IC": ic["metodo"]})
    tab = pd.DataFrame(linhas)
    salvar_tabela(tab, "parte2_tempos_grafo_completo",
                  f"Tempo de execução no grafo completo (n = {G.n:,}, m = {G.m:,}); {reps} repetições, IC de 95%.".replace(",", "."),
                  "tab:tempos", casas={"Repetições": 0}, alinhamento="llrrrrrl")
    res["tempos_grafo_completo"] = linhas

    # ---------- escala em subgrafos BFS (n < 30 repetições -> t-Student)
    tamanhos = [1000, 2000, 4000, 6000, 8000, G.n] if not ctx.rapido else [1000, 3000, G.n]
    esc = []
    for nalvo in tamanhos:
        H = subgrafo_bfs(Gw, nalvo, seed=SEED)
        fH = _funcoes(H)
        _, rod_H, _ = A.bellman_ford(H.n, A.arcos_nao_direcionados(H.u.tolist(), H.v.tolist(), H.w.tolist()), 0)
        for nome, f in fH.items():
            t = cronometrar(f, reps=reps_esc)
            ic = intervalo_confianca(t)
            esc.append({"Algoritmo": nome, "n": H.n, "m": H.m, "rodadas_BF": rod_H, "Repetições": ic["n"],
                        "Média (ms)": 1e3 * ic["media"], "DP (ms)": 1e3 * ic["dp"],
                        "IC95 inf. (ms)": 1e3 * ic["ic_inf"], "IC95 sup. (ms)": 1e3 * ic["ic_sup"],
                        "Método do IC": ic["metodo"]})
    esc = pd.DataFrame(esc)
    salvar_tabela(esc.drop(columns=["rodadas_BF"]), "parte2_tempos_escala",
                  f"Tempo em subgrafos BFS crescentes; {reps_esc} repetições por tamanho, IC de 95% (t-Student).",
                  "tab:escala", casas={"n": 0, "m": 0, "Repetições": 0}, alinhamento="lrrrrrrrl")

    # ajustes t ≈ c·f(n,m) e t ≈ a·n + b·m
    ajustes = []
    for nome in funcs:
        d = esc[esc["Algoritmo"] == nome]
        f = np.array([_f_teorica(nome, r.n, r.m, r.rodadas_BF) for r in d.itertuples()])
        aj = ajuste_proporcional(d["Média (ms)"], f)
        lin = ajuste_linear_nm(d["Média (ms)"], d["n"], d["m"])
        ajustes.append({"Algoritmo": nome, "f(n,m)": COMPLEXIDADE[nome], "c": aj["c"], "R² (c·f)": aj["r2"],
                        "a (ms/vértice)": lin["a"], "b (ms/aresta)": lin["b"], "R² (a·n+b·m)": lin["r2"]})
    aj = pd.DataFrame(ajustes)
    salvar_tabela(aj, "parte2_ajustes", "Ajuste dos tempos medidos à complexidade teórica.", "tab:ajustes",
                  casas={"c": 6, "R² (c·f)": 4, "a (ms/vértice)": 6, "b (ms/aresta)": 6, "R² (a·n+b·m)": 4},
                  alinhamento="llrrrrr")
    res["escala"] = esc.to_dict("records")
    res["ajustes"] = ajustes

    fig, axs = plt.subplots(2, 4, figsize=(15, 6.5))
    for ax, nome in zip(axs.ravel(), funcs):
        d = esc[esc["Algoritmo"] == nome]
        x = d["n"].to_numpy()
        ax.errorbar(x, d["Média (ms)"], yerr=[d["Média (ms)"] - d["IC95 inf. (ms)"], d["IC95 sup. (ms)"] - d["Média (ms)"]],
                    fmt="o", ms=4, capsize=3, color=CORES["azul"], label="medido (IC95)")
        c = aj.loc[aj["Algoritmo"] == nome, "c"].item()
        f = np.array([_f_teorica(nome, r.n, r.m, r.rodadas_BF) for r in d.itertuples()])
        ax.plot(x, c * f, "--", color=CORES["laranja"], label=f"c·{COMPLEXIDADE[nome].split(' [')[0][2:-1]}")
        ax.set(title=nome, xlabel="n (vértices do subgrafo)", ylabel="tempo (ms)")
        ax.legend(fontsize=7)
    salvar(fig, "p2_escala_tempos")

    res["floyd_warshall"] = _floyd(ctx, reps_esc)
    salvar_json(res, "parte2")
    salvar_json(tempos_brutos, "parte2_tempos_brutos")
    return res


def _floyd(ctx, reps):
    Gw = ctx.Gw
    # validação: FW (numpy e Python) = Dijkstra de todas as fontes num subgrafo
    H = subgrafo_bfs(Gw, 150, seed=SEED)
    D = A.matriz_pesos(H.n, H.u, H.v, H.w)
    F1 = A.floyd_warshall_numpy(D)
    F2 = np.array(A.floyd_warshall_python(D))
    wadj = H.adj_ponderada()
    DJ = np.array([A.dijkstra(wadj, s)[0] for s in range(H.n)])
    valid = {"n": H.n, "numpy_igual_python": bool(np.allclose(F1, F2)),
             "igual_dijkstra_todas_fontes": bool(np.allclose(F1, DJ))}
    tam_np = [200, 400, 600, 800, 1000, 1200] if not ctx.rapido else [200, 400]
    tam_py = [50, 100, 150, 200, 250] if not ctx.rapido else [50, 100]
    linhas = []
    for nome, tams, fn in (("Floyd-Warshall (numpy)", tam_np, A.floyd_warshall_numpy),
                           ("Floyd-Warshall (Python)", tam_py, A.floyd_warshall_python)):
        for nalvo in tams:
            H = subgrafo_bfs(Gw, nalvo, seed=SEED)
            D = A.matriz_pesos(H.n, H.u, H.v, H.w)
            t = cronometrar(lambda: fn(D), reps=reps)
            ic = intervalo_confianca(t)
            linhas.append({"Versão": nome, "n": H.n, "m": H.m, "Repetições": ic["n"],
                           "Média (ms)": 1e3 * ic["media"], "DP (ms)": 1e3 * ic["dp"],
                           "IC95 inf. (ms)": 1e3 * ic["ic_inf"], "IC95 sup. (ms)": 1e3 * ic["ic_sup"],
                           "Método do IC": ic["metodo"]})
    df = pd.DataFrame(linhas)
    salvar_tabela(df, "parte2_floyd_warshall",
                  "Floyd-Warshall em subgrafos BFS (pesos de embeddedness), IC de 95% (t-Student).",
                  "tab:floyd", casas={"n": 0, "m": 0, "Repetições": 0}, alinhamento="lrrrrrrrl")
    N = ctx.G.n
    extra = {}
    fig, ax = plt.subplots(figsize=(6, 4))
    for nome, cor in (("Floyd-Warshall (numpy)", CORES["azul"]), ("Floyd-Warshall (Python)", CORES["laranja"])):
        d = df[df["Versão"] == nome]
        aj = ajuste_proporcional(d["Média (ms)"] / 1e3, d["n"].to_numpy(float) ** 3)
        extra[nome] = {"c_s_por_n3": aj["c"], "r2": aj["r2"], "estimativa_grafo_completo_s": aj["c"] * N ** 3}
        ax.errorbar(d["n"], d["Média (ms)"] / 1e3, yerr=(d["IC95 sup. (ms)"] - d["Média (ms)"]) / 1e3, fmt="o",
                    color=cor, capsize=3, label=nome)
        xs = np.linspace(d["n"].min(), d["n"].max(), 50)
        ax.plot(xs, aj["c"] * xs ** 3, "--", color=cor, label=f"c·n³ (R² = {aj['r2']:.3f})")
    ax.set(xlabel="n", ylabel="tempo (s)", yscale="log", title="Floyd-Warshall: tempo × n")
    ax.legend(fontsize=7)
    salvar(fig, "p2_floyd_warshall")
    return {"validacao": valid, "medicoes": linhas, "ajustes": extra,
            "memoria_matriz_grafo_completo_MB": N * N * 8 / 1e6,
            "operacoes_grafo_completo": N ** 3}


# ============================================================== PARTE III
def small_world(ctx: Contexto) -> dict:
    G = ctx.G
    dist = ctx.cache.get("dist") or todos_para_todos(G, g=ctx.ig)
    cl = ctx.cache.get("clus") or M.clusterizacao(G.adj)
    C, L = cl["clusterizacao_media"], dist["caminho_medio"]
    kmed = 2 * G.m / G.n
    n_er = 2 if ctx.rapido else 5
    n_cf = 1 if ctx.rapido else 3
    kf = 300 if ctx.rapido else 1000
    inst = []
    for i in range(n_er):
        E = Mo.erdos_renyi_gnm(G.n, G.m, seed=SEED + i)
        inst.append({"modelo": "Erdős–Rényi G(n,m)", "seed": SEED + i,
                     "C": M.clusterizacao(E.adj)["clusterizacao_media"],
                     "L": caminho_medio_amostrado(E, k=kf, seed=i), "componentes": M.componentes(E.adj)["n_componentes"]})
    for i in range(n_cf):
        T = Mo.troca_dupla_arestas(G, trocas_por_aresta=5, seed=SEED + i)
        inst.append({"modelo": "Configuração (trocas duplas, graus preservados)", "seed": SEED + i,
                     "C": M.clusterizacao(T.adj)["clusterizacao_media"],
                     "L": caminho_medio_amostrado(T, k=kf, seed=i), "componentes": M.componentes(T.adj)["n_componentes"],
                     "trocas_realizadas": T.trocas[0]})
    df = pd.DataFrame(inst)
    resumo = []
    out = {"C": C, "L": L, "k_medio": kmed, "instancias": inst,
           "analitico_er": Mo.referencias_analiticas_er(G.n, kmed), "fontes_amostradas_L_nulos": kf}
    for nome, d in df.groupby("modelo", sort=False):
        sw = Mo.small_world(C, L, d["C"].mean(), d["L"].mean(), kmed)
        sw.update({"instancias": len(d), "C_rand_dp": float(d["C"].std(ddof=1)) if len(d) > 1 else 0.0,
                   "L_rand_dp": float(d["L"].std(ddof=1)) if len(d) > 1 else 0.0})
        out[nome] = sw
        resumo.append({"Referência": nome, "C": C, "C_ref": sw["C_rand"], "C/C_ref": sw["gamma_C_sobre_Crand"],
                       "L": L, "L_ref": sw["L_rand"], "L/L_ref": sw["lambda_L_sobre_Lrand"],
                       "σ": sw["sigma"], "ω": sw["omega"], "Small-world?": sw["small_world"]})
    salvar_tabela(pd.DataFrame(resumo), "parte3_small_world",
                  "Small-world: grafo observado × referências aleatórias (implementação própria).",
                  "tab:smallworld", casas=4, alinhamento="lrrrrrrrrl")
    df.to_csv(RES / "parte3_small_world_instancias.csv", index=False)

    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    nomes = ["Orkut G27"] + [r["Referência"].split(" (")[0] for r in resumo]
    ax[0].bar(nomes, [C] + [r["C_ref"] for r in resumo], color=[CORES["azul"], CORES["cinza"], CORES["verde"]])
    ax[0].set(yscale="log", ylabel="clusterização média C", title="C (escala log)")
    ax[1].bar(nomes, [L] + [r["L_ref"] for r in resumo], color=[CORES["azul"], CORES["cinza"], CORES["verde"]])
    ax[1].set(ylabel="caminho médio L", title="L")
    for a in ax:
        a.tick_params(axis="x", labelsize=8)
    salvar(fig, "p3_small_world")
    salvar_json(out, "parte3_small_world")
    return out


def lei_potencia(ctx: Contexto) -> dict:
    G, nos = ctx.G, ctx.nos
    sub = LP.analisar(graus=G.graus())
    comp = pd.read_csv(PROC / "orkut_completo_distribuicao_graus.csv")
    full = LP.analisar(k=comp["k"], c=comp["n_vertices"], min_cauda=100)
    out = {}
    linhas = []
    for nome, r in (("Subgrafo G27 (LCC, 9.917 vértices)", sub), ("Orkut completo (3.072.441 vértices)", full)):
        out[nome] = {k: v for k, v in r.items() if not k.startswith("_")}
        linhas.append({"Grafo": nome, "γ regressão bruta": r["regressao_bruta"]["gamma"],
                       "R² bruta": r["regressao_bruta"]["r2"], "γ bins log": r["regressao_bins_log"]["gamma"],
                       "R² bins": r["regressao_bins_log"]["r2"], "γ CCDF": r["regressao_ccdf"]["gamma"],
                       "γ MLE": r["mle"]["gamma"], "k_min": r["mle"]["kmin"], "cauda": r["mle"]["n_cauda"],
                       "KS": r["mle"]["ks"], "vs exp.": r["vs_exponencial"]["favorece"],
                       "p (Vuong)": r["vs_exponencial"]["p_valor"], "vs Poisson": r["vs_poisson"]["favorece"]})
    salvar_tabela(pd.DataFrame(linhas), "parte3_lei_potencia",
                  "Estimativas do expoente da lei de potência P(k) ~ k^(-γ).", "tab:powerlaw",
                  casas={"γ regressão bruta": 3, "R² bruta": 3, "γ bins log": 3, "R² bins": 3, "γ CCDF": 3,
                         "γ MLE": 3, "k_min": 0, "cauda": 0, "KS": 4, "p (Vuong)": 4},
                  alinhamento="p{3cm}rrrrrrrrrlrl")

    fig, axs = plt.subplots(2, 2, figsize=(11, 8))
    for linha, (nome, r) in enumerate((("Subgrafo G27", sub), ("Orkut completo", full))):
        k, c = r["_k"], r["_c"]
        p = c / c.sum()
        ax = axs[linha, 0]
        ax.loglog(k, p, ".", ms=3, color=CORES["cinza"], alpha=0.6, label="P(k) bruto")
        bx, by = r["_bins"]
        ax.loglog(bx, by, "o", color=CORES["azul"], label="bins logarítmicos")
        g, km = r["mle"]["gamma"], r["mle"]["kmin"]
        xs = np.geomspace(km, k.max(), 50)
        frac = c[k >= km].sum() / c.sum()
        from scipy.special import zeta
        ax.loglog(xs, frac * xs ** (-g) / zeta(g, km), "--", color=CORES["vermelho"],
                  label=f"MLE γ = {g:.2f} (k ≥ {km:.0f})")
        lam = r["vs_poisson"]["lambda"]
        kk = np.arange(1, k.max() + 1)
        pois = stats.poisson.pmf(kk, lam)
        ax.loglog(kk[pois > 1e-12], pois[pois > 1e-12], ":", color=CORES["verde"], label=f"Poisson λ = {lam:.1f}")
        ax.set(xlabel="k", ylabel="P(k)", title=f"{nome}: distribuição de graus", ylim=(max(p.min() / 10, 1e-12), 1))
        ax.legend(fontsize=7)
        ax = axs[linha, 1]
        cc = LP.ccdf(k, c)
        ax.loglog(k, cc, "-", color=CORES["azul"], label="CCDF empírica")
        ax.loglog(xs, frac * zeta(g, xs) / zeta(g, km), "--", color=CORES["vermelho"],
                  label=f"lei de potência γ = {g:.2f}")
        ax.loglog(kk, stats.poisson.sf(kk - 1, lam), ":", color=CORES["verde"], label="Poisson")
        ax.set(xlabel="k", ylabel="P(K ≥ k)", title=f"{nome}: CCDF", ylim=(max(cc.min() / 10, 1e-9), 1.5))
        ax.legend(fontsize=7)
    salvar(fig, "p3_lei_potencia")

    # grau no subgrafo × grau original
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.loglog(nos["grau_original"], G.graus(), ".", ms=2, alpha=0.4, color=CORES["azul"])
    ax.loglog([1, 3e4], [1, 3e4], "k--", lw=0.8, label="y = x")
    ax.set(xlabel="grau no Orkut completo", ylabel="grau no subgrafo", title="Efeito da amostragem no grau")
    ax.legend()
    salvar(fig, "p3_grau_original_vs_subgrafo")
    out["correlacao_spearman_grau_original_vs_sub"] = float(stats.spearmanr(nos["grau_original"], G.graus())[0])
    out["fracao_media_grau_preservado"] = float((G.graus() / nos["grau_original"]).mean())
    salvar_json(out, "parte3_lei_potencia")
    return out


def centralidades(ctx: Contexto) -> dict:
    if "central" in ctx.cache:
        return ctx.cache["central"]
    G = ctx.G
    k = 200 if ctx.rapido else 1000
    pr, it = pagerank(G)
    bt_amostra, t_bt = _cronometrado(lambda: brandes(G.adj, fontes=list(
        np.random.default_rng(SEED).choice(G.n, k, replace=False))))
    bt_exato_ig = np.array(ctx.ig.betweenness(directed=False)) / ((G.n - 1) * (G.n - 2) / 2)
    pr_ig = np.array(ctx.ig.pagerank(damping=0.85))
    top = 100
    a = set(np.argsort(-bt_amostra)[:top]); b = set(np.argsort(-bt_exato_ig)[:top])
    c = {"pagerank": pr, "betweenness": bt_amostra, "grau": G.graus().astype(float),
         "info": {"pagerank_iteracoes": it, "pagerank_confere_igraph_maxdif": float(np.abs(pr - pr_ig).max()),
                  "brandes_fontes": k, "brandes_segundos": t_bt,
                  "spearman_brandes_amostrado_vs_exato": float(stats.spearmanr(bt_amostra, bt_exato_ig)[0]),
                  "sobreposicao_top100_amostrado_vs_exato": len(a & b) / top}}
    ctx.cache["central"] = c
    return c


def robustez(ctx: Contexto) -> dict:
    G = ctx.G
    r = R.r_remocao(G.n)
    reps = 5 if ctx.rapido else 30
    cen = centralidades(ctx)
    base = R.metricas_apos_remocao(G, [], seed=0)
    aleat = R.remocao_aleatoria(G, r, reps=reps, seed=SEED)
    ataques = {nome: R.remocao_ranking(G, cen[chave], r, seed=0)
               for nome, chave in (("Grau", "grau"), ("Betweenness", "betweenness"), ("PageRank", "pagerank"))}
    met = ["S", "c", "L", "E"]
    linhas = [{"Estratégia": "Sem remoção", **{m: base[m] for m in met}}]
    resumo_aleat = {}
    for m in met:
        resumo_aleat[m] = intervalo_confianca([x[m] for x in aleat])
    linhas.append({"Estratégia": f"Aleatória (média de {reps})", **{m: resumo_aleat[m]["media"] for m in met}})
    for nome, x in ataques.items():
        linhas.append({"Estratégia": f"Ataque: {nome}", **{m: x[m] for m in met}})
    df = pd.DataFrame(linhas).rename(columns={"S": "S (fração na LCC)", "c": "c (componentes)",
                                              "L": "L (caminho médio)", "E": "E (eficiência)"})
    salvar_tabela(df, "parte3_robustez", f"Robustez: remoção de r = ⌈0,05n⌉ = {r} vértices.".replace("⌈", "").replace("⌉", ""),
                  "tab:robustez", casas={"S (fração na LCC)": 4, "c (componentes)": 1, "L (caminho médio)": 3,
                                         "E (eficiência)": 4}, alinhamento="lrrrr")
    ica = pd.DataFrame([{"Métrica": m, "Média": v["media"], "DP": v["dp"], "IC95 inf.": v["ic_inf"],
                         "IC95 sup.": v["ic_sup"], "Método": v["metodo"]} for m, v in resumo_aleat.items()])
    salvar_tabela(ica, "parte3_robustez_aleatoria_ic", f"Remoção aleatória: {reps} repetições, IC de 95%.",
                  "tab:robustezic", casas=4, alinhamento="lrrrrl")

    fig, axs = plt.subplots(1, 4, figsize=(14, 3.6))
    tit = {"S": "A) fração na LCC (S)", "c": "B) nº de componentes (c)", "L": "C) caminho médio na LCC (L)",
           "E": "D) eficiência global (E)"}
    mk = {"Grau": ("v", CORES["vermelho"]), "Betweenness": ("s", CORES["roxo"]), "PageRank": ("D", CORES["laranja"])}
    for ax, m in zip(axs, met):
        ax.boxplot([x[m] for x in aleat], widths=0.5, showmeans=True)
        ax.axhline(base[m], color=CORES["cinza"], ls=":", label="intacto")
        for i, (nome, x) in enumerate(ataques.items()):
            ax.plot(1.35 + 0.1 * i, x[m], mk[nome][0], color=mk[nome][1], ms=7, label=f"ataque {nome}")
        ax.set_xticks([1]); ax.set_xticklabels([f"aleatória\n({reps} rep.)"])
        ax.set_title(tit[m], fontsize=9)
    axs[0].legend(fontsize=7, loc="lower left")
    salvar(fig, "p3_robustez_5pct")

    fracs = np.round(np.arange(0, 0.51, 0.025), 3)
    curvas = {"Aleatória": R.curva(G, None, fracs, reps=3 if not ctx.rapido else 1, seed=SEED)}
    for nome, chave in (("Grau", "grau"), ("Betweenness", "betweenness"), ("PageRank", "pagerank")):
        curvas[nome] = R.curva(G, np.argsort(-cen[chave], kind="stable"), fracs)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    cores = {"Aleatória": CORES["azul"], "Grau": CORES["vermelho"], "Betweenness": CORES["roxo"], "PageRank": CORES["laranja"]}
    for nome, cv in curvas.items():
        ax[0].plot([x["frac"] for x in cv], [x["S"] for x in cv], "-o", ms=3, color=cores[nome], label=nome)
        ax[1].plot([x["frac"] for x in cv], [x["c"] for x in cv], "-o", ms=3, color=cores[nome], label=nome)
    for a in ax:
        a.axvline(0.05, color="k", ls=":", lw=0.8)
    ax[0].set(xlabel="fração removida", ylabel="|LCC| / n original", title="Tamanho da maior componente")
    ax[1].set(xlabel="fração removida", ylabel="nº de componentes", yscale="log", title="Fragmentação")
    ax[0].legend(fontsize=8)
    salvar(fig, "p3_robustez_curvas")
    pd.DataFrame([{"estrategia": n, **x} for n, cv in curvas.items() for x in cv]).to_csv(
        RES / "parte3_robustez_curvas.csv", index=False)

    out = {"r": r, "n": G.n, "base": base, "aleatoria_ic": resumo_aleat, "aleatoria_reps": aleat,
           "ataques": ataques, "centralidades": cen["info"],
           "razao_S_ataque_grau_sobre_aleatoria": ataques["Grau"]["S"] / resumo_aleat["S"]["media"],
           "curvas": curvas}
    salvar_json(out, "parte3_robustez")
    return out


def descoberta(ctx: Contexto) -> dict:
    G, nos = ctx.G, ctx.nos
    cen = centralidades(ctx)
    cl = ctx.cache.get("clus") or M.clusterizacao(G.adj)
    graus = G.graus()
    loc = cl["local"]
    top = np.argsort(-graus, kind="stable")[:10]
    hubs = pd.DataFrame({"id Orkut": [str(x) for x in G.ids[top]], "grau (subgrafo)": graus[top],
                         "grau (Orkut completo)": nos["grau_original"].to_numpy()[top],
                         "betweenness": cen["betweenness"][top], "PageRank": cen["pagerank"][top],
                         "clusterização local": loc[top], "nº comunidades": nos["n_comunidades"].to_numpy()[top]})
    salvar_tabela(hubs, "parte3_hubs", "Os 10 vértices de maior grau no subgrafo.", "tab:hubs",
                  casas={"id Orkut": 0, "grau (subgrafo)": 0, "grau (Orkut completo)": 0, "betweenness": 4,
                         "PageRank": 5, "clusterização local": 3, "nº comunidades": 0})
    # estrutura comunitária: fração de arestas dentro de alguma comunidade comum
    coms = [set(map(int, str(s).split(";"))) for s in nos["comunidades"]]
    intra = np.array([bool(coms[a] & coms[b]) for a, b in zip(G.u.tolist(), G.v.tolist())])
    ncom = nos["n_comunidades"].to_numpy()
    sp = lambda x, y: float(stats.spearmanr(x, y)[0])
    # concentração: fração das arestas que tocam o top 1% / 5% por grau
    ordem = np.argsort(-graus)
    conc = {}
    for f in (0.01, 0.05):
        sel = np.zeros(G.n, bool); sel[ordem[:math.ceil(f * G.n)]] = True
        conc[f"arestas_tocando_top{int(f*100)}pct"] = float((sel[G.u] | sel[G.v]).mean())
    # clusterização × grau (hierarquia: C(k) ~ k^-1 em redes hierárquicas)
    ks = np.unique(graus)
    ck = np.array([loc[graus == k].mean() for k in ks])
    okk = (ks >= 5) & (ck > 0)
    reg = stats.linregress(np.log10(ks[okk]), np.log10(ck[okk]))
    r = R.r_remocao(G.n)
    top_g = np.argsort(-graus, kind="stable")[:r]
    top_b = np.argsort(-cen["betweenness"], kind="stable")[:r]
    out = {"hubs": hubs.to_dict("records"),
           "r": r, "sobreposicao_top_r_grau_betweenness": len(set(top_g) & set(top_b)),
           "frac_multicom_top_r_grau": float((ncom[top_g] > 1).mean()),
           "frac_multicom_top_r_betweenness": float((ncom[top_b] > 1).mean()),
           "clusterizacao_media_top_r_grau": float(loc[top_g].mean()),
           "clusterizacao_media_top_r_betweenness": float(loc[top_b].mean()),
           "razao_betweenness_hubs_top10_sobre_media": float(cen["betweenness"][top].mean() / cen["betweenness"].mean()),
           "frac_arestas_intra_comunidade": float(intra.mean()),
           "frac_vertices_em_mais_de_1_comunidade": float((ncom > 1).mean()),
           "spearman_grau_vs_clusterizacao": sp(graus, loc),
           "spearman_ncomunidades_vs_betweenness": sp(ncom, cen["betweenness"]),
           "spearman_grau_vs_betweenness": sp(graus, cen["betweenness"]),
           "betweenness_media_1_comunidade": float(cen["betweenness"][ncom == 1].mean()),
           "betweenness_media_2mais_comunidades": float(cen["betweenness"][ncom >= 2].mean()),
           "C_de_k_expoente": float(-reg.slope), "C_de_k_r2": float(reg.rvalue ** 2), **conc}
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ax[0].loglog(ks, ck, "o", ms=3, color=CORES["azul"])
    xs = np.geomspace(5, ks.max(), 20)
    ax[0].loglog(xs, 10 ** reg.intercept * xs ** reg.slope, "--", color=CORES["vermelho"],
                 label=f"C(k) ~ k^-{-reg.slope:.2f} (R²={reg.rvalue**2:.2f})")
    ax[0].set(xlabel="grau k", ylabel="clusterização média C(k)", title="Clusterização por grau")
    ax[0].legend(fontsize=8)
    grupos = [cen["betweenness"][ncom == i] for i in range(1, min(ncom.max(), 4) + 1)]
    ax[1].boxplot(grupos, showfliers=False)
    ax[1].set_xticklabels([str(i) if i < 4 else "4+" for i in range(1, len(grupos) + 1)])
    ax[1].set(xlabel="nº de comunidades do vértice", ylabel="betweenness", title="Vértices de fronteira entre comunidades")
    salvar(fig, "p3_descoberta")
    salvar_json(out, "parte3_descoberta")
    return out


def modelos_classicos(ctx: Contexto) -> dict:
    G = ctx.G
    n, m = G.n, G.m
    kmed = 2 * m / n
    m0 = round(m / n)
    kws = int(round(kmed / 2) * 2)
    kf = 300 if ctx.rapido else 500
    dist = ctx.cache.get("dist") or todos_para_todos(G, g=ctx.ig)
    cl = ctx.cache.get("clus") or M.clusterizacao(G.adj)
    grafos = {"Orkut G27": G, "Erdős–Rényi": Mo.erdos_renyi_gnm(n, m, seed=SEED),
              f"Barabási–Albert (m={m0})": Mo.barabasi_albert(n, m0, seed=SEED),
              f"Watts–Strogatz (k={kws}, p=0,1)": Mo.watts_strogatz(n, kws, 0.1, seed=SEED)}
    linhas, dists = [], {}
    for nome, H in grafos.items():
        g = H.graus()
        if H is G:
            C, L = cl["clusterizacao_media"], dist["caminho_medio"]
        else:
            C, L = M.clusterizacao(H.adj)["clusterizacao_media"], caminho_medio_amostrado(H, k=kf, seed=1)
        kk, cc = LP._tabela(graus=g)
        pl = LP.escolher_kmin(kk, cc)
        vu = LP.comparar_exponencial(kk, cc, pl["kmin"], pl["gamma"])
        cauda = "potência" if vu["favorece"] == "lei de potência" and vu["significativo"] else (
            "exponencial" if vu["significativo"] else "inconclusivo")
        linhas.append({"Modelo": nome, "n": H.n, "m": H.m, "grau médio": g.mean(), "grau máx.": int(g.max()),
                       "C": C, "L": L, "γ MLE": pl["gamma"], "k_min": pl["kmin"], "cauda (Vuong)": cauda})
        dists[nome] = g
    df = pd.DataFrame(linhas)
    salvar_tabela(df, "parte3_modelos", "Comparação com modelos clássicos (mesmos n e m aproximados).",
                  "tab:modelos", casas={"n": 0, "m": 0, "grau médio": 2, "grau máx.": 0, "C": 4, "L": 3,
                                        "γ MLE": 2, "k_min": 0}, alinhamento="lrrrrrrrrl")
    fig, ax = plt.subplots(figsize=(6, 4.2))
    for (nome, g), cor in zip(dists.items(), [CORES["azul"], CORES["cinza"], CORES["vermelho"], CORES["verde"]]):
        k, pk, cc = M.distribuicao_graus(g)
        ax.loglog(k, cc, "-", color=cor, lw=2 if nome == "Orkut G27" else 1.2, label=nome)
    ax.set(xlabel="grau k", ylabel="P(K ≥ k)", title="CCDF de graus: Orkut × modelos")
    ax.legend(fontsize=8)
    salvar(fig, "p3_modelos_ccdf")
    out = {"tabela": linhas}
    salvar_json(out, "parte3_modelos")
    return out


# ============================================================== CONSOLIDAÇÃO
def consolidar() -> dict:
    etapas = ["parte1", "parte2", "parte3_small_world", "parte3_lei_potencia", "parte3_robustez",
              "parte3_descoberta", "parte3_modelos"]
    tudo = {"tratamento": json.loads((PROC / "resumo_tratamento.json").read_text(encoding="utf-8"))}
    for e in etapas:
        try:
            tudo[e] = ler_json(e)
        except FileNotFoundError:
            pass
    # versão enxuta (sem listas por repetição)
    enx = json.loads(json.dumps(tudo))
    enx.get("parte3_robustez", {}).pop("aleatoria_reps", None)
    enx.get("parte3_robustez", {}).pop("curvas", None)
    enx.get("parte2", {}).pop("escala", None)
    (RES / "resultados.json").write_text(json.dumps(enx, indent=2, ensure_ascii=False), encoding="utf-8")

    p1, p2 = tudo["parte1"], tudo["parte2"]
    sw, pl = tudo["parte3_small_world"], tudo["parte3_lei_potencia"]
    rb, ds, mo = tudo["parte3_robustez"], tudo["parte3_descoberta"], tudo["parte3_modelos"]
    tr = tudo["tratamento"]
    er = sw["Erdős–Rényi G(n,m)"]
    cf = sw["Configuração (trocas duplas, graus preservados)"]
    sub = pl["Subgrafo G27 (LCC, 9.917 vértices)"]
    full = pl["Orkut completo (3.072.441 vértices)"]
    tempos = {t["Algoritmo"]: t for t in p2["tempos_grafo_completo"]}
    fw = p2["floyd_warshall"]
    v = {
        "nV": p1["vertices"], "nE": p1["arestas"], "grauMin": p1["grau_min"], "grauMax": p1["grau_max"],
        "grauMedio": p1["grau_medio"], "grauMediano": (p1["grau_mediano"], 1), "densidade": (p1["densidade"], 5),
        "nComp": p1["n_componentes"], "diametro": p1["diametro"], "raio": p1["raio"],
        "caminhoMedio": (p1["caminho_medio"], 3), "triangulos": p1["triangulos"],
        "clusterMedia": (p1["clusterizacao_media"], 4), "transitividade": (p1["transitividade"], 4),
        "origV": tr["auditoria_original"]["vertices_distintos"], "origE": tr["auditoria_original"]["linhas_de_aresta"],
        "origLacos": tr["auditoria_original"]["auto_lacos"],
        "origDuplicadas": tr["auditoria_original"]["linhas_duplicadas_(u,v)_ou_(v,u)"],
        "nComunidadesSel": tr["selecao"]["k_comunidades"], "sementeComunidade": tr["selecao"]["comunidade_semente"],
        "eulerImpares": p2["execucao"]["euler"]["n_impares"],
        "nPontes": p2["execucao"]["tarjan"]["pontes"], "nArticulacoes": p2["execucao"]["tarjan"]["pontos_articulacao"],
        "bfRodadas": p2["execucao"]["bellman_ford"]["rodadas_ate_convergir"],
        "agmPeso": (p2["execucao"]["agm"]["peso_kruskal"], 3),
        "tBFS": (tempos["BFS"]["Média (ms)"], 2), "tDFS": (tempos["DFS"]["Média (ms)"], 2),
        "tDijkstra": (tempos["Dijkstra"]["Média (ms)"], 2), "tBellmanFord": (tempos["Bellman-Ford"]["Média (ms)"], 2),
        "tTarjan": (tempos["Tarjan (pontes/articulações)"]["Média (ms)"], 2),
        "tPrim": (tempos["Prim"]["Média (ms)"], 2), "tKruskal": (tempos["Kruskal"]["Média (ms)"], 2),
        "fwMemoriaMB": (fw["memoria_matriz_grafo_completo_MB"], 0),
        "fwEstimativaNumpyH": (fw["ajustes"]["Floyd-Warshall (numpy)"]["estimativa_grafo_completo_s"] / 3600, 1),
        "fwEstimativaPythonDias": (fw["ajustes"]["Floyd-Warshall (Python)"]["estimativa_grafo_completo_s"] / 86400, 2),
        "swCrandER": (er["C_rand"], 5), "swLrandER": (er["L_rand"], 3), "swSigmaER": (er["sigma"], 1),
        "swOmegaER": (er["omega"], 3), "swCrandConf": (cf["C_rand"], 4), "swLrandConf": (cf["L_rand"], 3),
        "swSigmaConf": (cf["sigma"], 2), "swOmegaConf": (cf["omega"], 3),
        "plGammaSub": (sub["mle"]["gamma"], 2), "plKminSub": (sub["mle"]["kmin"], 0),
        "plGammaFull": (full["mle"]["gamma"], 2), "plKminFull": (full["mle"]["kmin"], 0),
        "plCaudaFull": full["mle"]["n_cauda"],
        "robR": rb["r"], "robSaleat": (rb["aleatoria_ic"]["S"]["media"], 4),
        "robSgrau": (rb["ataques"]["Grau"]["S"], 4), "robSbet": (rb["ataques"]["Betweenness"]["S"], 4),
        "robSpr": (rb["ataques"]["PageRank"]["S"], 4), "robCaleat": (rb["aleatoria_ic"]["c"]["media"], 1),
        "robCgrau": rb["ataques"]["Grau"]["c"], "robCbet": rb["ataques"]["Betweenness"]["c"],
        "robLaleat": (rb["aleatoria_ic"]["L"]["media"], 3), "robLgrau": (rb["ataques"]["Grau"]["L"], 3),
        "robLbet": (rb["ataques"]["Betweenness"]["L"], 3),
        "fracIntraComunidade": (100 * ds["frac_arestas_intra_comunidade"], 1),
        "robForaAleat": (np.mean([x["n_restante"] - x["maior"] for x in rb["aleatoria_reps"]]), 1),
        "robForaGrau": rb["ataques"]["Grau"]["n_restante"] - rb["ataques"]["Grau"]["maior"],
        "robForaBet": rb["ataques"]["Betweenness"]["n_restante"] - rb["ataques"]["Betweenness"]["maior"],
        "robForaPr": rb["ataques"]["PageRank"]["n_restante"] - rb["ataques"]["PageRank"]["maior"],
        "sobreposicaoTopR": ds["sobreposicao_top_r_grau_betweenness"],
        "sobreposicaoTopRPct": (100 * ds["sobreposicao_top_r_grau_betweenness"] / ds["r"], 0),
        "multicomTopGrau": (100 * ds["frac_multicom_top_r_grau"], 0),
        "multicomTopBet": (100 * ds["frac_multicom_top_r_betweenness"], 0),
        "clusTopGrau": (ds["clusterizacao_media_top_r_grau"], 3),
        "clusTopBet": (ds["clusterizacao_media_top_r_betweenness"], 3),
        "betHubsRazao": (ds["razao_betweenness_hubs_top10_sobre_media"], 0),
        "ckExpoente": (ds["C_de_k_expoente"], 2),
    }
    escrever_macros(v)

    md = ["# Resultados consolidados — G27 com-Orkut", "",
          f"Máquina: {p2['maquina']['cpu']} ({p2['maquina']['nucleos_logicos']} núcleos lógicos), "
          f"Python {p2['maquina']['python']}, {p2['maquina']['so']}.", "",
          "## Parte I — estrutura", "",
          f"- |V| = {p1['vertices']:,}, |E| = {p1['arestas']:,}; grau mín./máx./médio/mediano = "
          f"{p1['grau_min']}/{p1['grau_max']}/{p1['grau_medio']:.2f}/{p1['grau_mediano']:.0f}; densidade = {p1['densidade']:.5f}",
          f"- {p1['n_componentes']} componente(s); diâmetro = {p1['diametro']}, raio = {p1['raio']}, "
          f"caminho médio = {p1['caminho_medio']:.3f} (exato, {p1['pares_ordenados']:,} pares ordenados)",
          f"- triângulos = {p1['triangulos']:,}; clusterização média = {p1['clusterizacao_media']:.4f}; "
          f"transitividade = {p1['transitividade']:.4f}", "",
          "## Parte II — algoritmos", "",
          f"- Euler: {p2['execucao']['euler']['n_impares']:,} vértices de grau ímpar → não euleriano nem semi-euleriano",
          f"- Tarjan: {p2['execucao']['tarjan']['pontes']} pontes e {p2['execucao']['tarjan']['pontos_articulacao']} pontos de articulação "
          f"(confere com networkx: {p2['execucao']['tarjan']['confere_networkx']})",
          f"- Dijkstra com peso unitário = BFS: {p2['execucao']['dijkstra_unitario_igual_bfs']}; "
          f"Bellman-Ford = Dijkstra: {p2['execucao']['bellman_ford']['igual_dijkstra']} "
          f"({p2['execucao']['bellman_ford']['rodadas_ate_convergir']} rodadas de {p2['execucao']['bellman_ford']['limite_teorico_n_menos_1']:,} possíveis)",
          f"- AGM: Prim = Kruskal = {p2['execucao']['agm']['peso_kruskal']:.4f} (igraph confere: {p2['execucao']['agm']['confere_igraph']})",
          f"- Floyd-Warshall: validado (numpy = Python = Dijkstra) em n = {fw['validacao']['n']}; no grafo completo exigiria "
          f"{fw['memoria_matriz_grafo_completo_MB']:.0f} MB e ~{fw['ajustes']['Floyd-Warshall (numpy)']['estimativa_grafo_completo_s']/3600:.1f} h (numpy) "
          f"ou ~{fw['ajustes']['Floyd-Warshall (Python)']['estimativa_grafo_completo_s']/86400:.2f} dias (Python puro)", "",
          "| Algoritmo | Complexidade | média (ms) | DP (ms) | IC95 (ms) | método |", "|---|---|---|---|---|---|"]
    md += [f"| {t['Algoritmo']} | {t['Complexidade']} | {t['Média (ms)']:.3f} | {t['DP (ms)']:.3f} | "
           f"[{t['IC95 inf. (ms)']:.3f}; {t['IC95 sup. (ms)']:.3f}] | {t['Método do IC']} |" for t in p2["tempos_grafo_completo"]]
    md += ["", "## Parte III — redes complexas", "",
           f"- Small-world vs ER: C/C_rand = {er['gamma_C_sobre_Crand']:.1f}, L/L_rand = {er['lambda_L_sobre_Lrand']:.3f}, "
           f"σ = {er['sigma']:.1f}, ω = {er['omega']:.3f} → small-world: {er['small_world']}",
           f"- Small-world vs configuração (graus preservados): C/C_conf = {cf['gamma_C_sobre_Crand']:.2f}, "
           f"L/L_conf = {cf['lambda_L_sobre_Lrand']:.3f}, σ = {cf['sigma']:.2f}",
           f"- Lei de potência (MLE): subgrafo γ = {sub['mle']['gamma']:.2f} (k ≥ {sub['mle']['kmin']:.0f}, cauda {sub['mle']['n_cauda']}), "
           f"vs exponencial: {sub['vs_exponencial']['favorece']} (p = {sub['vs_exponencial']['p_valor']:.3g}); "
           f"Orkut completo γ = {full['mle']['gamma']:.2f} (k ≥ {full['mle']['kmin']:.0f}, cauda {full['mle']['n_cauda']:,}), "
           f"vs exponencial: {full['vs_exponencial']['favorece']} (p = {full['vs_exponencial']['p_valor']:.3g})",
           f"- Robustez (r = {rb['r']}): S aleatória = {rb['aleatoria_ic']['S']['media']:.4f} "
           f"[{rb['aleatoria_ic']['S']['ic_inf']:.4f}; {rb['aleatoria_ic']['S']['ic_sup']:.4f}]; "
           + "; ".join(f"ataque {k}: S = {x['S']:.4f}, c = {x['c']}, L = {x['L']:.3f}" for k, x in rb["ataques"].items()),
           f"- Descoberta: {100*ds['frac_arestas_intra_comunidade']:.1f}% das arestas ficam dentro de alguma comunidade comum; "
           f"C(k) ~ k^-{ds['C_de_k_expoente']:.2f}; Spearman(nº comunidades, betweenness) = {ds['spearman_ncomunidades_vs_betweenness']:.3f}", "",
           "### Modelos clássicos", "", "| Modelo | m | grau máx. | C | L | γ MLE | cauda (Vuong) |", "|---|---|---|---|---|---|---|"]
    md += [f"| {r['Modelo']} | {r['m']:,} | {r['grau máx.']} | {r['C']:.4f} | {r['L']:.3f} | {r['γ MLE']:.2f} | {r['cauda (Vuong)']} |" for r in mo["tabela"]]
    ra, at = rb["aleatoria_ic"], rb["ataques"]
    fora = {k: x["n_restante"] - x["maior"] for k, x in at.items()}
    fora_al = float(np.mean([x["n_restante"] - x["maior"] for x in rb["aleatoria_reps"]]))
    md += ["", "## Descoberta principal", "",
           "**Os hubs não são os pontos críticos desta região do Orkut: eles ficam no núcleo redundante onde as comunidades "
           "se sobrepõem, e quem segura a periferia são outros vértices.** "
           f"(i) A clusterização é {cf['gamma_C_sobre_Crand']:.0f}× a de um aleatório com a MESMA sequência de graus e "
           f"{100*ds['frac_arestas_intra_comunidade']:.1f}% das arestas ligam membros de uma mesma comunidade: o agrupamento vem das comunidades. "
           f"(ii) {100*ds['frac_multicom_top_r_grau']:.0f}% dos {ds['r']} vértices de maior grau pertencem a 2+ comunidades "
           f"(média geral: {100*ds['frac_vertices_em_mais_de_1_comunidade']:.0f}%), ou seja, os hubs estão na zona de sobreposição, que tem muitos caminhos alternativos; "
           f"removê-los desprende só {fora['Grau']} vértices da LCC (aleatória: {fora_al:.1f}). "
           f"(iii) Só {ds['sobreposicao_top_r_grau_betweenness']} dos {ds['r']} vértices de maior betweenness estão também no top de grau; "
           f"removê-los desprende {fora['Betweenness']} vértices ({at['Betweenness']['c']} componentes) e alonga L em "
           f"{100*(at['Betweenness']['L']/rb['base']['L']-1):.0f}% (grau: {100*(at['Grau']['L']/rb['base']['L']-1):.0f}%). "
           f"Mesmo assim a rede segue muito robusta (S ≥ {min(x['S'] for x in at.values()):.3f}). "
           f"(iv) A lei de potência existe no Orkut completo (γ = {full['mle']['gamma']:.2f}, k ≥ {full['mle']['kmin']:.0f}) mas some no subgrafo "
           f"(cauda exponencial): a amostragem por comunidades preserva C e L, mas corta a cauda de graus."]
    (RES / "resultados.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    return enx


def executar_tudo(rapido: bool = False):
    ctx = preparar(rapido)
    etapas = [("Parte I", parte1), ("Parte II", parte2), ("Small-world", small_world),
              ("Lei de potência", lei_potencia), ("Robustez", robustez), ("Descoberta", descoberta),
              ("Modelos clássicos", modelos_classicos)]
    for nome, f in etapas:
        t0 = time.time()
        f(ctx)
        print(f"[ok] {nome} ({time.time() - t0:.0f} s)", flush=True)
    consolidar()
    print("[ok] consolidação -> resultados/resultados.json, resultados.md, valores.tex")
    return ctx
