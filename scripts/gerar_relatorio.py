#!/usr/bin/env python3
"""Monta a pasta relatorio/ a partir de resultados/ (nenhum número digitado à mão).

  * relatorio/valores.tex        <- cópia de resultados/valores.tex
  * relatorio/valores_extra.tex  <- macros adicionais derivadas dos JSON/CSV
  * relatorio/tabelas/*.tex      <- tabelas no formato do artigo (booktabs)
  * relatorio/figuras/*.pdf      <- figuras vetoriais

Uso:  python scripts/gerar_relatorio.py [--compilar]
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ))
from src.relatorio import esc, fmt_num, valor_macro  # noqa: E402

RES = RAIZ / "resultados"
REL = RAIZ / "relatorio"
TAB = REL / "tabelas"
FIGD = REL / "figuras"
MATRICULA = "224215935"


def macro(nome, v, casas=2):
    return rf"\newcommand{{\{nome}}}{{{valor_macro(v, casas)}}}"


def tabela(nome, df, caption, label, casas=None, larga=False, alinh=None, nota=None):
    casas = casas or {}
    cols = list(df.columns)
    alinh = alinh or ("l" + "r" * (len(cols) - 1))
    env = "table*" if larga else "table"
    L = [rf"\begin{{{env}}}[!htb]", r"\centering", rf"\caption{{{caption}}}", rf"\label{{{label}}}",
         r"\begin{adjustbox}{max width=\linewidth}", rf"\begin{{tabular}}{{{alinh}}}", r"\toprule",
         " & ".join(esc(c) for c in cols) + r" \\", r"\midrule"]
    for _, row in df.iterrows():
        cel = []
        for c in cols:
            v = row[c]
            cel.append(esc(v) if isinstance(v, str) else esc(fmt_num(v, casas.get(c, 2))))
        L.append(" & ".join(cel) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}", r"\end{adjustbox}"]
    if nota:
        L.append(rf"\par\smallskip{{\footnotesize {nota}\par}}")
    L.append(rf"\end{{{env}}}")
    (TAB / f"{nome}.tex").write_text("\n".join(L) + "\n", encoding="utf-8")


def main(compilar: bool):
    TAB.mkdir(parents=True, exist_ok=True)
    FIGD.mkdir(parents=True, exist_ok=True)
    shutil.copy(RES / "valores.tex", REL / "valores.tex")
    for f in (RAIZ / "figuras").glob("*.pdf"):
        shutil.copy(f, FIGD / f.name)
    R = json.loads((RES / "resultados.json").read_text(encoding="utf-8"))
    p1, p2 = R["parte1"], R["parte2"]
    sw, pl, rb = R["parte3_small_world"], R["parte3_lei_potencia"], R["parte3_robustez"]
    ds, mo, tr = R["parte3_descoberta"], R["parte3_modelos"], R["tratamento"]
    er, cf = sw["Erdős–Rényi G(n,m)"], sw["Configuração (trocas duplas, graus preservados)"]
    sub, full = pl["Subgrafo G27 (LCC, 9.917 vértices)"], pl["Orkut completo (3.072.441 vértices)"]
    ex, fw = p2["execucao"], p2["floyd_warshall"]
    aj = {a["Algoritmo"]: a for a in p2["ajustes"]}
    tempos = {t["Algoritmo"]: t for t in p2["tempos_grafo_completo"]}
    ra = rb["aleatoria_ic"]
    hist = p1["hist_distancias"]
    pares = p1["pares_ordenados"]

    M = []
    add = lambda n, v, c=2: M.append(macro(n, v, c))
    # máquina
    add("cpu", p2["maquina"]["cpu"]); add("nucleos", p2["maquina"]["nucleos_logicos"])
    add("pyver", p2["maquina"]["python"])
    # original / tratamento
    a = tr["auditoria_original"]
    add("origGrauMax", a["grau_maximo"]); add("origGrauMedio", a["grau_medio"]); add("origGrauMediano", a["grau_mediano"], 0)
    add("origIdMax", a["id_maximo"]); add("origSegLeitura", a["segundos_leitura"], 0)
    add("indV", tr["induzido"]["vertices"]); add("indE", tr["induzido"]["arestas"])
    add("grauPreservado", 100 * pl["fracao_media_grau_preservado"], 0)
    add("spearmanGrauOrig", pl["correlacao_spearman_grau_original_vs_sub"], 3)
    # Parte I
    add("desvioGrau", p1["grau_desvio"]); add("centroTam", p1["centro_tamanho"]); add("periferiaTam", p1["periferia_tamanho"])
    add("paresOrdenados", pares)
    add("fracDistTres", 100 * hist[3] / pares, 1); add("fracDistQuatro", 100 * hist[4] / pares, 1)
    add("fracAteQuatro", 100 * sum(hist[1:5]) / pares, 1)
    add("segDistancias", p1["segundos_distancias"], 1); add("segTriangulos", p1["segundos_triangulos"], 2)
    # pesos
    add("wMin", p2["pesos"]["w_min"], 4); add("wMediana", p2["pesos"]["w_mediana"], 4)
    add("wFracUm", 100 * p2["pesos"]["frac_w_igual_1"], 1)
    # Parte II
    add("bfsExc", ex["bfs"]["excentricidade_fonte"]); add("bfsNivelQuatro", ex["bfs"]["vertices_por_nivel"][4])
    add("dfsProf", ex["dfs"]["profundidade_max_arvore"])
    add("dijDistMax", ex["dijkstra_ponderado"]["dist_max"], 3); add("dijSaltosMedios", ex["dijkstra_ponderado"]["saltos_medios_no_caminho_ponderado"], 2)
    add("dijSaltosMax", ex["dijkstra_ponderado"]["saltos_max"])
    add("bfNegN", ex["bellman_ford_demo_ciclo_negativo"]["n"]); add("bfNegRodadas", ex["bellman_ford_demo_ciclo_negativo"]["rodadas"])
    add("pontesFolha", ex["tarjan"]["pontes_com_folha"])
    add("agmFracUm", 100 * ex["agm"]["frac_arestas_peso_1_na_agm"], 1)
    for nome, chave in (("BFS", "BFS"), ("DFS", "DFS"), ("Euler", "Euler (verificação)"), ("Dijkstra", "Dijkstra"),
                        ("BellmanFord", "Bellman-Ford"), ("Tarjan", "Tarjan (pontes/articulações)"), ("Prim", "Prim"),
                        ("Kruskal", "Kruskal")):
        t = tempos[chave]
        add(f"t{nome}Med", t["Média (ms)"]); add(f"t{nome}Dp", t["DP (ms)"])
        add(f"t{nome}Inf", t["IC95 inf. (ms)"]); add(f"t{nome}Sup", t["IC95 sup. (ms)"])
        add(f"rdois{nome}", aj[chave]["R² (c·f)"], 3)
    add("rdoisMin", min(a["R² (c·f)"] for a in p2["ajustes"]), 3)
    add("fwValN", fw["validacao"]["n"])
    fn, fp = fw["ajustes"]["Floyd-Warshall (numpy)"], fw["ajustes"]["Floyd-Warshall (Python)"]
    add("fwRdoisNumpy", fn["r2"], 4); add("fwRdoisPython", fp["r2"], 4)
    add("fwEstNumpyMin", fn["estimativa_grafo_completo_s"] / 60, 0); add("fwEstPythonH", fp["estimativa_grafo_completo_s"] / 3600, 1)
    med = {(m["Versão"], m["n"]): m["Média (ms)"] for m in fw["medicoes"]}
    add("fwRazao", med[("Floyd-Warshall (Python)", 200)] / med[("Floyd-Warshall (numpy)", 200)], 0)
    add("fwOps", fw["operacoes_grafo_completo"] / 1e12, 2)
    # Parte III
    add("swC", sw["C"], 4); add("swL", sw["L"], 3); add("swKmedio", sw["k_medio"])
    add("swGammaER", er["gamma_C_sobre_Crand"], 0); add("swLambdaER", er["lambda_L_sobre_Lrand"], 3)
    add("swGammaConf", cf["gamma_C_sobre_Crand"], 1); add("swLambdaConf", cf["lambda_L_sobre_Lrand"], 3)
    add("swClatt", er["C_latt"], 3); add("swCrandAnal", sw["analitico_er"]["C_rand_analitico"], 5)
    add("swLrandAnal", sw["analitico_er"]["L_rand_analitico"], 3); add("swFontes", sw["fontes_amostradas_L_nulos"])
    add("swTrocas", [i for i in sw["instancias"] if "trocas_realizadas" in i][0]["trocas_realizadas"])
    for tag, r in (("Sub", sub), ("Full", full)):
        add(f"plBruta{tag}", r["regressao_bruta"]["gamma"]); add(f"plBrutaRdois{tag}", r["regressao_bruta"]["r2"])
        add(f"plBins{tag}", r["regressao_bins_log"]["gamma"]); add(f"plBinsRdois{tag}", r["regressao_bins_log"]["r2"])
        add(f"plCcdf{tag}", r["regressao_ccdf"]["gamma"]); add(f"plCcdfRdois{tag}", r["regressao_ccdf"]["r2"])
        add(f"plKs{tag}", r["mle"]["ks"], 3); add(f"plErro{tag}", r["mle"]["erro_padrao"], 3)
        add(f"plCauda{tag}x", r["mle"]["n_cauda"]); add(f"plVuong{tag}", r["vs_exponencial"]["z_vuong"], 1)
        add(f"plKsPoisson{tag}", r["vs_poisson"]["ks_poisson"], 2)
    add("robCbetRazao", rb["ataques"]["Betweenness"]["c"] / ra["c"]["media"], 1)
    add("robLbase", rb["base"]["L"], 3); add("robEbase", rb["base"]["E"], 3)
    add("robEaleat", ra["E"]["media"], 3); add("robEgrau", rb["ataques"]["Grau"]["E"], 3)
    add("robEbet", rb["ataques"]["Betweenness"]["E"], 3); add("robEpr", rb["ataques"]["PageRank"]["E"], 3)
    add("robCpr", rb["ataques"]["PageRank"]["c"]); add("robLpr", rb["ataques"]["PageRank"]["L"], 3)
    add("robSaleatInf", ra["S"]["ic_inf"], 4); add("robSaleatSup", ra["S"]["ic_sup"], 4)
    add("robCaleatInf", ra["c"]["ic_inf"], 1); add("robCaleatSup", ra["c"]["ic_sup"], 1)
    add("robLaumBet", 100 * (rb["ataques"]["Betweenness"]["L"] / rb["base"]["L"] - 1), 0)
    add("robLaumGrau", 100 * (rb["ataques"]["Grau"]["L"] / rb["base"]["L"] - 1), 0)
    add("robLaumAleat", 100 * (ra["L"]["media"] / rb["base"]["L"] - 1), 1)
    c = rb["centralidades"]
    add("brandesFontes", c["brandes_fontes"]); add("brandesSeg", c["brandes_segundos"], 0)
    add("brandesSpearman", c["spearman_brandes_amostrado_vs_exato"], 3)
    add("brandesTopCem", 100 * c["sobreposicao_top100_amostrado_vs_exato"], 0)
    add("prIter", c["pagerank_iteracoes"])
    add("fracMultiCom", 100 * ds["frac_vertices_em_mais_de_1_comunidade"], 1)
    add("spGrauBet", ds["spearman_grau_vs_betweenness"], 2); add("spComBet", ds["spearman_ncomunidades_vs_betweenness"], 2)
    add("spGrauClus", ds["spearman_grau_vs_clusterizacao"], 2)
    add("betRazaoCom", ds["betweenness_media_2mais_comunidades"] / ds["betweenness_media_1_comunidade"], 1)
    add("ckRdois", ds["C_de_k_r2"], 2)
    add("topUmPct", 100 * ds["arestas_tocando_top1pct"], 1); add("topCincoPct", 100 * ds["arestas_tocando_top5pct"], 1)
    mt = {r["Modelo"].split(" (")[0]: r for r in mo["tabela"]}
    add("baGrauMax", mt["Barabási–Albert"]["grau máx."]); add("baGamma", mt["Barabási–Albert"]["γ MLE"])
    add("baC", mt["Barabási–Albert"]["C"], 4); add("wsC", mt["Watts–Strogatz"]["C"], 3); add("wsL", mt["Watts–Strogatz"]["L"], 3)
    add("erGrauMax", mt["Erdős–Rényi"]["grau máx."]); add("wsGrauMax", mt["Watts–Strogatz"]["grau máx."])
    add("razaoBFDij", tempos["Bellman-Ford"]["Média (ms)"] / tempos["Dijkstra"]["Média (ms)"], 1)
    esc_df = pd.read_csv(RES / "parte2_tempos_escala.csv").groupby("n")["m"].first()
    add("mnMin", (esc_df / esc_df.index).iloc[0], 1); add("mnMax", (esc_df / esc_df.index).iloc[-1], 1)
    hb = pd.read_csv(RES / "parte3_hubs.csv", dtype={"id Orkut": str})
    exc = hb.nlargest(2, "betweenness")
    resto = hb.drop(exc.index)["betweenness"].median()
    M.append(rf"\newcommand{{\hubExcA}}{{{exc.iloc[0]['id Orkut']}}}")
    M.append(rf"\newcommand{{\hubExcB}}{{{exc.iloc[1]['id Orkut']}}}")
    add("hubExcRazaoMin", exc["betweenness"].min() / resto, 0); add("hubExcRazaoMax", exc["betweenness"].max() / resto, 0)
    add("hubClusMin", hb.drop(exc.index)["clusterização local"].min()); add("hubClusMax", hb.drop(exc.index)["clusterização local"].max())
    add("pontesNaoFolha", ex["tarjan"]["pontes"] - ex["tarjan"]["pontes_com_folha"])
    (REL / "valores_extra.tex").write_text("% Gerado por scripts/gerar_relatorio.py\n" + "\n".join(M) + "\n", encoding="utf-8")

    # ------------------------------------------------------------------ tabelas
    t = pd.read_csv(RAIZ / "data" / "processed" / "tabela_tratamentos.csv")
    tabela("tratamentos", t, "Tratamentos aplicados ao com-Orkut, com as contagens após cada passo.", "tab:tratamentos",
           larga=True, alinh="p{3.2cm}p{8.6cm}p{5.2cm}")
    t = pd.read_csv(RES / "parte1_medidas.csv")
    t["Valor"] = [fmt_num(float(v), 4) if "." in str(v) else fmt_num(int(float(v))) for v in t["Valor"]]
    tabela("parte1", t, r"Medidas estruturais da LCC tratada (distâncias exatas sobre todos os pares).", "tab:parte1")
    t = pd.read_csv(RES / "parte1_distancias.csv")
    tabela("distancias", t, "Distribuição das distâncias (em saltos) entre os pares ordenados de vértices.", "tab:distancias",
           casas={"Distância": 0, "Pares ordenados": 0, "Fração (%)": 2})
    t = pd.read_csv(RES / "parte1_calculadas.csv")
    tabela("calculadas", t, "Medidas calculadas e não calculadas, com o método ou a justificativa.", "tab:calculadas",
           larga=True, alinh="p{4.6cm}cp{11cm}")
    t = pd.read_csv(RES / "parte2_tempos_grafo_completo.csv").drop(columns=["Repetições"])
    t["Método do IC"] = t["Método do IC"].str.replace("Normal (z)", "z", regex=False)
    tabela("tempos", t, rf"Tempo no grafo completo (n = {fmt_num(p1['vertices'])}, m = {fmt_num(p1['arestas'])}): "
           r"30 repetições após 1 aquecimento, IC de 95\% Normal (z = 1,96).", "tab:tempos", larga=True,
           alinh="llrrrrl", casas={"Média (ms)": 2, "DP (ms)": 2, "IC95 inf. (ms)": 2, "IC95 sup. (ms)": 2})
    e = pd.read_csv(RES / "parte2_tempos_escala.csv")
    pv = e.pivot_table(index="n", columns="Algoritmo", values="Média (ms)", sort=False).reset_index()
    mm = e.groupby("n")["m"].first().reset_index()
    pv = mm.merge(pv, on="n")
    pv = pv.rename(columns={"Tarjan (pontes/articulações)": "Tarjan", "Euler (verificação)": "Euler"})
    pv = pv[["n", "m", "BFS", "DFS", "Euler", "Dijkstra", "Bellman-Ford", "Tarjan", "Prim", "Kruskal"]]
    tabela("escala", pv, r"Tempo médio (ms) em subgrafos BFS crescentes; 10 repetições por tamanho. Os desvios-padrão e "
           r"os IC de 95\% (t-Student, 9 g.l.) estão em \texttt{resultados/parte2\_tempos\_escala.csv}.", "tab:escala",
           larga=True, casas={"n": 0, "m": 0})
    a = pd.read_csv(RES / "parte2_ajustes.csv")
    a = a[["Algoritmo", "f(n,m)", "c", "R² (c·f)", "R² (a·n+b·m)"]]
    tabela("ajustes", a, r"Ajuste de mínimos quadrados $t \approx c\cdot f(n,m)$ aos tempos médios da Tabela~\ref{tab:escala} "
           r"($t$ em ms) e ajuste alternativo $t \approx a n + b m$.", "tab:ajustes",
           casas={"c": 6, "R² (c·f)": 4, "R² (a·n+b·m)": 4}, alinh="llrrr")
    f = pd.read_csv(RES / "parte2_floyd_warshall.csv")
    f["Versão"] = f["Versão"].str.replace("Floyd-Warshall ", "", regex=False)
    f = f[["Versão", "n", "m", "Média (ms)", "DP (ms)", "IC95 inf. (ms)", "IC95 sup. (ms)"]]
    tabela("floyd", f, r"Floyd-Warshall em subgrafos BFS (pesos de \emph{embeddedness}); 10 repetições, IC de 95\% t-Student.",
           "tab:floyd", casas={"n": 0, "m": 0, "Média (ms)": 1, "DP (ms)": 1, "IC95 inf. (ms)": 1, "IC95 sup. (ms)": 1})
    s = pd.read_csv(RES / "parte3_small_world.csv")
    s["Referência"] = ["Erdős–Rényi G(n,m)", "Configuração (mesmos graus)"]
    s["Small-world?"] = ["sim" if str(x) == "True" else "não" for x in s["Small-world?"]]
    tabela("smallworld", s, r"Small-world: grafo observado contra as duas referências aleatórias (5 instâncias ER e 3 de "
           r"configuração; $L_{\text{ref}}$ estimado com 1.000 fontes por instância).", "tab:smallworld", larga=True,
           casas={"C": 4, "C_ref": 4, "C/C_ref": 1, "L": 3, "L_ref": 3, "L/L_ref": 3, "σ": 1, "ω": 3})
    linhas = []
    for nome, chave in (("Regressão log-log, pontos brutos: γ (R²)", "bruta"), ("Regressão log-log, bins logarítmicos: γ (R²)", "bins"),
                        ("Regressão na CCDF: γ = 1 − inclinação (R²)", "ccdf")):
        k = {"bruta": "regressao_bruta", "bins": "regressao_bins_log", "ccdf": "regressao_ccdf"}[chave]
        linhas.append({"Método": nome, "Subgrafo G27": f"{fmt_num(sub[k]['gamma'])} ({fmt_num(sub[k]['r2'])})",
                       "Orkut completo": f"{fmt_num(full[k]['gamma'])} ({fmt_num(full[k]['r2'])})"})
    linhas += [
        {"Método": "MLE discreta: γ ± erro-padrão", "Subgrafo G27": f"{fmt_num(sub['mle']['gamma'])} ± {fmt_num(sub['mle']['erro_padrao'])}",
         "Orkut completo": f"{fmt_num(full['mle']['gamma'])} ± {fmt_num(full['mle']['erro_padrao'], 3)}"},
        {"Método": "k_min (mínima distância KS) / vértices na cauda", "Subgrafo G27": f"{fmt_num(int(sub['mle']['kmin']))} / {fmt_num(sub['mle']['n_cauda'])}",
         "Orkut completo": f"{fmt_num(int(full['mle']['kmin']))} / {fmt_num(full['mle']['n_cauda'])}"},
        {"Método": "Distância KS na cauda", "Subgrafo G27": fmt_num(sub["mle"]["ks"], 3), "Orkut completo": fmt_num(full["mle"]["ks"], 3)},
        {"Método": "Vuong vs. exponencial: z (p)", "Subgrafo G27": f"{fmt_num(sub['vs_exponencial']['z_vuong'])} ({fmt_num(sub['vs_exponencial']['p_valor'], 4)})",
         "Orkut completo": f"{fmt_num(full['vs_exponencial']['z_vuong'], 1)} (< 1e-200)"},
        {"Método": "Modelo favorecido na cauda", "Subgrafo G27": sub["vs_exponencial"]["favorece"], "Orkut completo": full["vs_exponencial"]["favorece"]},
        {"Método": "KS global da Poisson (λ = grau médio)", "Subgrafo G27": fmt_num(sub["vs_poisson"]["ks_poisson"]),
         "Orkut completo": fmt_num(full["vs_poisson"]["ks_poisson"])},
    ]
    tabela("powerlaw", pd.DataFrame(linhas), r"Estimativas do expoente de $P(k) \sim C k^{-\gamma}$ no subgrafo e no Orkut completo.",
           "tab:powerlaw", larga=True, alinh="lcc")
    r = pd.read_csv(RES / "parte3_robustez.csv")
    tabela("robustez", r, rf"Robustez: remoção de $r = \lceil 0{{,}}05 n\rceil = {rb['r']}$ vértices. Aleatória: média de 30 "
           r"repetições; ataques: ranking fixo calculado no grafo intacto. $L$ e $E$ estimados com 200 fontes.", "tab:robustez",
           casas={"S (fração na LCC)": 4, "c (componentes)": 1, "L (caminho médio)": 3, "E (eficiência)": 4})
    r = pd.read_csv(RES / "parte3_robustez_aleatoria_ic.csv").drop(columns=["Método"])
    tabela("robustezic", r, r"Remoção aleatória: média, desvio-padrão e IC de 95\% (Normal, $n = 30$).", "tab:robustezic", casas=
           {"Média": 4, "DP": 4, "IC95 inf.": 4, "IC95 sup.": 4})
    h = pd.read_csv(RES / "parte3_hubs.csv", dtype={"id Orkut": str})
    tabela("hubs", h, "Os 10 vértices de maior grau no subgrafo (betweenness normalizada, estimada com 1.000 fontes).", "tab:hubs",
           larga=True, casas={"grau (subgrafo)": 0, "grau (Orkut completo)": 0, "betweenness": 4, "PageRank": 5,
                              "clusterização local": 3, "nº comunidades": 0})
    m = pd.read_csv(RES / "parte3_modelos.csv")
    tabela("modelos", m, r"Orkut G27 contra modelos clássicos gerados com o mesmo $n$ e $m$ aproximado. $L$ dos modelos "
           r"estimado com 500 fontes; γ e $k_{\min}$ pela MLE discreta; última coluna: teste de Vuong na cauda.", "tab:modelos",
           larga=True, casas={"n": 0, "m": 0, "grau médio": 2, "grau máx.": 0, "C": 4, "L": 3, "γ MLE": 2, "k_min": 0})
    print("relatorio/ atualizado:", len(list(TAB.glob('*.tex'))), "tabelas,", len(list(FIGD.glob('*.pdf'))), "figuras,",
          len(M), "macros extras")

    if compilar:
        r = subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", "main.tex"], cwd=REL,
                           capture_output=True, text=True, errors="replace")
        if r.returncode != 0:
            print(r.stdout[-3000:])
            sys.exit(r.returncode)
        shutil.copy(REL / "main.pdf", REL / f"{MATRICULA}.pdf")
        print(f"compilado -> relatorio/{MATRICULA}.pdf")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--compilar", action="store_true")
    main(ap.parse_args().compilar)
