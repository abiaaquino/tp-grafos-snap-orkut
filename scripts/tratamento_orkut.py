#!/usr/bin/env python3
"""
TP "Analisando Grafos do SNAP" — G27: com-Orkut (não direcionado).

Auditoria do grafo original + redução metodologicamente defensável.

Etapas (uma única leitura em streaming do .txt.gz, em blocos via pandas):
  1. Auditoria exata do arquivo bruto: linhas de comentário, arestas (linhas),
     vértices distintos, auto-laços, arestas duplicadas (inclusive u-v e v-u),
     grau máximo/médio. A contagem de duplicatas é exata e com pouca RAM:
     cada aresta vira uma chave canônica int64 (min<<32 | max), gravada em
     N_BUCKETS arquivos temporários por hash; depois cada bucket é ordenado
     sozinho com numpy (≈ 60 MB por vez).
  2. Seleção de comunidades ground-truth (com-orkut.top5000.cmty.txt.gz):
     semente aleatória FIXA escolhe uma comunidade inicial; depois, de forma
     gulosa, adiciona-se a comunidade com maior sobreposição relativa
     |C ∩ S| / |C| com o conjunto já selecionado (vizinhança coerente no
     espaço de comunidades). Coleta-se um superconjunto de vértices.
  3. Na mesma leitura, guardam-se as arestas induzidas por esse superconjunto.
  4. Avaliam-se os prefixos da lista de comunidades e escolhe-se o maior
     prefixo cuja MAIOR COMPONENTE CONEXA (LCC) do subgrafo induzido simples
     tenha no máximo ALVO_MAX vértices (e, idealmente, >= ALVO_MIN).
  5. Grafo final = LCC simples (sem laços, sem duplicatas), não direcionado.

Uso:
  python scripts/tratamento_orkut.py
  python scripts/tratamento_orkut.py --seed 42 --alvo-min 5000 --alvo-max 10000
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
import os
import platform
import random
import shutil
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx

RAIZ = Path(__file__).resolve().parents[1]
RAW = RAIZ / "data" / "raw"
PROC = RAIZ / "data" / "processed"

ARQ_GRAFO = RAW / "com-orkut.ungraph.txt.gz"
ARQ_CMTY = RAW / "com-orkut.top5000.cmty.txt.gz"

# Contagens oficiais da página do SNAP (https://snap.stanford.edu/data/com-Orkut.html)
ESPERADO_V = 3_072_441
ESPERADO_E = 117_185_083

N_BUCKETS = 32
CHUNK = 5_000_000

log = logging.getLogger("tratamento")


# ----------------------------------------------------------------------------
# Funções puras (testáveis com arrays pequenos)
# ----------------------------------------------------------------------------
def chave_canonica(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Chave int64 única por aresta não direcionada: (min << 32) | max."""
    lo = np.minimum(u, v).astype(np.int64)
    hi = np.maximum(u, v).astype(np.int64)
    return (lo << 32) | hi


def decodificar_chave(k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return (k >> 32).astype(np.int64), (k & 0xFFFFFFFF).astype(np.int64)


def ler_comunidades(caminho: Path) -> list[list[int]]:
    comunidades = []
    with gzip.open(caminho, "rt") as fh:
        for linha in fh:
            linha = linha.strip()
            if not linha or linha.startswith("#"):
                continue
            comunidades.append(sorted({int(x) for x in linha.split()}))
    return comunidades


def selecionar_comunidades(comunidades: list[list[int]], seed: int,
                           alvo_superconjunto: int, tam_max_cmty: int):
    """Expansão gulosa por sobreposição a partir de uma comunidade-semente.

    Retorna a lista ordenada de índices de comunidades e o conjunto união.
    Comunidades maiores que tam_max_cmty são ignoradas (evita que uma única
    comunidade estoure o alvo).
    """
    elegiveis = [i for i, c in enumerate(comunidades) if 2 <= len(c) <= tam_max_cmty]
    rng = random.Random(seed)
    semente = rng.choice(elegiveis)

    no_para_cmty = defaultdict(list)
    for i in elegiveis:
        for n in comunidades[i]:
            no_para_cmty[n].append(i)

    ordem = [semente]
    usados = {semente}
    S = set(comunidades[semente])
    while len(S) < alvo_superconjunto:
        # sobreposição |C ∩ S| para comunidades que tocam S
        inter = defaultdict(int)
        for n in S:
            for c in no_para_cmty.get(n, ()):
                if c not in usados:
                    inter[c] += 1
        candidatos = [(inter[c] / len(comunidades[c]), inter[c], -c)
                      for c in inter if len(comunidades[c]) > inter[c]]  # traz nó novo
        for c in inter:  # comunidades já contidas em S: marca como usadas
            if len(comunidades[c]) == inter[c]:
                usados.add(c)
        if candidatos:
            _, _, negc = max(candidatos)
            prox = -negc
        else:  # nenhuma comunidade sobreposta restante: nova semente aleatória
            restantes = [c for c in elegiveis if c not in usados]
            if not restantes:
                break
            prox = rng.choice(restantes)
        ordem.append(prox)
        usados.add(prox)
        S.update(comunidades[prox])
    return semente, ordem, S


def avaliar_prefixos(ordem, comunidades, arestas_sup: np.ndarray,
                     alvo_min: int, alvo_max: int):
    """Para cada prefixo k da lista de comunidades, calcula |V|, |E| induzidos
    e o tamanho da LCC. Escolhe o maior k com LCC <= alvo_max."""
    u_all, v_all = decodificar_chave(arestas_sup)
    max_id = int(max(u_all.max(initial=0), v_all.max(initial=0),
                     max(max(comunidades[c]) for c in ordem))) + 1
    # posição (índice de prefixo) em que cada nó entra no conjunto
    entrada = np.full(max_id, np.iinfo(np.int32).max, dtype=np.int32)
    for k, c in enumerate(ordem):
        idx = np.asarray(comunidades[c])
        entrada[idx] = np.minimum(entrada[idx], k)
    k_aresta = np.maximum(entrada[u_all], entrada[v_all])  # prefixo em que a aresta surge

    tabela = []
    G = nx.Graph()
    escolhido = None
    for k, c in enumerate(ordem):
        G.add_nodes_from(comunidades[c])
        sel = k_aresta == k
        G.add_edges_from(zip(u_all[sel].tolist(), v_all[sel].tolist()))
        lcc = max(nx.connected_components(G), key=len)
        tabela.append({"k_comunidades": k + 1, "comunidade_adicionada": int(c),
                       "V_induzido": G.number_of_nodes(), "E_induzido": G.number_of_edges(),
                       "n_componentes": nx.number_connected_components(G),
                       "V_lcc": len(lcc)})
        if len(lcc) <= alvo_max:
            escolhido = k + 1
        else:
            break
    if escolhido is None:
        escolhido = 1
    return escolhido, pd.DataFrame(tabela), entrada


# ----------------------------------------------------------------------------
# Leitura em streaming
# ----------------------------------------------------------------------------
def contar_comentarios(caminho: Path) -> tuple[int, list[str]]:
    n, linhas = 0, []
    with gzip.open(caminho, "rt") as fh:
        for linha in fh:
            if linha.startswith("#"):
                n += 1
                linhas.append(linha.rstrip())
            else:
                break
    return n, linhas


def passada_unica(caminho: Path, mascara_sup: np.ndarray, tmpdir: Path):
    """Lê o arquivo em blocos e devolve as estatísticas de auditoria e as
    chaves das arestas induzidas pelo superconjunto (sem laços)."""
    buckets = [open(tmpdir / f"b{i:02d}.bin", "wb") for i in range(N_BUCKETS)]
    grau = np.zeros(4_000_000, dtype=np.int64)
    n_linhas = n_laços = n_u_maior_v = 0
    max_id = -1
    min_id = None
    ordenado = True
    ultimo = None
    induzidas = []
    t0 = time.time()
    leitor = pd.read_csv(caminho, sep=r"\s+", comment="#", header=None,
                         names=["u", "v"], dtype=np.int64, chunksize=CHUNK,
                         engine="c")
    for i, df in enumerate(leitor):
        u = df["u"].to_numpy()
        v = df["v"].to_numpy()
        n_linhas += len(u)
        laço = u == v
        n_laços += int(laço.sum())
        n_u_maior_v += int((u > v).sum())
        m = int(max(u.max(), v.max()))
        mi = int(min(u.min(), v.min()))
        max_id = max(max_id, m)
        min_id = mi if min_id is None else min(min_id, mi)
        if m >= len(grau):
            grau = np.concatenate([grau, np.zeros(m + 1 - len(grau) + 1_000_000, np.int64)])
        # grau "bruto" (cada linha conta 1 para cada ponta; laço conta 2)
        grau += np.bincount(u, minlength=len(grau))[: len(grau)]
        grau += np.bincount(v, minlength=len(grau))[: len(grau)]

        k = chave_canonica(u, v)
        # o arquivo está ordenado lexicograficamente por (u, v)?
        kk = (u.astype(np.int64) << 32) | v.astype(np.int64)
        if ordenado:
            if ultimo is not None and kk[0] < ultimo:
                ordenado = False
            elif np.any(np.diff(kk) < 0):
                ordenado = False
        ultimo = int(kk[-1])

        b = (k ^ (k >> 17)) % N_BUCKETS
        for j in range(N_BUCKETS):
            k[b == j].tofile(buckets[j])

        # arestas induzidas pelo superconjunto
        dentro = np.zeros(len(u), dtype=bool)
        ok = (u < len(mascara_sup)) & (v < len(mascara_sup))
        dentro[ok] = mascara_sup[u[ok]] & mascara_sup[v[ok]]
        dentro &= ~laço
        if dentro.any():
            induzidas.append(k[dentro])
        log.info("  bloco %3d | linhas lidas: %12s | %.0fs", i + 1, f"{n_linhas:,}", time.time() - t0)
    for fh in buckets:
        fh.close()

    # duplicatas exatas, bucket a bucket
    n_unicas = 0
    for j in range(N_BUCKETS):
        arr = np.fromfile(tmpdir / f"b{j:02d}.bin", dtype=np.int64)
        arr.sort()
        if len(arr):
            n_unicas += int(1 + np.count_nonzero(arr[1:] != arr[:-1]))
        del arr
        (tmpdir / f"b{j:02d}.bin").unlink()

    grau = grau[: max_id + 1]
    presentes = grau > 0
    est = {
        "linhas_de_aresta": n_linhas,
        "vertices_distintos": int(presentes.sum()),
        "id_minimo": min_id,
        "id_maximo": max_id,
        "auto_lacos": n_laços,
        "linhas_com_u_maior_que_v": n_u_maior_v,
        "arestas_nao_direcionadas_unicas_incl_lacos": n_unicas,
        "linhas_duplicadas_(u,v)_ou_(v,u)": n_linhas - n_unicas,
        "arquivo_ordenado_por_(u,v)": ordenado,
        "grau_maximo": int(grau.max()),
        "grau_medio": float(grau[presentes].mean()),
        "grau_mediano": float(np.median(grau[presentes])),
        "segundos_leitura": round(time.time() - t0, 1),
    }
    sup = np.unique(np.concatenate(induzidas)) if induzidas else np.array([], np.int64)
    return est, grau, sup


# ----------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--grafo", type=Path, default=ARQ_GRAFO)
    ap.add_argument("--comunidades", type=Path, default=ARQ_CMTY)
    ap.add_argument("--saida", type=Path, default=PROC)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--alvo-min", type=int, default=5000)
    ap.add_argument("--alvo-max", type=int, default=10000)
    ap.add_argument("--fator-superconjunto", type=float, default=1.6,
                    help="superconjunto = fator × alvo-max (margem p/ escolher o prefixo)")
    ap.add_argument("--tam-max-cmty", type=int, default=2000)
    ap.add_argument("--tmp", type=Path, default=None,
                    help="pasta temporária p/ os buckets (~1 GB). Padrão: temp do sistema")
    ap.add_argument("--sem-verificacao", action="store_true",
                    help="não falhar se |V|,|E| diferirem dos valores do SNAP")
    a = ap.parse_args()

    a.saida.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S",
                        handlers=[logging.StreamHandler(sys.stdout),
                                  logging.FileHandler(a.saida / "log_tratamento.txt", "w", "utf-8")])
    t_ini = time.time()
    log.info("=== TP Grafos SNAP — G27 com-Orkut — tratamento ===")
    log.info("Python %s | numpy %s | pandas %s | networkx %s | %s",
             platform.python_version(), np.__version__, pd.__version__, nx.__version__, platform.platform())
    log.info("Arquivo do grafo: %s (%s bytes)", a.grafo.name, f"{a.grafo.stat().st_size:,}")
    log.info("Arquivo de comunidades: %s (%s bytes)", a.comunidades.name, f"{a.comunidades.stat().st_size:,}")
    log.info("Parâmetros: seed=%d alvo=[%d, %d] tam_max_cmty=%d", a.seed, a.alvo_min, a.alvo_max, a.tam_max_cmty)

    # --- comunidades e superconjunto ---
    n_com, cab = contar_comentarios(a.grafo)
    log.info("Cabeçalho do arquivo (%d linhas de comentário):", n_com)
    for l in cab:
        log.info("    %s", l)
    comunidades = ler_comunidades(a.comunidades)
    tams = np.array([len(c) for c in comunidades])
    log.info("Comunidades top5000 lidas: %d (tamanho mín=%d, mediana=%.0f, máx=%d)",
             len(comunidades), tams.min(), np.median(tams), tams.max())
    alvo_sup = int(a.fator_superconjunto * a.alvo_max)
    semente, ordem, S = selecionar_comunidades(comunidades, a.seed, alvo_sup, a.tam_max_cmty)
    log.info("Comunidade-semente (seed=%d): #%d (%d nós)", a.seed, semente, len(comunidades[semente]))
    log.info("Superconjunto: %d comunidades, %d vértices", len(ordem), len(S))
    mascara = np.zeros(max(S) + 1, dtype=bool)
    mascara[np.fromiter(S, dtype=np.int64)] = True

    # --- leitura única em streaming ---
    tmpbase = a.tmp or Path(tempfile.gettempdir())
    tmpbase.mkdir(parents=True, exist_ok=True)
    tmpdir = Path(tempfile.mkdtemp(prefix="orkut_buckets_", dir=tmpbase))
    log.info("Leitura em streaming (blocos de %s linhas); buckets temporários em %s", f"{CHUNK:,}", tmpdir)
    try:
        est, grau_orig, sup = passada_unica(a.grafo, mascara, tmpdir)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    est["linhas_de_comentario"] = n_com
    log.info("Auditoria do grafo original:")
    for k, v in est.items():
        log.info("    %-45s %s", k, f"{v:,}" if isinstance(v, int) and not isinstance(v, bool) else v)

    ok_v = est["vertices_distintos"] == ESPERADO_V
    ok_e = est["linhas_de_aresta"] == ESPERADO_E
    log.info("Conferência com o SNAP: |V| %s (esperado %s) | |E| %s (esperado %s)",
             "OK" if ok_v else "DIVERGE", f"{ESPERADO_V:,}", "OK" if ok_e else "DIVERGE", f"{ESPERADO_E:,}")

    # --- escolha do prefixo e LCC ---
    log.info("Arestas induzidas no superconjunto (simples, sem laços): %s", f"{len(sup):,}")
    k_esc, tab_pref, entrada = avaliar_prefixos(ordem, comunidades, sup, a.alvo_min, a.alvo_max)
    sel_cmty = ordem[:k_esc]
    nos_ind = set().union(*(comunidades[c] for c in sel_cmty))
    u_all, v_all = decodificar_chave(sup)
    m = (entrada[u_all] < k_esc) & (entrada[v_all] < k_esc)
    G = nx.Graph()
    G.add_nodes_from(nos_ind)
    G.add_edges_from(zip(u_all[m].tolist(), v_all[m].tolist()))
    n_ind, m_ind, c_ind = G.number_of_nodes(), G.number_of_edges(), nx.number_connected_components(G)
    lcc_nos = max(nx.connected_components(G), key=len)
    H = G.subgraph(lcc_nos).copy()
    n_lcc, m_lcc = H.number_of_nodes(), H.number_of_edges()
    log.info("Prefixo escolhido: %d comunidades -> induzido |V|=%s |E|=%s (%d componentes) -> LCC |V|=%s |E|=%s",
             k_esc, f"{n_ind:,}", f"{m_ind:,}", c_ind, f"{n_lcc:,}", f"{m_lcc:,}")
    if not (a.alvo_min <= n_lcc <= a.alvo_max):
        log.warning("ATENÇÃO: LCC fora do intervalo alvo [%d, %d]. Ajuste --seed/--tam-max-cmty.",
                    a.alvo_min, a.alvo_max)

    # --- saídas ---
    arestas = pd.DataFrame(sorted((min(x, y), max(x, y)) for x, y in H.edges()), columns=["src", "dst"])
    arestas.to_csv(a.saida / "orkut_g27_lcc_arestas.csv", index=False)
    with gzip.open(a.saida / "orkut_g27_lcc.edgelist.gz", "wt") as fh:
        fh.write(f"# G27 com-Orkut — LCC do subgrafo induzido por {k_esc} comunidades (seed={a.seed})\n")
        fh.write(f"# Nodes: {n_lcc} Edges: {m_lcc} (não direcionado, simples)\n")
        for s, d in arestas.itertuples(index=False):
            fh.write(f"{s}\t{d}\n")

    pertence = defaultdict(list)
    for c in sel_cmty:
        for n in comunidades[c]:
            if n in lcc_nos:
                pertence[n].append(c)
    nos = pd.DataFrame({"node_id": sorted(lcc_nos)})
    nos["grau_original"] = grau_orig[nos["node_id"]]
    nos["grau_subgrafo"] = [H.degree(n) for n in nos["node_id"]]
    nos["comunidade_principal"] = [pertence[n][0] for n in nos["node_id"]]
    nos["n_comunidades"] = [len(pertence[n]) for n in nos["node_id"]]
    nos["comunidades"] = [";".join(map(str, pertence[n])) for n in nos["node_id"]]
    nos.to_csv(a.saida / "orkut_g27_lcc_nos.csv", index=False)

    cm = pd.DataFrame([{"ordem": i + 1, "indice_comunidade_top5000": c,
                        "tamanho_original": len(comunidades[c]),
                        "nos_na_lcc": sum(1 for n in comunidades[c] if n in lcc_nos)}
                       for i, c in enumerate(sel_cmty)])
    cm.to_csv(a.saida / "orkut_g27_comunidades_selecionadas.csv", index=False)
    tab_pref.to_csv(a.saida / "orkut_g27_prefixos_avaliados.csv", index=False)

    graus = np.array([d for _, d in H.degree()])
    dens = 2 * m_lcc / (n_lcc * (n_lcc - 1))
    clus = nx.average_clustering(H)
    final = {"vertices": n_lcc, "arestas": m_lcc, "grau_medio": round(float(graus.mean()), 3),
             "grau_max": int(graus.max()), "grau_min": int(graus.min()),
             "densidade": dens, "clusterizacao_media": round(clus, 4),
             "componentes": nx.number_connected_components(H), "auto_lacos": nx.number_of_selfloops(H)}
    log.info("Grafo final: %s", final)

    # tabela de tratamentos
    tr = [
        ("Leitura do bruto", f"{a.grafo.name}, leitura em streaming (blocos de {CHUNK:,} linhas)",
         f"{n_com} linhas de comentário; {est['linhas_de_aresta']:,} linhas de aresta"),
        ("Auditoria de vértices", "vértices distintos que aparecem em alguma aresta",
         f"{est['vertices_distintos']:,} (SNAP: {ESPERADO_V:,}; ids {est['id_minimo']}–{est['id_maximo']:,})"),
        ("Auditoria de arestas", "linhas de aresta vs. valor oficial", f"{est['linhas_de_aresta']:,} (SNAP: {ESPERADO_E:,})"),
        ("Auto-laços", "linhas com u = v", f"{est['auto_lacos']:,}"),
        ("Arestas duplicadas", "pares repetidos, considerando (u,v) ≡ (v,u)",
         f"{est['linhas_duplicadas_(u,v)_ou_(v,u)']:,}"),
        ("Orientação", "linhas com u > v (arquivo grava cada aresta uma vez?)", f"{est['linhas_com_u_maior_que_v']:,}"),
        ("Seleção de comunidades", f"semente aleatória fixa (seed={a.seed}) + expansão gulosa por sobreposição |C∩S|/|C|; "
         f"comunidades com > {a.tam_max_cmty} nós ignoradas",
         f"semente #{semente}; {k_esc} comunidades selecionadas"),
        ("Subgrafo induzido", "todas as arestas do grafo original com as duas pontas no conjunto",
         f"|V|={n_ind:,}; |E|={m_ind:,}; {c_ind} componentes"),
        ("Simplificação", "remoção de auto-laços e duplicatas; grafo não direcionado, sem pesos",
         f"{nx.number_of_selfloops(H)} laços; arestas únicas por par"),
        ("Maior componente conexa", "mantida só a LCC (distâncias, diâmetro, small-world)",
         f"|V|={n_lcc:,}; |E|={m_lcc:,}"),
    ]
    tdf = pd.DataFrame(tr, columns=["Tratamento", "Descrição", "Resultado"])
    tdf.to_csv(a.saida / "tabela_tratamentos.csv", index=False)
    esc = lambda s: (str(s).replace("\\", "\\textbackslash{}").replace("&", "\\&").replace("%", "\\%")
                     .replace("_", "\\_").replace("#", "\\#").replace("|", "$|$").replace("∩", "$\\cap$")
                     .replace("≡", "$\\equiv$").replace("–", "--").replace(">", "$>$"))
    with open(a.saida / "tabela_tratamentos.tex", "w", encoding="utf-8") as fh:
        fh.write("\\begin{table}[ht]\n\\centering\n\\caption{Tratamentos aplicados ao grafo G27 (com-Orkut).}\n"
                 "\\label{tab:tratamentos}\n\\small\n\\begin{tabular}{p{3cm}p{6cm}p{4.5cm}}\n\\hline\n"
                 "Tratamento & Descrição & Resultado \\\\\n\\hline\n")
        for r in tr:
            fh.write(" & ".join(esc(x) for x in r) + " \\\\\n")
        fh.write("\\hline\n\\end{tabular}\n\\end{table}\n")

    resumo = {"grafo": "G27 com-Orkut", "fonte": "https://snap.stanford.edu/data/com-Orkut.html",
              "parametros": vars(a) | {"grafo": str(a.grafo.name), "comunidades": str(a.comunidades.name),
                                        "saida": str(a.saida.relative_to(RAIZ)) if a.saida.is_relative_to(RAIZ) else str(a.saida),
                                        "tmp": None},
              "auditoria_original": est, "conferencia_snap": {"V_ok": ok_v, "E_ok": ok_e},
              "selecao": {"comunidade_semente": semente, "k_comunidades": k_esc,
                          "comunidades": [int(c) for c in sel_cmty]},
              "induzido": {"vertices": n_ind, "arestas": m_ind, "componentes": c_ind},
              "final_lcc": final, "segundos_total": round(time.time() - t_ini, 1)}
    with open(a.saida / "resumo_tratamento.json", "w", encoding="utf-8") as fh:
        json.dump(resumo, fh, indent=2, ensure_ascii=False, default=str)

    md = [f"# Resumo do tratamento — G27 com-Orkut\n",
          f"Gerado por `scripts/tratamento_orkut.py` em {time.strftime('%Y-%m-%d %H:%M')} "
          f"({resumo['segundos_total']:.0f} s). Log completo em `log_tratamento.txt`.\n",
          "## Auditoria do grafo original\n", "| Medida | Valor |", "|---|---|"]
    md += [f"| {k} | {v:,} |" if isinstance(v, int) and not isinstance(v, bool) else f"| {k} | {v} |"
           for k, v in est.items()]
    md += [f"\nConferência com a página do SNAP: |V| **{'OK' if ok_v else 'DIVERGE'}**, "
           f"|E| **{'OK' if ok_e else 'DIVERGE'}**.\n",
           "## Tratamentos\n", "| Tratamento | Descrição | Resultado |", "|---|---|---|"]
    md += [f"| {t} | {d} | {r} |" for t, d, r in tr]
    md += ["\n## Grafo final (data/processed)\n", "| Medida | Valor |", "|---|---|"]
    md += [f"| {k} | {v} |" for k, v in final.items()]
    md += ["\n## Justificativa metodológica\n",
           "- O grafo completo (117 M arestas) não cabe confortavelmente em memória com NetworkX e inviabiliza "
           "algoritmos O(n²)/O(n³) do TP; por isso foi reduzido.",
           "- Usou-se um **subgrafo induzido por comunidades ground-truth** do próprio SNAP (top5000), e não "
           "amostragem aleatória de arestas: o subgrafo induzido preserva **todas** as arestas entre os vértices "
           "escolhidos, então graus locais, triângulos e estrutura comunitária são mantidos.",
           "- As comunidades foram encadeadas por **sobreposição** (|C∩S|/|C|) a partir de uma semente aleatória "
           "**fixa**, gerando uma região coesa e reprodutível da rede; a LCC garante conectividade para "
           "distâncias, diâmetro e small-world.",
           "- **Viés conhecido:** a amostra é centrada em comunidades densas, então clusterização e densidade "
           "tendem a ser maiores que no Orkut completo, e o grau no subgrafo é menor que o original "
           "(coluna `grau_original` em `orkut_g27_lcc_nos.csv` permite comparar).",
           "- O Orkut é não direcionado e sem pesos: não há versão dirigida a manter; algoritmos ponderados "
           "(Dijkstra, AGM) precisarão de pesos declarados explicitamente na análise.",
           "\n## Arquivos\n",
           "- `orkut_g27_lcc_arestas.csv` — arestas (src,dst), ids originais do Orkut",
           "- `orkut_g27_lcc.edgelist.gz` — mesma lista em formato edgelist (NetworkX `read_edgelist`)",
           "- `orkut_g27_lcc_nos.csv` — nós, grau original, grau no subgrafo e comunidades",
           "- `orkut_g27_comunidades_selecionadas.csv`, `orkut_g27_prefixos_avaliados.csv` — rastreabilidade da seleção",
           "- `tabela_tratamentos.csv` / `.tex`, `resumo_tratamento.json`, `log_tratamento.txt`"]
    (a.saida / "resumo_tratamento.md").write_text("\n".join(md) + "\n", encoding="utf-8")

    for f in sorted(a.saida.iterdir()):
        log.info("  saída: %-45s %10s bytes", f.name, f"{f.stat().st_size:,}")
    log.info("Concluído em %.0f s.", time.time() - t_ini)
    if not (ok_v and ok_e) and not a.sem_verificacao:
        log.error("Contagens divergem do SNAP (use --sem-verificacao para ignorar).")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
