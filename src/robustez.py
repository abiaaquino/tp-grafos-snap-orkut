"""Parte III — robustez: remoção aleatória × remoção dos mais centrais.

n = |V|, r = ⌈0,05·n⌉. Métricas após a remoção:
  A) S = |maior componente| / (n − r)       (fração do que restou na LCC)
  B) c = número de componentes conexas
  C) L = caminho médio na maior componente (amostrado por fontes, igraph)
  D) E = eficiência global média 1/d(u,v) (amostrada; pares desconexos valem 0)
Ataques usam ranking FIXO calculado no grafo intacto (grau, betweenness, PageRank).
"""
from __future__ import annotations

import math

import numpy as np

from .algoritmos import componentes_conexas
from .tratamento import subgrafo_induzido


def r_remocao(n: int, frac: float = 0.05) -> int:
    return math.ceil(frac * n)


def metricas_apos_remocao(G, removidos, k_fontes: int = 200, seed: int = 0, distancias: bool = True) -> dict:
    import igraph as ig
    rem = set(int(x) for x in removidos)
    restantes = [x for x in range(G.n) if x not in rem]
    H = subgrafo_induzido(G, restantes)
    rot, tam = componentes_conexas(H.adj)
    n_rest = H.n
    maior = max(tam) if tam else 0
    out = {"removidos": len(rem), "n_restante": n_rest, "S": maior / n_rest if n_rest else 0.0,
           "S_sobre_n_original": maior / G.n, "c": len(tam), "maior": maior}
    if not distancias:
        return out
    g = ig.Graph(n=H.n, edges=np.column_stack([H.u, H.v]).tolist())
    rng = np.random.default_rng(seed)
    fontes = rng.choice(H.n, size=min(k_fontes, H.n), replace=False)
    D = np.asarray(g.distances(source=fontes.tolist()), dtype=float)
    inv = np.where(np.isfinite(D) & (D > 0), 1.0 / np.where(D > 0, D, 1), 0.0)
    out["E"] = float(inv.sum() / (len(fontes) * (H.n - 1)))
    rot = np.asarray(rot)
    lcc_id = int(np.argmax(tam))
    na_lcc = rot == lcc_id
    Dl = D[na_lcc[fontes]][:, na_lcc]
    out["L"] = float(Dl[Dl > 0].mean()) if Dl.size else float("nan")
    return out


def remocao_aleatoria(G, r: int, reps: int = 30, seed: int = 0, **kw):
    rng = np.random.default_rng(seed)
    return [metricas_apos_remocao(G, rng.choice(G.n, size=r, replace=False), seed=seed + i, **kw)
            for i in range(reps)]


def remocao_ranking(G, scores, r: int, **kw) -> dict:
    ordem = np.argsort(-np.asarray(scores), kind="stable")
    return metricas_apos_remocao(G, ordem[:r], **kw)


def curva(G, ordem_ou_none, fracs, reps: int = 1, seed: int = 0):
    """S e c em função da fração removida (sem distâncias, para ser rápido)."""
    res = []
    rng = np.random.default_rng(seed)
    for f in fracs:
        r = int(round(f * G.n))
        vals = []
        for i in range(reps):
            rem = rng.choice(G.n, size=r, replace=False) if ordem_ou_none is None else ordem_ou_none[:r]
            vals.append(metricas_apos_remocao(G, rem, distancias=False))
        res.append({"frac": f, "S": float(np.mean([v["S_sobre_n_original"] for v in vals])),
                    "c": float(np.mean([v["c"] for v in vals]))})
    return res
