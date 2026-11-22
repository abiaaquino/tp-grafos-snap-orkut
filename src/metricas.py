"""Parte I — medidas estruturais (implementação própria)."""
from __future__ import annotations

import numpy as np

from .algoritmos import componentes_conexas


def resumo_graus(graus: np.ndarray) -> dict:
    g = np.asarray(graus)
    return {"grau_min": int(g.min()), "grau_max": int(g.max()), "grau_medio": float(g.mean()),
            "grau_mediano": float(np.median(g)), "grau_desvio": float(g.std(ddof=1))}


def distribuicao_graus(graus: np.ndarray):
    """(k, P(k)) só para k >= 1 (k = 0 quebra o log-log) e CCDF P(K >= k)."""
    g = np.asarray(graus)
    g = g[g > 0]
    k, c = np.unique(g, return_counts=True)
    pk = c / c.sum()
    ccdf = 1.0 - np.concatenate([[0.0], np.cumsum(pk)[:-1]])
    return k, pk, ccdf


def densidade(n: int, m: int) -> float:
    return 2.0 * m / (n * (n - 1)) if n > 1 else 0.0


def componentes(adj) -> dict:
    rot, tam = componentes_conexas(adj)
    tam = sorted(tam, reverse=True)
    return {"n_componentes": len(tam), "tamanhos": tam, "maior": tam[0] if tam else 0, "rotulos": rot}


def triangulos_por_vertice(adj) -> np.ndarray:
    """Contagem de triângulos com orientação por grau (cada triângulo é contado
    uma vez a partir do vértice de menor posto): O(m·√m)."""
    n = len(adj)
    grau = [len(a) for a in adj]
    posto = sorted(range(n), key=lambda x: (grau[x], x))
    pos = [0] * n
    for i, x in enumerate(posto):
        pos[x] = i
    fora = [set(y for y in adj[x] if pos[y] > pos[x]) for x in range(n)]
    t = np.zeros(n, dtype=np.int64)
    for x in range(n):
        fx = fora[x]
        for y in fx:
            comuns = fx & fora[y]
            if comuns:
                c = len(comuns)
                t[x] += c
                t[y] += c
                for z in comuns:
                    t[z] += 1
    return t


def clusterizacao(adj, tri: np.ndarray | None = None) -> dict:
    tri = triangulos_por_vertice(adj) if tri is None else tri
    k = np.array([len(a) for a in adj], dtype=float)
    pares = k * (k - 1) / 2
    local = np.divide(tri, pares, out=np.zeros_like(pares), where=pares > 0)
    return {"triangulos": int(tri.sum() // 3), "clusterizacao_media": float(local.mean()),
            "transitividade": float(tri.sum() / pares.sum()) if pares.sum() else 0.0,
            "local": local, "tri_por_vertice": tri}
