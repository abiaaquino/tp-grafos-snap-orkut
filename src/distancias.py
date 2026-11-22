"""Distâncias todos-para-todos (diâmetro, raio, caminho médio).

Fazer n BFS em Python puro no grafo de ~10 mil vértices e 173 mil arestas
levaria minutos; usamos o igraph (BFS em C) **em blocos de fontes**, sem nunca
materializar a matriz n×n inteira (que teria 9.917² × 8 B ≈ 787 MB = 750 MiB). A
implementação própria (BFS em Python, ``algoritmos.bfs``) é usada para validar
os resultados num conjunto de fontes e nos testes.
"""
from __future__ import annotations

import numpy as np

from .algoritmos import bfs


def para_igraph(G):
    import igraph as ig
    g = ig.Graph(n=G.n, edges=np.column_stack([G.u, G.v]).tolist(), directed=False)
    if G.w is not None:
        g.es["weight"] = G.w.tolist()
    return g


def todos_para_todos(G, bloco: int = 500, fontes=None, g=None) -> dict:
    """Excentricidades, soma das distâncias e histograma das distâncias a
    partir das ``fontes`` (padrão: todos os vértices -> resultado exato)."""
    g = para_igraph(G) if g is None else g
    fontes = np.arange(G.n) if fontes is None else np.asarray(fontes)
    exc = np.zeros(len(fontes), dtype=np.int64)
    soma, pares, inalc = 0.0, 0, 0
    hist = np.zeros(1, dtype=np.int64)
    for i in range(0, len(fontes), bloco):
        f = fontes[i:i + bloco]
        D = np.asarray(g.distances(source=f.tolist(), mode="all"), dtype=float)
        fin = np.isfinite(D)
        inalc += int((~fin).sum())
        Df = np.where(fin, D, -1).astype(np.int64)
        exc[i:i + len(f)] = Df.max(axis=1)
        vals = Df[Df > 0]
        soma += float(vals.sum())
        pares += int(vals.size)
        h = np.bincount(vals)
        if len(h) > len(hist):
            hist = np.pad(hist, (0, len(h) - len(hist)))
        hist[:len(h)] += h
    return {"excentricidades": exc, "diametro": int(exc.max()), "raio": int(exc.min()),
            "caminho_medio": soma / pares if pares else float("nan"),
            "pares_ordenados": pares, "pares_inalcancaveis": inalc,
            "hist_distancias": hist, "exato": len(fontes) == G.n}


def caminho_medio_amostrado(G, k: int = 300, seed: int = 0, g=None) -> float:
    rng = np.random.default_rng(seed)
    fontes = rng.choice(G.n, size=min(k, G.n), replace=False)
    return todos_para_todos(G, fontes=fontes, g=g)["caminho_medio"]


def validar_com_bfs_propria(G, fontes) -> bool:
    """Confere excentricidades do igraph com a BFS própria nas fontes dadas."""
    r = todos_para_todos(G, fontes=np.asarray(fontes))
    for f, e in zip(fontes, r["excentricidades"]):
        d, _ = bfs(G.adj, int(f))
        if max(d) != e:
            return False
    return True
