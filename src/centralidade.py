"""Centralidades próprias: PageRank (iteração de potência) e betweenness de
Brandes (exato ou amostrado por fontes, escala n/k)."""
from __future__ import annotations

from collections import deque
import random

import numpy as np


def pagerank(G, alpha: float = 0.85, tol: float = 1e-10, max_iter: int = 200):
    """Iteração de potência em CSR. Grafo não direcionado: cada aresta vale nos
    dois sentidos. Vértices sem vizinhos redistribuem uniformemente."""
    n = G.n
    grau = G.graus().astype(float)
    origem = np.repeat(np.arange(n), np.diff(G.indptr))
    destino = G.indices
    pr = np.full(n, 1.0 / n)
    for it in range(max_iter):
        contrib = np.divide(pr, grau, out=np.zeros(n), where=grau > 0)
        novo = np.bincount(destino, weights=contrib[origem], minlength=n)
        pendente = pr[grau == 0].sum()
        novo = alpha * (novo + pendente / n) + (1 - alpha) / n
        if np.abs(novo - pr).sum() < tol:
            pr = novo
            break
        pr = novo
    return pr / pr.sum(), it + 1


def brandes(adj, fontes=None, normalizar: bool = True):
    """Betweenness de Brandes (não ponderado). Com ``fontes`` amostradas,
    o acumulado é multiplicado por n/k (estimador não viesado)."""
    n = len(adj)
    fontes = range(n) if fontes is None else fontes
    k = len(fontes)
    cb = [0.0] * n
    for s in fontes:
        pilha = []
        pred = [[] for _ in range(n)]
        sigma = [0] * n
        sigma[s] = 1
        dist = [-1] * n
        dist[s] = 0
        fila = deque([s])
        while fila:
            x = fila.popleft()
            pilha.append(x)
            dx = dist[x] + 1
            sx = sigma[x]
            for y in adj[x]:
                if dist[y] < 0:
                    dist[y] = dx
                    fila.append(y)
                if dist[y] == dx:
                    sigma[y] += sx
                    pred[y].append(x)
        delta = [0.0] * n
        while pilha:
            y = pilha.pop()
            coef = (1.0 + delta[y]) / sigma[y]
            for x in pred[y]:
                delta[x] += sigma[x] * coef
            if y != s:
                cb[y] += delta[y]
    cb = np.array(cb) * (n / k) / 2.0          # não direcionado: cada par contado 2x
    if normalizar and n > 2:
        cb /= (n - 1) * (n - 2) / 2.0
    return cb


def brandes_amostrado(adj, k: int, seed: int = 42):
    fontes = random.Random(seed).sample(range(len(adj)), k)
    return brandes(adj, fontes=fontes)
