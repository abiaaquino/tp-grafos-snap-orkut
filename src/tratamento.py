"""Estrutura de grafo usada em todo o projeto e funções de carga/subamostragem.

O grafo é não direcionado e simples. Os vértices são renumerados para 0..n-1
(o id original do Orkut fica em ``Grafo.ids``). Guardamos três representações:
  * lista de arestas (u < v) em arrays numpy  -> geradores, Kruskal, Bellman-Ford
  * lista de adjacência (list[list[int]])     -> algoritmos em Python puro
  * CSR (indptr, indices)                     -> rotinas vetorizadas (PageRank etc.)
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from pathlib import Path
import random

import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parents[1]
PROC = RAIZ / "data" / "processed"


@dataclass
class Grafo:
    n: int
    u: np.ndarray            # extremidade menor de cada aresta
    v: np.ndarray            # extremidade maior
    ids: np.ndarray          # id original de cada vértice 0..n-1
    adj: list = field(repr=False)
    indptr: np.ndarray = field(repr=False)
    indices: np.ndarray = field(repr=False)
    w: np.ndarray | None = field(default=None, repr=False)   # peso por aresta (opcional)

    @property
    def m(self) -> int:
        return len(self.u)

    def graus(self) -> np.ndarray:
        return np.diff(self.indptr)

    def adj_ponderada(self) -> list:
        """wadj[x] = [(y, peso), ...]; usa peso 1 se o grafo não tiver pesos."""
        w = self.w if self.w is not None else np.ones(self.m)
        wadj = [[] for _ in range(self.n)]
        for a, b, p in zip(self.u.tolist(), self.v.tolist(), w.tolist()):
            wadj[a].append((b, p))
            wadj[b].append((a, p))
        return wadj


def de_arestas(u, v, n: int | None = None, ids=None, w=None) -> Grafo:
    """Monta o Grafo a partir de arestas já em índices 0..n-1. Remove laços e
    duplicatas (mantendo o primeiro peso)."""
    u = np.asarray(u, dtype=np.int64)
    v = np.asarray(v, dtype=np.int64)
    lo, hi = np.minimum(u, v), np.maximum(u, v)
    ok = lo != hi
    lo, hi = lo[ok], hi[ok]
    if w is not None:
        w = np.asarray(w, dtype=float)[ok]
    if n is None:
        n = int(hi.max()) + 1 if len(hi) else 0
    chave = lo * n + hi
    _, pos = np.unique(chave, return_index=True)
    pos.sort()
    lo, hi = lo[pos], hi[pos]
    if w is not None:
        w = w[pos]
    origem = np.concatenate([lo, hi])
    destino = np.concatenate([hi, lo])
    ordem = np.lexsort((destino, origem))
    origem, destino = origem[ordem], destino[ordem]
    indptr = np.zeros(n + 1, dtype=np.int64)
    np.add.at(indptr, origem + 1, 1)
    indptr = np.cumsum(indptr)
    indices = destino.astype(np.int64)
    adj = [indices[indptr[i]:indptr[i + 1]].tolist() for i in range(n)]
    ids = np.arange(n) if ids is None else np.asarray(ids)
    return Grafo(n=n, u=lo, v=hi, ids=ids, adj=adj, indptr=indptr, indices=indices, w=w)


def carregar_lcc(pasta: Path = PROC) -> tuple[Grafo, pd.DataFrame]:
    """Lê o dataset tratado (LCC do subgrafo induzido por comunidades)."""
    ar = pd.read_csv(pasta / "orkut_g27_lcc_arestas.csv")
    nos = pd.read_csv(pasta / "orkut_g27_lcc_nos.csv").sort_values("node_id").reset_index(drop=True)
    ids = nos["node_id"].to_numpy()
    idx = pd.Series(np.arange(len(ids)), index=ids)
    G = de_arestas(idx[ar["src"]].to_numpy(), idx[ar["dst"]].to_numpy(), n=len(ids), ids=ids)
    return G, nos


def vizinhos_comuns_por_aresta(G: Grafo) -> np.ndarray:
    conj = [set(a) for a in G.adj]
    return np.array([len(conj[a] & conj[b]) for a, b in zip(G.u.tolist(), G.v.tolist())], dtype=np.int64)


def pesos_embeddedness(G: Grafo) -> np.ndarray:
    """Peso derivado dos dados: w(u,v) = 1 / (1 + |N(u) ∩ N(v)|).

    Arestas muito "embutidas" (muitos amigos em comum, laços fortes) ficam
    curtas; arestas-ponte entre grupos ficam longas (w = 1). Serve para dar
    sentido a Dijkstra, Bellman-Ford, Floyd-Warshall e à AGM num grafo que é
    originalmente não ponderado."""
    return 1.0 / (1.0 + vizinhos_comuns_por_aresta(G))


def com_pesos(G: Grafo, w) -> Grafo:
    return Grafo(G.n, G.u, G.v, G.ids, G.adj, G.indptr, G.indices, np.asarray(w, float))


def subgrafo_induzido(G: Grafo, vertices) -> Grafo:
    vertices = np.asarray(sorted(set(int(x) for x in vertices)), dtype=np.int64)
    mapa = np.full(G.n, -1, dtype=np.int64)
    mapa[vertices] = np.arange(len(vertices))
    ok = (mapa[G.u] >= 0) & (mapa[G.v] >= 0)
    w = G.w[ok] if G.w is not None else None
    return de_arestas(mapa[G.u[ok]], mapa[G.v[ok]], n=len(vertices), ids=G.ids[vertices], w=w)


def ordem_bfs(G: Grafo, raiz: int) -> list[int]:
    vis = bytearray(G.n)
    vis[raiz] = 1
    fila, ordem = deque([raiz]), []
    while fila:
        x = fila.popleft()
        ordem.append(x)
        for y in G.adj[x]:
            if not vis[y]:
                vis[y] = 1
                fila.append(y)
    return ordem


def subgrafo_bfs(G: Grafo, n_alvo: int, raiz: int | None = None, seed: int = 42) -> Grafo:
    """Subgrafo induzido pelos n_alvo primeiros vértices de uma BFS (conexo por
    construção). Raiz aleatória com semente fixa, para não começar sempre no hub."""
    if raiz is None:
        raiz = random.Random(seed).randrange(G.n)
    return subgrafo_induzido(G, ordem_bfs(G, raiz)[:n_alvo])
