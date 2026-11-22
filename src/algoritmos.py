"""Parte II — algoritmos clássicos, implementação própria (Python puro / numpy).

Convenções: ``adj`` é lista de adjacência (list[list[int]]); ``wadj`` é lista
ponderada (list[list[(vizinho, peso)]]); ``arestas`` é lista de (u, v, peso).
Complexidades (n = |V|, m = |E|):
  BFS, DFS, componentes, Euler, Tarjan .............. O(n + m)
  Dijkstra (heap binário) ........................... O((n + m) log n)
  Bellman-Ford ...................................... O(n·m)  (efetivo O(k·m), k = rodadas)
  Floyd-Warshall .................................... O(n³) tempo, O(n²) memória
  Prim (heap) ....................................... O(m log n)
  Kruskal (união-busca) ............................. O(m log m)
"""
from __future__ import annotations

from collections import deque
import heapq
import math

import numpy as np

INF = math.inf


# ---------------------------------------------------------------- buscas
def bfs(adj, s: int):
    """Retorna (distância em saltos, pai). -1 = não alcançado."""
    n = len(adj)
    dist = [-1] * n
    pai = [-1] * n
    dist[s] = 0
    fila = deque([s])
    while fila:
        x = fila.popleft()
        dx = dist[x] + 1
        for y in adj[x]:
            if dist[y] < 0:
                dist[y] = dx
                pai[y] = x
                fila.append(y)
    return dist, pai


def dfs(adj, s: int):
    """DFS iterativa. Retorna (ordem de descoberta, pai)."""
    n = len(adj)
    vis = bytearray(n)
    pai = [-1] * n
    ordem = []
    pilha = [(s, iter(adj[s]))]
    vis[s] = 1
    ordem.append(s)
    while pilha:
        x, it = pilha[-1]
        for y in it:
            if not vis[y]:
                vis[y] = 1
                pai[y] = x
                ordem.append(y)
                pilha.append((y, iter(adj[y])))
                break
        else:
            pilha.pop()
    return ordem, pai


def componentes_conexas(adj):
    """Rótulo de componente por vértice (BFS). Retorna (rotulo, tamanhos)."""
    n = len(adj)
    rot = [-1] * n
    tamanhos = []
    for s in range(n):
        if rot[s] >= 0:
            continue
        c = len(tamanhos)
        rot[s] = c
        fila = deque([s])
        t = 0
        while fila:
            x = fila.popleft()
            t += 1
            for y in adj[x]:
                if rot[y] < 0:
                    rot[y] = c
                    fila.append(y)
        tamanhos.append(t)
    return rot, tamanhos


# ---------------------------------------------------------------- Euler
def verificar_euler(adj):
    """Grafo não direcionado: euleriano ⇔ conexo (nos vértices com grau > 0) e
    todos os graus pares; semi-euleriano ⇔ conexo e exatamente 2 graus ímpares."""
    graus = [len(a) for a in adj]
    impares = [i for i, g in enumerate(graus) if g % 2]
    com_aresta = [i for i, g in enumerate(graus) if g > 0]
    conexo = True
    if com_aresta:
        dist, _ = bfs(adj, com_aresta[0])
        conexo = all(dist[i] >= 0 for i in com_aresta)
    return {
        "conexo": conexo,
        "n_impares": len(impares),
        "euleriano": conexo and len(impares) == 0,
        "semi_euleriano": conexo and len(impares) == 2,
        "impares_exemplo": impares[:10],
    }


def circuito_euleriano(adj):
    """Hierholzer (iterativo). Exige grafo euleriano ou semi-euleriano."""
    info = verificar_euler(adj)
    if not (info["euleriano"] or info["semi_euleriano"]):
        return None
    n = len(adj)
    usados = set()
    ptr = [0] * n
    inicio = next((i for i in range(n) if len(adj[i]) % 2), next(i for i in range(n) if adj[i]))
    pilha, caminho = [inicio], []
    while pilha:
        x = pilha[-1]
        while ptr[x] < len(adj[x]):
            y = adj[x][ptr[x]]
            ptr[x] += 1
            e = (x, y) if x < y else (y, x)
            if e not in usados:
                usados.add(e)
                pilha.append(y)
                break
        else:
            caminho.append(pilha.pop())
    return caminho[::-1]


# ---------------------------------------------------------------- caminhos mínimos
def dijkstra(wadj, s: int):
    n = len(wadj)
    dist = [INF] * n
    pai = [-1] * n
    dist[s] = 0.0
    heap = [(0.0, s)]
    while heap:
        d, x = heapq.heappop(heap)
        if d > dist[x]:
            continue
        for y, w in wadj[x]:
            nd = d + w
            if nd < dist[y]:
                dist[y] = nd
                pai[y] = x
                heapq.heappush(heap, (nd, y))
    return dist, pai


def bellman_ford(n: int, arcos, s: int):
    """``arcos`` = lista de arcos dirigidos (u, v, w); para grafo não direcionado
    passe os dois sentidos. Parada antecipada quando uma rodada não relaxa nada.
    Retorna (dist, rodadas_executadas, tem_ciclo_negativo)."""
    dist = [INF] * n
    dist[s] = 0.0
    rodadas = 0
    for _ in range(n - 1):
        rodadas += 1
        mudou = False
        for a, b, w in arcos:
            da = dist[a]
            if da + w < dist[b]:
                dist[b] = da + w
                mudou = True
        if not mudou:
            break
    neg = any(dist[a] + w < dist[b] for a, b, w in arcos if dist[a] < INF)
    return dist, rodadas, neg


def arcos_nao_direcionados(u, v, w):
    out = []
    for a, b, p in zip(list(u), list(v), list(w)):
        out.append((a, b, p))
        out.append((b, a, p))
    return out


def matriz_pesos(n: int, u, v, w) -> np.ndarray:
    D = np.full((n, n), np.inf)
    np.fill_diagonal(D, 0.0)
    D[u, v] = np.minimum(D[u, v], w)
    D[v, u] = np.minimum(D[v, u], w)
    return D


def floyd_warshall_numpy(D: np.ndarray) -> np.ndarray:
    """Vetorizado: n iterações de uma atualização O(n²) feita em C pelo numpy."""
    D = D.copy()
    for k in range(D.shape[0]):
        np.minimum(D, D[:, k:k + 1] + D[k:k + 1, :], out=D)
    return D


def floyd_warshall_python(D) -> list:
    """Versão em Python puro (três laços), para comparar com a vetorizada."""
    D = [list(l) for l in (D.tolist() if hasattr(D, "tolist") else D)]
    n = len(D)
    for k in range(n):
        Dk = D[k]
        for i in range(n):
            Di = D[i]
            dik = Di[k]
            if dik == INF:
                continue
            for j in range(n):
                nd = dik + Dk[j]
                if nd < Di[j]:
                    Di[j] = nd
    return D


# ---------------------------------------------------------------- Tarjan
def tarjan_pontes_articulacoes(adj):
    """Tarjan (low-link) iterativo para grafo não direcionado.
    Retorna (lista de pontes (u<v), conjunto de pontos de articulação)."""
    n = len(adj)
    disc = [-1] * n
    low = [0] * n
    pontes, artic = [], set()
    t = 0
    for r in range(n):
        if disc[r] >= 0:
            continue
        disc[r] = low[r] = t
        t += 1
        filhos_raiz = 0
        pilha = [(r, -1, iter(adj[r]))]
        while pilha:
            x, p, it = pilha[-1]
            avancou = False
            for y in it:
                if disc[y] < 0:
                    disc[y] = low[y] = t
                    t += 1
                    pilha.append((y, x, iter(adj[y])))
                    if x == r:
                        filhos_raiz += 1
                    avancou = True
                    break
                elif y != p:
                    if disc[y] < low[x]:
                        low[x] = disc[y]
            if avancou:
                continue
            pilha.pop()
            if p >= 0:
                if low[x] < low[p]:
                    low[p] = low[x]
                if low[x] > disc[p]:
                    pontes.append((min(p, x), max(p, x)))
                if p != r and low[x] >= disc[p]:
                    artic.add(p)
        if filhos_raiz >= 2:
            artic.add(r)
    return pontes, artic


# ---------------------------------------------------------------- AGM
class UniaoBusca:
    def __init__(self, n):
        self.p = list(range(n))
        self.r = [0] * n

    def achar(self, x):
        p = self.p
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x

    def unir(self, a, b):
        a, b = self.achar(a), self.achar(b)
        if a == b:
            return False
        if self.r[a] < self.r[b]:
            a, b = b, a
        self.p[b] = a
        if self.r[a] == self.r[b]:
            self.r[a] += 1
        return True


def kruskal(n: int, u, v, w):
    ordem = sorted(range(len(w)), key=w.__getitem__)
    ub = UniaoBusca(n)
    arv, total = [], 0.0
    for i in ordem:
        if ub.unir(u[i], v[i]):
            arv.append((u[i], v[i], w[i]))
            total += w[i]
            if len(arv) == n - 1:
                break
    return arv, total


def prim(wadj, s: int = 0):
    """Prim com heap binário (versão "lazy"). Em grafo desconexo devolve a AGM
    da componente de s."""
    n = len(wadj)
    na_arvore = bytearray(n)
    arv, total = [], 0.0
    heap = [(0.0, s, -1)]
    while heap:
        w, x, p = heapq.heappop(heap)
        if na_arvore[x]:
            continue
        na_arvore[x] = 1
        if p >= 0:
            arv.append((p, x, w))
            total += w
        for y, wy in wadj[x]:
            if not na_arvore[y]:
                heapq.heappush(heap, (wy, y, x))
    return arv, total
