"""Parte III — modelos nulos/clássicos e small-world (implementação própria).

Geradores com mesmo n e (aprox.) mesmo m do grafo observado:
  * Erdős–Rényi G(n, m)        -> referência aleatória clássica
  * troca dupla de arestas     -> aleatório que PRESERVA a sequência de graus
  * Barabási–Albert            -> ligação preferencial (cauda pesada)
  * Watts–Strogatz             -> anel regular + religação com prob. p

Small-world (Watts & Strogatz 1998; Humphries & Gurney 2008; Telesford 2011):
  L(G) ≈ L_rand  e  C(G) ≫ C_rand
  σ = (C/C_rand) / (L/L_rand)  > 1  (tipicamente ≫ 1)
  ω = L_rand/L − C/C_latt      ∈ [−1, 1]; ≈ 0 small-world, → −1 reticulado, → 1 aleatório
O cálculo de σ, ω e dos critérios é feito aqui; nenhuma função pronta de
small-world de biblioteca é chamada.
"""
from __future__ import annotations

import math
import random

import numpy as np

from .tratamento import de_arestas


def erdos_renyi_gnm(n: int, m: int, seed: int = 0):
    rng = np.random.default_rng(seed)
    chaves = set()
    while len(chaves) < m:
        falta = m - len(chaves)
        a = rng.integers(0, n, size=int(falta * 1.1) + 10)
        b = rng.integers(0, n, size=len(a))
        lo, hi = np.minimum(a, b), np.maximum(a, b)
        for x, y in zip(lo[lo != hi].tolist(), hi[lo != hi].tolist()):
            chaves.add(x * n + y)
            if len(chaves) == m:
                break
    k = np.fromiter(chaves, dtype=np.int64)
    return de_arestas(k // n, k % n, n=n)


def troca_dupla_arestas(G, trocas_por_aresta: float = 5.0, seed: int = 0):
    """Configuração por trocas duplas (a-b, c-d) -> (a-d, c-b), rejeitando laços
    e arestas múltiplas. Preserva exatamente o grau de cada vértice."""
    rng = random.Random(seed)
    n = G.n
    E = list(zip(G.u.tolist(), G.v.tolist()))
    S = {a * n + b for a, b in E}
    m = len(E)
    alvo = int(trocas_por_aresta * m)
    feitas = tent = 0
    while feitas < alvo and tent < 20 * alvo:
        tent += 1
        i, j = rng.randrange(m), rng.randrange(m)
        if i == j:
            continue
        a, b = E[i]
        c, d = E[j]
        if rng.random() < 0.5:
            c, d = d, c
        if a == d or c == b or a == c or b == d:
            continue
        k1 = a * n + d if a < d else d * n + a
        k2 = c * n + b if c < b else b * n + c
        if k1 in S or k2 in S:
            continue
        S.discard(E[i][0] * n + E[i][1])
        S.discard(E[j][0] * n + E[j][1])
        S.add(k1)
        S.add(k2)
        E[i] = (min(a, d), max(a, d))
        E[j] = (min(c, b), max(c, b))
        feitas += 1
    arr = np.array(E)
    H = de_arestas(arr[:, 0], arr[:, 1], n=n)
    H.trocas = (feitas, tent)
    return H


def barabasi_albert(n: int, m0: int, seed: int = 0):
    """Cada novo vértice liga-se a m0 vértices escolhidos com prob. ∝ grau
    (lista de extremidades repetidas). Começa de um clique de m0+1 vértices."""
    rng = random.Random(seed)
    u, v, rep = [], [], []
    for a in range(m0 + 1):
        for b in range(a + 1, m0 + 1):
            u.append(a); v.append(b); rep += [a, b]
    for x in range(m0 + 1, n):
        alvos = set()
        while len(alvos) < m0:
            alvos.add(rep[rng.randrange(len(rep))])
        for y in alvos:
            u.append(x); v.append(y); rep += [x, y]
    return de_arestas(u, v, n=n)


def watts_strogatz(n: int, k: int, p: float, seed: int = 0):
    """Anel com cada vértice ligado aos k/2 vizinhos de cada lado; cada aresta
    tem a ponta final religada com prob. p (sem laços/duplicatas)."""
    rng = random.Random(seed)
    S = set()
    for x in range(n):
        for j in range(1, k // 2 + 1):
            y = (x + j) % n
            S.add((min(x, y), max(x, y)))
    for (x, y) in sorted(S):
        if rng.random() < p:
            for _ in range(20):
                z = rng.randrange(n)
                e = (min(x, z), max(x, z))
                if z != x and e not in S:
                    S.discard((x, y))
                    S.add(e)
                    break
    arr = np.array(sorted(S))
    return de_arestas(arr[:, 0], arr[:, 1], n=n)


def clusterizacao_reticulado(k_medio: float) -> float:
    """C de um anel regular com grau k: 3(k−2) / (4(k−1))."""
    k = max(2.0, k_medio)
    return 3 * (k - 2) / (4 * (k - 1))


def small_world(C: float, L: float, C_rand: float, L_rand: float, k_medio: float) -> dict:
    C_latt = clusterizacao_reticulado(k_medio)
    gamma = C / C_rand
    lam = L / L_rand
    sigma = gamma / lam
    omega = L_rand / L - C / C_latt
    return {"C": C, "L": L, "C_rand": C_rand, "L_rand": L_rand, "C_latt": C_latt,
            "gamma_C_sobre_Crand": gamma, "lambda_L_sobre_Lrand": lam,
            "sigma": sigma, "omega": omega,
            "criterio_L_proximo": bool(lam <= 1.5), "criterio_C_muito_maior": bool(gamma >= 10),
            "small_world": bool(sigma > 1 and lam <= 1.5 and gamma >= 10)}


def referencias_analiticas_er(n: int, k_medio: float) -> dict:
    return {"C_rand_analitico": k_medio / n, "L_rand_analitico": math.log(n) / math.log(k_medio)}
