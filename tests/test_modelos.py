import numpy as np
from src import modelos as Mo
from src.metricas import clusterizacao


def test_gnm_e_trocas(aleatorio):
    _, G = aleatorio
    E = Mo.erdos_renyi_gnm(500, 2000, seed=1)
    assert E.n == 500 and E.m == 2000
    T = Mo.troca_dupla_arestas(G, 3, seed=1)
    assert np.array_equal(np.sort(T.graus()), np.sort(G.graus())) and np.array_equal(T.graus(), G.graus())
    assert T.m == G.m


def test_ba_ws():
    B = Mo.barabasi_albert(300, 3, seed=1)
    assert B.n == 300 and B.graus().min() >= 3
    W = Mo.watts_strogatz(200, 6, 0.0, seed=1)
    assert W.m == 600
    assert abs(clusterizacao(W.adj)["clusterizacao_media"] - Mo.clusterizacao_reticulado(6)) < 1e-12


def test_small_world_criterios():
    r = Mo.small_world(C=0.4, L=3.0, C_rand=0.004, L_rand=2.8, k_medio=30)
    assert r["small_world"] and r["sigma"] > 1
    r2 = Mo.small_world(C=0.005, L=3.0, C_rand=0.004, L_rand=2.8, k_medio=30)
    assert not r2["small_world"]
