import networkx as nx
import numpy as np
from src import metricas as M
from src.distancias import todos_para_todos, validar_com_bfs_propria
from conftest import de_nx


def test_triangulos_clusterizacao(aleatorio):
    H, G = aleatorio
    c = M.clusterizacao(G.adj)
    assert c["triangulos"] == sum(nx.triangles(H).values()) // 3
    assert abs(c["clusterizacao_media"] - nx.average_clustering(H)) < 1e-12
    assert abs(c["transitividade"] - nx.transitivity(H)) < 1e-12


def test_distribuicao_sem_grau_zero():
    k, pk, cc = M.distribuicao_graus(np.array([0, 1, 1, 2, 3]))
    assert 0 not in k and abs(pk.sum() - 1) < 1e-12 and cc[0] == 1.0


def test_distancias():
    H = nx.connected_watts_strogatz_graph(200, 6, 0.1, seed=1)
    G = de_nx(H)
    r = todos_para_todos(G, bloco=37)
    assert r["diametro"] == nx.diameter(H) and r["raio"] == nx.radius(H)
    assert abs(r["caminho_medio"] - nx.average_shortest_path_length(H)) < 1e-9
    assert validar_com_bfs_propria(G, [0, 5, 50])


def test_densidade():
    assert M.densidade(4, 6) == 1.0
