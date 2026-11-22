import random
import networkx as nx
import numpy as np
from src import algoritmos as A
from src.tratamento import com_pesos
from conftest import de_nx


def test_bfs_dfs(aleatorio):
    H, G = aleatorio
    d, pai = A.bfs(G.adj, 0)
    ref = nx.single_source_shortest_path_length(H, 0)
    assert all(d[v] == ref.get(v, -1) for v in range(G.n))
    ordem, _ = A.dfs(G.adj, 0)
    assert set(ordem) == nx.node_connected_component(H, 0)


def test_componentes_desconexo():
    G = de_nx(nx.disjoint_union(nx.path_graph(5), nx.cycle_graph(3)))
    rot, tam = A.componentes_conexas(G.adj)
    assert sorted(tam) == [3, 5]


def test_euler():
    assert A.verificar_euler(de_nx(nx.cycle_graph(6)).adj)["euleriano"]
    assert A.verificar_euler(de_nx(nx.path_graph(4)).adj)["semi_euleriano"]
    assert not A.verificar_euler(de_nx(nx.star_graph(3)).adj)["semi_euleriano"]
    G = de_nx(nx.complete_graph(5))
    c = A.circuito_euleriano(G.adj)
    assert len(c) == G.m + 1 and c[0] == c[-1]


def _ponderado(seed=3):
    H = nx.gnm_random_graph(60, 200, seed=seed)
    rnd = random.Random(seed)
    for a, b in H.edges():
        H[a][b]["weight"] = rnd.uniform(0.1, 5)
    G = de_nx(H)
    w = np.array([H[a][b]["weight"] for a, b in zip(G.u, G.v)])
    return H, com_pesos(G, w)


def test_dijkstra_bellman_floyd_consistentes():
    H, G = _ponderado()
    d, _ = A.dijkstra(G.adj_ponderada(), 0)
    ref = nx.single_source_dijkstra_path_length(H, 0)
    assert all(abs(d[v] - ref[v]) < 1e-9 for v in ref)
    bf, rod, neg = A.bellman_ford(G.n, A.arcos_nao_direcionados(G.u, G.v, G.w), 0)
    assert not neg and rod < G.n
    assert all(abs(bf[v] - d[v]) < 1e-9 for v in range(G.n) if d[v] < A.INF)
    D = A.matriz_pesos(G.n, G.u, G.v, G.w)
    F1 = A.floyd_warshall_numpy(D)
    F2 = np.array(A.floyd_warshall_python(D))
    assert np.allclose(F1, F2)
    assert np.allclose(F1[0][np.isfinite(F1[0])], [d[v] for v in range(G.n) if d[v] < A.INF])


def test_bellman_ford_ciclo_negativo():
    arcos = [(0, 1, 1.0), (1, 2, -2.0), (2, 0, 0.5), (2, 3, 1.0)]
    _, _, neg = A.bellman_ford(4, arcos, 0)
    assert neg


def test_tarjan(aleatorio):
    H = nx.barbell_graph(5, 3)
    H.add_edge(0, 100)
    G = de_nx(H)
    pontes, art = A.tarjan_pontes_articulacoes(G.adj)
    Hn = nx.convert_node_labels_to_integers(H)
    assert {tuple(sorted(e)) for e in nx.bridges(Hn)} == set(pontes)
    assert set(nx.articulation_points(Hn)) == art
    H2, G2 = aleatorio
    p2, a2 = A.tarjan_pontes_articulacoes(G2.adj)
    assert set(nx.articulation_points(H2)) == a2
    assert {tuple(sorted(e)) for e in nx.bridges(H2)} == set(p2)


def test_prim_kruskal():
    H, G = _ponderado(5)
    H = H.subgraph(max(nx.connected_components(H), key=len)).copy()
    G = de_nx(H)
    w = np.array([H[a][b]["weight"] for a, b in zip(*[np.array(sorted(H.nodes()))[x] for x in (G.u, G.v)])])
    _, tk = A.kruskal(G.n, G.u.tolist(), G.v.tolist(), w.tolist())
    _, tp = A.prim(com_pesos(G, w).adj_ponderada(), 0)
    ref = nx.minimum_spanning_tree(H).size(weight="weight")
    assert abs(tk - tp) < 1e-9 and abs(tk - ref) < 1e-9
