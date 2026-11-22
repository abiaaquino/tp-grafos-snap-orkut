import networkx as nx
import numpy as np
from src.centralidade import pagerank, brandes


def test_pagerank(aleatorio):
    H, G = aleatorio
    pr, _ = pagerank(G)
    ref = nx.pagerank(H, alpha=0.85, tol=1e-12)
    assert np.allclose(pr, [ref[i] for i in range(G.n)], atol=1e-6)


def test_brandes_exato(aleatorio):
    H, G = aleatorio
    b = brandes(G.adj)
    ref = nx.betweenness_centrality(H, normalized=True)
    assert np.allclose(b, [ref[i] for i in range(G.n)], atol=1e-9)
