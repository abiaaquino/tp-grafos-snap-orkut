import sys
from pathlib import Path
import networkx as nx
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.tratamento import de_arestas  # noqa: E402


def de_nx(H):
    H = nx.convert_node_labels_to_integers(H)
    e = np.array(list(H.edges())) if H.number_of_edges() else np.zeros((0, 2), int)
    return de_arestas(e[:, 0], e[:, 1], n=H.number_of_nodes())


@pytest.fixture
def aleatorio():
    H = nx.gnm_random_graph(120, 400, seed=7)
    return H, de_nx(H)
