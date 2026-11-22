import math
import networkx as nx
import numpy as np
from src import robustez as R
from src.medicao import intervalo_confianca, ajuste_proporcional, ajuste_linear_nm
from src.relatorio import fmt_num, valor_macro, tabela_latex
import pandas as pd
from conftest import de_nx


def test_r_e_estrela():
    assert R.r_remocao(9917) == 496
    G = de_nx(nx.star_graph(20))
    m = R.metricas_apos_remocao(G, [0], k_fontes=5)
    assert m["c"] == 20 and m["S"] == 1 / 20


def test_ic():
    x = np.arange(30.0)
    ic = intervalo_confianca(x)
    assert ic["metodo"].startswith("Normal") and abs(ic["critico"] - 1.96) < 0.01
    ic2 = intervalo_confianca(x[:10])
    assert ic2["metodo"].startswith("t-Student") and ic2["critico"] > 2.2
    assert ic2["n"] == 10        # regressão: n é o nº de repetições, não o tamanho do grafo


def test_ajustes():
    n = np.array([1, 2, 3, 4.0]); m = np.array([2, 5, 5, 9.0])
    assert abs(ajuste_proporcional(2 * n, n)["c"] - 2) < 1e-12
    r = ajuste_linear_nm(3 * n + 0.5 * m, n, m)
    assert abs(r["a"] - 3) < 1e-9 and abs(r["b"] - 0.5) < 1e-9


def test_formatacao_br():
    assert fmt_num(1234567) == "1.234.567"
    assert fmt_num(2.5) == "2,50"
    assert fmt_num(True) == "sim"
    assert valor_macro(96.9, 1) == "96{,}9"
    assert valor_macro("10%") == r"10\%"
    t = tabela_latex(pd.DataFrame({"a": ["x_y"], "b": [1.5]}), "cap", "tab:x")
    assert r"x\_y" in t and "1,50" in t
