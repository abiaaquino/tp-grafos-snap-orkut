import numpy as np
from src import leipotencia as LP


def test_mle_recupera_gamma():
    rng = np.random.default_rng(0)
    x = np.floor(rng.pareto(1.5, 20000) + 1).astype(int) * 1  # contínua aproximada, gamma ~2.5
    r = LP.analisar(graus=x)
    assert 2.2 < r["mle"]["gamma"] < 2.8
    assert r["vs_exponencial"]["favorece"] == "lei de potência"


def test_poisson_nao_e_lei_de_potencia():
    rng = np.random.default_rng(1)
    x = rng.poisson(30, 5000)
    r = LP.analisar(graus=x)
    assert r["vs_poisson"]["favorece"] == "Poisson"
