"""Parte III — verificação de lei de potência P(k) ~ C·k^(−γ).

Três abordagens (todas próprias):
  1. regressão linear em log P(k) × log k: (a) pontos brutos; (b) bins logarítmicos;
  2. regressão na CCDF P(K ≥ k): inclinação = −(γ − 1);
  3. máxima verossimilhança discreta (Clauset, Shalizi & Newman 2009) com
     normalização exata pela zeta de Hurwitz ζ(γ, k_min) e k_min escolhido por
     mínima distância de Kolmogorov–Smirnov; razão de verossimilhança (Vuong)
     contra exponencial discreta na mesma cauda; comparação global com Poisson.
Entrada sempre como (valores k, contagens) — funciona para 10 mil ou 3 milhões
de vértices.
"""
from __future__ import annotations

import numpy as np
from scipy import optimize, special, stats


def _tabela(graus=None, k=None, c=None):
    if graus is not None:
        g = np.asarray(graus)
        k, c = np.unique(g[g > 0], return_counts=True)
    k, c = np.asarray(k, float), np.asarray(c, float)
    ok = k > 0
    return k[ok], c[ok]


def regressao_loglog(k, c) -> dict:
    p = c / c.sum()
    x, y = np.log10(k), np.log10(p)
    r = stats.linregress(x, y)
    return {"gamma": -r.slope, "intercepto": r.intercept, "r2": r.rvalue ** 2, "pontos": len(k)}


def regressao_bins_log(k, c, n_bins: int = 25) -> dict:
    bordas = np.unique(np.floor(np.logspace(0, np.log10(k.max() + 1), n_bins + 1)))
    total = c.sum()
    xs, ys = [], []
    for a, b in zip(bordas[:-1], bordas[1:]):
        sel = (k >= a) & (k < b)
        if c[sel].sum() > 0:
            xs.append(np.sqrt(a * (b - 1)) if b - 1 > a else a)
            ys.append(c[sel].sum() / total / (b - a))
    xs, ys = np.array(xs), np.array(ys)
    r = stats.linregress(np.log10(xs), np.log10(ys))
    return {"gamma": -r.slope, "intercepto": r.intercept, "r2": r.rvalue ** 2,
            "pontos": len(xs), "x": xs, "y": ys}


def ccdf(k, c):
    p = c / c.sum()
    return 1.0 - np.concatenate([[0.0], np.cumsum(p)[:-1]])


def regressao_ccdf(k, c) -> dict:
    cc = ccdf(k, c)
    r = stats.linregress(np.log10(k), np.log10(cc))
    return {"inclinacao": r.slope, "gamma": 1 - r.slope, "r2": r.rvalue ** 2}


def _loglik_pl(alpha, k, c, kmin):
    return -c.sum() * np.log(special.zeta(alpha, kmin)) - alpha * (c * np.log(k)).sum()


def mle_discreta(k, c, kmin) -> dict:
    sel = k >= kmin
    kt, ct = k[sel], c[sel]
    res = optimize.minimize_scalar(lambda a: -_loglik_pl(a, kt, ct, kmin),
                                   bounds=(1.01, 20.0), method="bounded")
    alpha = float(res.x)
    # KS entre CCDF empírica e do modelo, na cauda
    emp = ccdf(kt, ct)
    mod = special.zeta(alpha, kt) / special.zeta(alpha, kmin)
    ks = float(np.max(np.abs(emp - mod)))
    ntail = ct.sum()
    erro = (alpha - 1) / np.sqrt(ntail)
    return {"gamma": alpha, "kmin": float(kmin), "ks": ks, "n_cauda": int(ntail), "erro_padrao": float(erro),
            "loglik": float(_loglik_pl(alpha, kt, ct, kmin))}


def escolher_kmin(k, c, min_cauda: int = 50, max_candidatos: int = 120) -> dict:
    cauda = np.cumsum(c[::-1])[::-1]
    cand = k[cauda >= min_cauda]
    if len(cand) > max_candidatos:
        cand = np.unique(cand[np.unique(np.geomspace(1, len(cand), max_candidatos).astype(int) - 1)])
    melhores = [mle_discreta(k, c, km) for km in cand]
    best = min(melhores, key=lambda r: r["ks"])
    best["candidatos_testados"] = len(cand)
    return best


def comparar_exponencial(k, c, kmin, alpha) -> dict:
    """Vuong: log-razão por observação entre lei de potência e exponencial
    discreta P(k) = (1 − e^{−λ}) e^{−λ(k − kmin)} na mesma cauda."""
    sel = k >= kmin
    kt, ct = k[sel], c[sel]
    lam = np.log(1 + 1 / ((ct * kt).sum() / ct.sum() - kmin))
    lp = -np.log(special.zeta(alpha, kmin)) - alpha * np.log(kt)
    le = np.log(1 - np.exp(-lam)) - lam * (kt - kmin)
    d = lp - le
    n = ct.sum()
    R = float((ct * d).sum())
    media = R / n
    dp = np.sqrt((ct * (d - media) ** 2).sum() / n)
    z = R / (np.sqrt(n) * dp) if dp > 0 else np.inf
    p = float(2 * stats.norm.sf(abs(z)))
    return {"lambda_exp": float(lam), "R": R, "z_vuong": float(z), "p_valor": p,
            "favorece": "lei de potência" if R > 0 else "exponencial", "significativo": p < 0.05}


def comparar_poisson(k, c) -> dict:
    """Poisson com λ = ⟨k⟩ (o que um grafo ER produziria): KS global e razão
    de verossimilhança contra a lei de potência ajustada com k_min = 1."""
    n = c.sum()
    lam = (k * c).sum() / n
    emp = ccdf(k, c)
    mod = stats.poisson.sf(k - 1, lam)
    ll_pois = float((c * stats.poisson.logpmf(k, lam)).sum())
    pl = mle_discreta(k, c, 1)
    return {"lambda": float(lam), "ks_poisson": float(np.max(np.abs(emp - mod))),
            "loglik_poisson": ll_pois, "loglik_lei_potencia_kmin1": pl["loglik"],
            "gamma_kmin1": pl["gamma"], "favorece": "lei de potência" if pl["loglik"] > ll_pois else "Poisson"}


def analisar(graus=None, k=None, c=None, min_cauda: int = 50) -> dict:
    k, c = _tabela(graus, k, c)
    bruto = regressao_loglog(k, c)
    bins = regressao_bins_log(k, c)
    cc = regressao_ccdf(k, c)
    mle = escolher_kmin(k, c, min_cauda=min_cauda)
    exp = comparar_exponencial(k, c, mle["kmin"], mle["gamma"])
    poi = comparar_poisson(k, c)
    return {"regressao_bruta": bruto, "regressao_bins_log": {x: y for x, y in bins.items() if x not in ("x", "y")},
            "regressao_ccdf": cc, "mle": mle, "vs_exponencial": exp, "vs_poisson": poi,
            "_bins": (bins["x"], bins["y"]), "_k": k, "_c": c}
