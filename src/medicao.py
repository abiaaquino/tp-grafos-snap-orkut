"""Medição de tempo com intervalos de confiança.

* aquecimento descartado; ``gc.disable()`` durante cada execução;
* IC de 95%: Normal (z = 1,96) quando n >= 30; t-Student (n-1 g.l.) quando
  n < 30 e σ desconhecido;  IC = x̄ ± crítico · s/√n.
"""
from __future__ import annotations

import gc
import os
import platform
import time

import numpy as np
from scipy import stats


def cronometrar(fn, reps: int = 30, aquecimento: int = 1) -> np.ndarray:
    for _ in range(aquecimento):
        fn()
    tempos = []
    for _ in range(reps):
        gc.collect()
        gc.disable()
        try:
            t0 = time.perf_counter()
            fn()
            tempos.append(time.perf_counter() - t0)
        finally:
            gc.enable()
    return np.array(tempos)


def intervalo_confianca(x, conf: float = 0.95) -> dict:
    x = np.asarray(x, dtype=float)
    n = len(x)
    media = float(x.mean())
    dp = float(x.std(ddof=1)) if n > 1 else 0.0
    if n >= 30:
        crit, metodo = float(stats.norm.ppf(0.5 + conf / 2)), "Normal (z)"
    else:
        crit, metodo = float(stats.t.ppf(0.5 + conf / 2, df=n - 1)), f"t-Student ({n - 1} g.l.)"
    meia = crit * dp / np.sqrt(n) if n > 1 else 0.0
    return {"n": n, "media": media, "dp": dp, "metodo": metodo, "critico": crit,
            "ic_inf": media - meia, "ic_sup": media + meia, "meia_largura": meia}


def r2(y, yhat) -> float:
    y, yhat = np.asarray(y, float), np.asarray(yhat, float)
    ss = ((y - y.mean()) ** 2).sum()
    return float(1 - ((y - yhat) ** 2).sum() / ss) if ss > 0 else float("nan")


def ajuste_proporcional(t, f) -> dict:
    """Mínimos quadrados de t ≈ c · f (sem intercepto)."""
    t, f = np.asarray(t, float), np.asarray(f, float)
    c = float((f * t).sum() / (f * f).sum())
    return {"c": c, "r2": r2(t, c * f)}


def ajuste_linear_nm(t, n, m) -> dict:
    """t ≈ a·n + b·m (para algoritmos O(n+m); a razão m/n muda entre subgrafos)."""
    A = np.column_stack([np.asarray(n, float), np.asarray(m, float)])
    coef, *_ = np.linalg.lstsq(A, np.asarray(t, float), rcond=None)
    return {"a": float(coef[0]), "b": float(coef[1]), "r2": r2(t, A @ coef)}


def info_maquina() -> dict:
    cpu = platform.processor() or ""
    try:
        with open("/proc/cpuinfo") as fh:
            for l in fh:
                if l.startswith("model name"):
                    cpu = l.split(":", 1)[1].strip()
                    break
    except OSError:
        pass
    return {"python": platform.python_version(), "so": platform.platform(), "cpu": cpu,
            "nucleos_logicos": os.cpu_count()}
