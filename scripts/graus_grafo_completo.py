#!/usr/bin/env python3
"""Distribuição de graus do com-Orkut COMPLETO (3.072.441 vértices), em streaming.

Gera data/processed/orkut_completo_distribuicao_graus.csv (k, n_vertices), arquivo
pequeno que permite testar a lei de potência no grafo original sem versionar o bruto.
Uso: python scripts/graus_grafo_completo.py  (~75 s; exige data/raw/com-orkut.ungraph.txt.gz)
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parents[1]
BRUTO = RAIZ / "data" / "raw" / "com-orkut.ungraph.txt.gz"
SAIDA = RAIZ / "data" / "processed" / "orkut_completo_distribuicao_graus.csv"


def main() -> int:
    grau = np.zeros(3_100_000, dtype=np.int64)
    for df in pd.read_csv(BRUTO, sep=r"\s+", comment="#", header=None, names=["u", "v"],
                          dtype=np.int64, chunksize=5_000_000):
        for col in ("u", "v"):
            x = df[col].to_numpy()
            if x.max() >= len(grau):
                grau = np.concatenate([grau, np.zeros(x.max() + 1 - len(grau), np.int64)])
            grau += np.bincount(x, minlength=len(grau))
    grau = grau[grau > 0]
    k, c = np.unique(grau, return_counts=True)
    pd.DataFrame({"k": k, "n_vertices": c}).to_csv(SAIDA, index=False)
    print(f"vértices={len(grau):,} arestas={grau.sum() // 2:,} k_max={k.max():,} -> {SAIDA.name}")
    return 0 if len(grau) == 3_072_441 else 2


if __name__ == "__main__":
    sys.exit(main())
