#!/usr/bin/env python3
"""Confere o dataset tratado em data/processed contra o resumo do tratamento.
Retorna código 0 se tudo bate e 1 caso contrário."""
import json
import sys
from pathlib import Path

import pandas as pd

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ))
from src.algoritmos import componentes_conexas  # noqa: E402
from src.tratamento import PROC, carregar_lcc  # noqa: E402


def main() -> int:
    resumo = json.loads((PROC / "resumo_tratamento.json").read_text(encoding="utf-8"))
    ar = pd.read_csv(PROC / "orkut_g27_lcc_arestas.csv")
    nos = pd.read_csv(PROC / "orkut_g27_lcc_nos.csv")
    G, _ = carregar_lcc()
    _, tam = componentes_conexas(G.adj)
    aud = resumo["auditoria_original"]
    checks = {
        "original |V| = 3.072.441 (SNAP)": aud["vertices_distintos"] == 3_072_441,
        "original |E| = 117.185.083 (SNAP)": aud["linhas_de_aresta"] == 117_185_083,
        "|V| final = resumo": len(nos) == resumo["final_lcc"]["vertices"] == G.n,
        "|E| final = resumo": len(ar) == resumo["final_lcc"]["arestas"] == G.m,
        "sem auto-laços": bool((ar["src"] != ar["dst"]).all()),
        "sem arestas duplicadas": not ar.assign(a=ar.min(axis=1), b=ar.max(axis=1))[["a", "b"]].duplicated().any(),
        "todas as pontas estão em nos.csv": bool(set(ar["src"]) | set(ar["dst"]) == set(nos["node_id"])),
        "grafo conexo (1 componente)": len(tam) == 1,
        "grau no subgrafo <= grau original": bool((nos["grau_subgrafo"] <= nos["grau_original"]).all()),
    }
    dist = PROC / "orkut_completo_distribuicao_graus.csv"
    if dist.exists():
        d = pd.read_csv(dist)
        checks["distribuição completa soma 3.072.441 vértices"] = int(d["n_vertices"].sum()) == 3_072_441
        checks["distribuição completa soma 2|E|"] = int((d["k"] * d["n_vertices"]).sum()) == 2 * 117_185_083
    for nome, ok in checks.items():
        print(f"[{'OK' if ok else 'FALHA'}] {nome}")
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
