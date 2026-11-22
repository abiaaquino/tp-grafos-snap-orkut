#!/usr/bin/env python3
"""Executa todas as análises (Partes I–III) sem abrir o notebook.

    python scripts/run_all.py            # completo (~8–12 min)
    python scripts/run_all.py --rapido   # menos repetições, só para testar o fluxo
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import analise  # noqa: E402

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rapido", action="store_true", help="poucas repetições (não usar para o relatório)")
    analise.executar_tudo(rapido=ap.parse_args().rapido)
