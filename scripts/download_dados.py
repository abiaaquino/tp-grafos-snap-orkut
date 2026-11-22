#!/usr/bin/env python3
"""Baixa os arquivos brutos do com-Orkut do SNAP para data/raw/ e confere o tamanho.

    python scripts/download_dados.py

Os arquivos brutos (~450 MB) não são versionados (ver .gitignore).
"""
from pathlib import Path
import sys
import urllib.request

RAIZ = Path(__file__).resolve().parents[1]
RAW = RAIZ / "data" / "raw"
BASE = "https://snap.stanford.edu/data/bigdata/communities/"
ARQUIVOS = {  # nome -> tamanho esperado em bytes (conferido na cópia usada no TP)
    "com-orkut.ungraph.txt.gz": 447_251_958,
    "com-orkut.top5000.cmty.txt.gz": 3_124_707,
}


def main() -> int:
    RAW.mkdir(parents=True, exist_ok=True)
    ok = True
    for nome, esperado in ARQUIVOS.items():
        destino = RAW / nome
        if destino.exists() and destino.stat().st_size == esperado:
            print(f"[ok] {nome} já existe ({esperado:,} bytes)")
            continue
        print(f"baixando {BASE + nome} ...")
        urllib.request.urlretrieve(BASE + nome, destino)
        tam = destino.stat().st_size
        bate = tam == esperado
        ok &= bate
        print(f"[{'ok' if bate else 'DIVERGE'}] {nome}: {tam:,} bytes (esperado {esperado:,})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
