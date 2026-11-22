"""Formatação pt-BR, tabelas LaTeX, macros (valores.tex) e gravação de saídas."""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parents[1]
RES = RAIZ / "resultados"
FIG = RAIZ / "figuras"


def fmt_num(x, casas: int = 2) -> str:
    """1234567.891 -> '1.234.567,89'; inteiros sem casas; bool -> sim/não."""
    if isinstance(x, (bool, np.bool_)):
        return "sim" if x else "não"
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "--"
    if isinstance(x, (int, np.integer)):
        return f"{int(x):,}".replace(",", ".")
    if isinstance(x, (float, np.floating)):
        x = float(x)
        if x != 0 and (abs(x) < 10 ** (-casas) or abs(x) >= 1e9):
            m, e = f"{x:.{casas}e}".split("e")
            return f"{m.replace('.', ',')}e{int(e)}"
        s = f"{x:,.{casas}f}"
        return s.replace(",", "X").replace(".", ",").replace("X", ".")
    return str(x)


_ESC = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#",
        "_": r"\_", "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}", "^": r"\^{}"}


def esc(s) -> str:
    return "".join(_ESC.get(ch, ch) for ch in str(s))


def tabela_latex(df: pd.DataFrame, caption: str, label: str, casas: int | dict = 2,
                 alinhamento: str | None = None) -> str:
    cols = list(df.columns)
    alinhamento = alinhamento or ("l" + "r" * (len(cols) - 1))
    linhas = [r"\begin{table}[ht]", r"\centering", rf"\caption{{{esc(caption)}}}",
              rf"\label{{{label}}}", r"\small", rf"\begin{{tabular}}{{{alinhamento}}}", r"\hline",
              " & ".join(esc(c) for c in cols) + r" \\", r"\hline"]
    for _, row in df.iterrows():
        cel = []
        for c in cols:
            cs = casas.get(c, 2) if isinstance(casas, dict) else casas
            v = row[c]
            cel.append(esc(fmt_num(v, cs)) if not isinstance(v, str) else esc(v))
        linhas.append(" & ".join(cel) + r" \\")
    linhas += [r"\hline", r"\end{tabular}", r"\end{table}"]
    return "\n".join(linhas) + "\n"


def salvar_tabela(df: pd.DataFrame, nome: str, caption: str, label: str, casas=2, alinhamento=None):
    RES.mkdir(exist_ok=True)
    df.to_csv(RES / f"{nome}.csv", index=False)
    (RES / f"{nome}.tex").write_text(tabela_latex(df, caption, label, casas, alinhamento), encoding="utf-8")


def _nome_macro(nome: str) -> str:
    if not re.fullmatch(r"[A-Za-z]+", nome):
        raise ValueError(f"nome de macro inválido (só letras): {nome}")
    return nome


def valor_macro(v, casas: int = 2) -> str:
    s = esc(fmt_num(v, casas)) if not isinstance(v, str) else esc(v)
    return s.replace(",", "{,}")


def escrever_macros(valores: dict, caminho: Path | None = None) -> Path:
    """valores = {nome: valor | (valor, casas)} -> \\newcommand{\\nome}{valor}."""
    caminho = caminho or RES / "valores.tex"
    linhas = ["% Gerado automaticamente por src/relatorio.py — não editar à mão"]
    for nome, v in valores.items():
        casas = 2
        if isinstance(v, tuple):
            v, casas = v
        linhas.append(rf"\newcommand{{\{_nome_macro(nome)}}}{{{valor_macro(v, casas)}}}")
    caminho.write_text("\n".join(linhas) + "\n", encoding="utf-8")
    return caminho


class _Enc(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            return float(o)
        if isinstance(o, np.bool_):
            return bool(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        return super().default(o)


def salvar_json(obj, nome: str) -> Path:
    (RES / "json").mkdir(parents=True, exist_ok=True)
    p = RES / "json" / f"{nome}.json"
    p.write_text(json.dumps(obj, cls=_Enc, indent=2, ensure_ascii=False), encoding="utf-8")
    return p


def ler_json(nome: str):
    return json.loads((RES / "json" / f"{nome}.json").read_text(encoding="utf-8"))
