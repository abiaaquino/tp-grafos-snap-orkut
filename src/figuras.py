"""Gravação padronizada de figuras (PNG 300 dpi + PDF vetorial)."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .relatorio import FIG  # noqa: E402

CORES = {"azul": "#2a6fb0", "laranja": "#e07b28", "verde": "#3a9a5b", "vermelho": "#c23b3b",
         "roxo": "#7a4fa3", "cinza": "#6b6b6b"}

plt.rcParams.update({"figure.dpi": 100, "savefig.dpi": 300, "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.alpha": 0.3, "axes.formatter.use_locale": False})


def salvar(fig, nome: str):
    FIG.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / f"{nome}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{nome}.pdf", bbox_inches="tight")
    plt.close(fig)
    return FIG / f"{nome}.png"
