"""Estilo común de las figuras.

Paleta validada para daltonismo. Con la variable de entorno REPORTE=1 las figuras salen
sin títulos internos y con fuente Latin Modern, para el informe en LaTeX
(se guardan en figuras/reporte/).
"""
import os
import sys
from pathlib import Path

import matplotlib as mpl
from matplotlib import font_manager

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))
REPORTE = os.environ.get("REPORTE") == "1"
FIG = RAIZ / "figuras" / ("reporte" if REPORTE else "")
RES = RAIZ / "resultados"
FIG.mkdir(parents=True, exist_ok=True)
RES.mkdir(exist_ok=True)

COL = {
    "tierra": "#2563A8",
    "jupiter": "#B7791F",
    "saturno": "#6D5BB0",
    "sonda": "#C2410C",
    "alterno": "#0891B2",
    "tinta": "#1F2A37",
    "gris": "#8A94A3",
    "rejilla": "#E4E8EE",
}

for f in Path("/usr/share/fonts/truetype/google-fonts").glob("Poppins-*.ttf"):
    font_manager.fontManager.addfont(str(f))

mpl.rcParams.update({
    "font.family": ["Poppins", "DejaVu Sans"] if any(f.name == "Poppins" for f in font_manager.fontManager.ttflist)
                   else ["DejaVu Sans"],
    "font.size": 11,
    "axes.titlesize": 14,
    "axes.titleweight": "bold",
    "axes.edgecolor": COL["gris"],
    "axes.labelcolor": COL["tinta"],
    "xtick.color": COL["tinta"],
    "ytick.color": COL["tinta"],
    "axes.grid": True,
    "grid.color": COL["rejilla"],
    "grid.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "figure.dpi": 110,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
    "mathtext.fontset": "dejavusans",
})

if REPORTE:
    for f in Path("/usr/share/texmf/fonts/opentype/public/lm").glob("lmroman10-*.otf"):
        font_manager.fontManager.addfont(str(f))
    mpl.rcParams.update({
        "font.family": ["Latin Modern Roman", "DejaVu Serif"],
        "font.size": 13, "axes.titlesize": 13, "axes.titleweight": "normal",
        "axes.labelsize": 13, "mathtext.fontset": "cm",
    })


def figura(*a, figsize=None, **k):
    """plt.subplots; en modo reporte reduce el tamaño para que el texto quede legible."""
    if REPORTE and figsize is not None:
        figsize = (figsize[0] * 0.72, figsize[1] * 0.72)
    return plt.subplots(*a, figsize=figsize, **k)


def titulo(obj, texto, **kw):
    """Título de la figura; en el informe se omite (lo reemplaza el caption de LaTeX)."""
    if REPORTE:
        return
    (obj.suptitle if hasattr(obj, "suptitle") else obj.set_title)(texto, **kw)


def nota(obj, x, y, texto, **kw):
    """Nota dentro de la figura; en el informe se omite."""
    if not REPORTE:
        obj.text(x, y, texto, **kw)
