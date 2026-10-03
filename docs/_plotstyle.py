"""Plot style shared by the pages of the pyscal documentation.

This module is not part of pyscal. The documentation pages import it in
hidden cells, so that every figure uses the same colours and layout.
"""
import os

from pychromatic import Multiplot, Palette

DOCS = os.path.dirname(os.path.abspath(__file__))
EXAMPLES = os.path.join(DOCS, "..", "examples")

DARK = "#363636"
_p = Palette("tableau10")

# one colour per structure, the same in every figure
COLOURS = {
    "fcc": _p.blue.hex,
    "bcc": _p.orange.hex,
    "hcp": _p.green.hex,
    "ico": _p.purple.hex,
    "liquid": _p.grey.hex,
    "others": _p.grey.hex,
    "diamond": _p.red.hex,
    "hex. diamond": _p.teal.hex,
}
# one colour and marker per neighbor method
METHODS = {
    "fixed cutoff": (_p.red.hex, "o"),
    "adaptive": (_p.teal.hex, "s"),
    "SANN": (_p.brown.hex, "^"),
    "Voronoi": (_p.pink.hex, "D"),
}
BLUE, ORANGE, GREEN, RED, PURPLE, GREY, TEAL = (
    _p.blue.hex, _p.orange.hex, _p.green.hex, _p.red.hex, _p.purple.hex,
    _p.grey.hex, _p.teal.hex)


def figure(columns=1, rows=1, width=500, ratio=0.45, **kwargs):
    """A Multiplot of the width used throughout the documentation."""
    return Multiplot(width=width, ratio=ratio, columns=columns, rows=rows, **kwargs)


def label(ax, text):
    """Panel label at the top left, as in "(a)  q6"."""
    ax.set_title(text, loc="left", fontsize=11)


def note(ax, text, loc="upper left"):
    """A short note inside a panel."""
    x, ha = (0.03, "left") if "left" in loc else (0.97, "right")
    y, va = (0.96, "top") if "upper" in loc else (0.04, "bottom")
    ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va, fontsize=9, color=DARK)
