"""One look for every report figure: colours, fonts, line weights, and a save helper.

Colour follows the design, the same in every figure, so a reader who learns "blue is the
slot design" once can carry it through the report. The three design hues were checked as a
set for colour-blind separation (dataviz palette validator: all adjacent pairs pass). The
conventional / standard-stack reference is grey, because it is what we compare against.
Aqua sits below 3:1 contrast on white, so aqua lines always get a direct label.
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

FIG_DIR = Path(__file__).resolve().parents[2] / "figures"

COLOR = {
    "standard": "#898781",   # grey: the conventional reference
    "undercut": "#1baf7a",   # aqua
    "slot": "#2a78d6",       # blue: the headline design
    "vacuum": "#eb6834",     # orange
    "dual": "#eda100",       # yellow: the dual-use (steer + lens) layout, used rarely
}
INK = "#0b0b0b"        # primary text
INK2 = "#52514e"       # secondary text, direct labels
MUTED = "#898781"      # axis ticks and reference marks
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
HEAT = "#e34948"       # heat-flow arrows in the schematics
BLUE_RAMP = ["#b7d3f6", "#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]  # ordered series

# Materials for the cross-section schematics
MAT = {"si": "#b9b7ae", "oxide": "#dfe9f3", "core": "#3a3a38", "heater": "#d9a521", "gas": "#ffffff"}

FULL = 6.4   # figure width for a full text column, inches
HALF = 3.15  # half column


def apply():
    """Set matplotlib defaults for the report figures."""
    mpl.rcParams.update({
        "font.family": ["Segoe UI", "DejaVu Sans"],
        "mathtext.fontset": "dejavusans",
        "font.size": 8.5,
        "axes.titlesize": 9,
        "axes.labelsize": 8.5,
        "axes.labelcolor": INK2,
        "axes.edgecolor": AXIS,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "axes.axisbelow": True,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelcolor": INK2,
        "ytick.labelcolor": INK2,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "lines.linewidth": 1.6,
        "lines.markersize": 4.5,
        "legend.frameon": False,
        "legend.fontsize": 8,
        "legend.labelcolor": INK2,
        "figure.dpi": 110,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
        "pdf.fonttype": 42,          # embed TrueType so the PDF text stays editable in LaTeX
    })


def panel_label(ax, letter, x=-0.02, y=1.02):
    """(a), (b), ... in the top-left corner, outside the plot."""
    ax.text(x, y, f"({letter})", transform=ax.transAxes, ha="right", va="bottom",
            fontsize=9.5, fontweight="bold", color=INK)


def save(fig, name):
    """Write figures/<name>.pdf (for the LaTeX report) and .png (for slides; 300 dpi, so a
    full-width figure is about 1900 px wide, enough for a 1080p slide)."""
    FIG_DIR.mkdir(exist_ok=True)
    fig.savefig(FIG_DIR / f"{name}.pdf")
    fig.savefig(FIG_DIR / f"{name}.png", dpi=300)
    plt.close(fig)
    print(f"  wrote figures/{name}.pdf and .png")
