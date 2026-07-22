"""Shared matplotlib style for calorimeter figures.

Palette: Okabe-Ito subset, CVD-validated (dataviz validator: all checks pass on
light surface). Identity is never color-alone: every figure uses distinct markers
and a legend.

save_fig() writes into PLOTS_ROOT, which defaults to <g4cal>/plots but is
overridable via the G4CAL_PLOTS_ROOT env var (cal-fpga points it at its own
plots/ so its phase-folder convention keeps working through the shim).
"""
import os
import shutil
from pathlib import Path

import matplotlib as mpl

PLOTS_ROOT = Path(os.environ.get(
    "G4CAL_PLOTS_ROOT", Path(__file__).resolve().parent.parent / "plots"))


def save_fig(fig, name, phase, current_name=None):
    """Save a figure twice: into its phase archive (plots/<phase>/, kept forever so
    progress/regress across iterations stays traceable) and into plots/current/
    (always the latest canonical version, overwritten). Iteration-tagged figures
    pass an untagged current_name so plots/current/ holds exactly one live copy.
    Returns the archive path."""
    arch = PLOTS_ROOT / phase
    cur = PLOTS_ROOT / "current"
    arch.mkdir(parents=True, exist_ok=True)
    cur.mkdir(parents=True, exist_ok=True)
    fig.savefig(arch / name)
    shutil.copyfile(arch / name, cur / (current_name or name))
    return arch / name

C_SIM = "#0072B2"      # simulation, full response chain ON
C_OFF = "#D55E00"      # chain (b)+(c) OFF (reco on edep_true)
C_ALT = "#009E73"      # second preset in overlays
C_REF = "#555555"      # harvested reference curves

def apply():
    mpl.rcParams.update({
        "figure.dpi": 130,
        "font.size": 10,
        "axes.grid": True,
        "grid.color": "#dddddd",
        "grid.linewidth": 0.6,
        "axes.axisbelow": True,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "errorbar.capsize": 2.5,
        "lines.linewidth": 1.6,
        "lines.markersize": 5,
        "savefig.bbox": "tight",
    })
