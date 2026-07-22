#!/usr/bin/env python
"""CI sanity check for the analysis path, no Geant4 needed.

Feeds synthetic per-energy Gaussian samples (widths drawn from a known
a/sqrt(E) (+) b/E (+) c model) through the same fitting the resolution analysis
uses, asserts the model is recovered, and writes one plot artifact. Run:

    python scripts/ci_check.py
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from g4cal import plotstyle
from g4cal.fitting import fit_gauss_iterative, fit_resolution, resolution_model

TRUE = (0.030, 0.010, 0.007)          # a, b, c
ENERGIES = [1, 2, 3, 5, 8, 12, 16, 20]
RNG = np.random.default_rng(1)


def main():
    E, res, err = [], [], []
    for e in ENERGIES:
        sigma_rel = resolution_model(e, *TRUE)
        sample = RNG.normal(e, sigma_rel * e, 20000)   # reco energy around truth
        fit = fit_gauss_iterative(sample)
        assert fit is not None, f"Gaussian fit failed at E={e}"
        mu, sig, mu_err, sig_err = fit
        E.append(e); res.append(sig / mu); err.append(sig_err / mu)

    popt, perr, chi2, ndf = fit_resolution(np.array(E), np.array(res), np.array(err))
    a, b, c = popt
    print(f"recovered a={a:.4f} b={b:.4f} c={c:.4f}  (true {TRUE})")
    # loose tolerances: this checks the pipeline wiring, not statistics
    assert abs(a - TRUE[0]) < 0.006, f"a off: {a}"
    assert abs(c - TRUE[2]) < 0.006, f"c off: {c}"

    plotstyle.apply()
    fig, ax = plt.subplots(figsize=(5, 3.5))
    xs = np.linspace(1, 20, 100)
    ax.plot(xs, 100 * resolution_model(xs, a, b, c), color=plotstyle.C_SIM,
            label=f"fit: {a*100:.2f}%/√E ⊕ {b*100:.2f}%/E ⊕ {c*100:.2f}%")
    ax.errorbar(E, 100 * np.array(res), yerr=100 * np.array(err), fmt="o",
                color=plotstyle.C_SIM, label="synthetic points")
    ax.set_xlabel("E [GeV]"); ax.set_ylabel(r"$\sigma_E/E$ [%]")
    ax.set_title("CI self-check: resolution-fit round-trip")
    ax.legend()
    plotstyle.save_fig(fig, "ci_selfcheck.png", "ci")
    print("OK: analysis path recovers the model; artifact written to plots/ci/")


if __name__ == "__main__":
    sys.exit(main())
