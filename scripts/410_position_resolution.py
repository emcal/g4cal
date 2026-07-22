#!/usr/bin/env python
"""Position-resolution analysis (companion deliverable).

(1) sigma_x(reco_x - true_x) vs energy, e- on 3x3 and 10x10, with a p/sqrt(E) (+) q fit.
(2) S-shape at 5 GeV: mean and sigma of (reco_x - true_x) vs impact position across the
    central crystal (bias near cell edges is expected).

Outputs: plots/position_resolution_vs_E.png, plots/position_sshape_5gev.png,
         reports/position_fits.json, fit numbers on stdout.
"""
import json
import sys
from pathlib import Path

import duckdb
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).parent))
from g4cal import plotstyle
from g4cal.config import load_config
from g4cal.store import store_glob
from g4cal.fitting import fit_gauss_iterative


def dx_and_local(cfg, preset, particle, energy):
    """(reco_x - true_x, local impact x within the central cell) arrays."""
    det = cfg.detector.presets[preset]
    pitch = float(cfg.detector.pitch_mm)
    cx0 = (det.nx // 2 - 0.5 * (det.nx - 1)) * pitch  # central crystal center, mm
    q = f"""
        SELECT e.reco_x - e.x AS dx, e.x - {cx0} AS xloc
        FROM read_parquet('{store_glob(cfg, "cal_events")}') e
        JOIN read_parquet('{store_glob(cfg, "runs")}') r USING (run_id)
        WHERE r.preset = ? AND r.particle = ? AND r.e_min = r.e_max
          AND abs(r.e_min - ?) < 1e-9 AND r.run_id NOT LIKE 'calib-%' AND r.run_id NOT LIKE 'disp-%'
          AND e.n_clusters >= 1 AND isfinite(e.reco_x)
    """
    df = duckdb.execute(q, [preset, particle, energy]).df()
    return df["dx"].to_numpy(), df["xloc"].to_numpy()


def sigma_model(E, p, q):
    return np.sqrt(p**2 / E + q**2)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--override", nargs="*", default=[],
                    help="OmegaConf dotlist, e.g. analysis.min_events=100")
    plotstyle.apply()
    cfg = load_config(ap.parse_args().override)
    plots, reports = Path(cfg.paths.plots), Path(cfg.paths.reports)
    plots.mkdir(exist_ok=True)
    reports.mkdir(exist_ok=True)
    out = {}

    # (1) sigma_x vs E
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for preset, color, marker in (("p10x10", plotstyle.C_SIM, "o"),
                                  ("p3x3", plotstyle.C_ALT, "^")):
        Es, sigs, sig_errs = [], [], []
        for e in cfg.production.energies_gev:
            dx, _ = dx_and_local(cfg, preset, "e-", float(e))
            if len(dx) < cfg.analysis.min_events:
                continue
            fit = fit_gauss_iterative(dx, cfg.analysis.fit_nsigma,
                                      cfg.analysis.fit_iterations)
            if fit is None:
                continue
            _, sig, _, sig_err = fit
            Es.append(float(e)); sigs.append(sig); sig_errs.append(sig_err)
        if len(Es) < 3:
            print(f"  {preset}: not enough points")
            continue
        popt, pcov = curve_fit(sigma_model, Es, sigs, p0=[3.0, 0.5],
                               sigma=sig_errs, absolute_sigma=True, maxfev=5000)
        perr = np.sqrt(np.diag(pcov))
        print(f"{preset} e-: sigma_x = {popt[0]:.2f}+-{perr[0]:.2f} mm/sqrt(E) "
              f"(+) {popt[1]:.2f}+-{perr[1]:.2f} mm")
        out[f"{preset}_sigma_x_vs_E"] = {
            "p_mm": popt[0], "p_err": perr[0], "q_mm": popt[1], "q_err": perr[1],
            "E": Es, "sigma_mm": list(map(float, sigs)),
            "sigma_err_mm": list(map(float, sig_errs))}
        ax.errorbar(Es, sigs, yerr=sig_errs, fmt=marker, color=color, zorder=3,
                    label=f"{preset[1:]}: ({popt[0]:.2f}±{perr[0]:.2f})/√E ⊕ "
                          f"({popt[1]:.2f}±{perr[1]:.2f}) mm")
        Ec = np.linspace(0.8, 21, 200)
        ax.plot(Ec, sigma_model(Ec, *popt), "-", color=color, alpha=0.85)
    ax.set_xlabel("kinetic energy [GeV]")
    ax.set_ylabel(r"$\sigma_x$(reco $-$ true) [mm]")
    ax.set_title("e-, position resolution vs energy (central-crystal illumination)")
    ax.set_xlim(0, 21)
    ax.set_ylim(0, None)
    ax.legend(fontsize=8.5)
    plotstyle.save_fig(fig, "position_resolution_vs_E.png", "phase1")
    plt.close(fig)

    # (2) S-shape at 5 GeV, 10x10
    dx, xloc = dx_and_local(cfg, "p10x10", "e-", 5.0)
    pitch = float(cfg.detector.pitch_mm)
    edges = np.linspace(-pitch / 2, pitch / 2, 21)
    centers = 0.5 * (edges[1:] + edges[:-1])
    mean_b, sig_b = np.full(len(centers), np.nan), np.full(len(centers), np.nan)
    for i in range(len(centers)):
        sel = dx[(xloc >= edges[i]) & (xloc < edges[i + 1])]
        if len(sel) < 100:
            continue
        fit = fit_gauss_iterative(sel, cfg.analysis.fit_nsigma,
                                  cfg.analysis.fit_iterations, bins=60)
        if fit:
            mean_b[i], sig_b[i] = fit[0], fit[1]
    out["sshape_5gev_10x10"] = {"x_local_mm": list(map(float, centers)),
                                "bias_mm": list(map(float, mean_b)),
                                "sigma_mm": list(map(float, sig_b))}

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.4, 6.2), sharex=True,
                                   gridspec_kw={"hspace": 0.08})
    ax1.axhline(0, color="#bbbbbb", lw=0.8)
    ax1.plot(centers, mean_b, "o-", color=plotstyle.C_SIM)
    ax1.set_ylabel(r"$\langle$reco $-$ true$\rangle_x$ [mm]")
    ax1.set_title("10x10 e- 5 GeV, S-shape across the central crystal")
    ax2.plot(centers, sig_b, "o-", color=plotstyle.C_SIM)
    ax2.set_ylabel(r"$\sigma_x$ [mm]")
    ax2.set_xlabel("impact x relative to crystal center [mm]")
    ax2.set_ylim(0, None)
    plotstyle.save_fig(fig, "position_sshape_5gev.png", "phase1")
    plt.close(fig)

    (reports / "position_fits.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {reports / 'position_fits.json'} and figures in {plots}")


if __name__ == "__main__":
    main()
