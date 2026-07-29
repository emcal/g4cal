#!/usr/bin/env python
"""Energy-resolution analysis (primary validation deliverable).

Per (preset, particle): Gaussian-fit reconstructed energy at each fixed energy point,
fit sigma_E/E = a/sqrt(E) (+) b/E (+) c (quadrature), overlay the harvested measured
model. Done twice: full response chain ON (reco_e) and chain (b)+(c) OFF (reco on
edep_true, same events) which should land near the intrinsic Geant part.

Outputs: plots/energy_resolution_<preset>_<particle>.png,
         plots/energy_resolution_3x3_vs_10x10_e-.png (leakage overlay),
         reports/resolution_fits.json, fit numbers on stdout.
"""
import json
import sys
from pathlib import Path

import duckdb
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from g4cal import plotstyle
from g4cal.config import load_config
from g4cal.store import store_glob
from g4cal.fitting import fit_gauss_iterative, fit_resolution, resolution_model, res_terms


def energy_points(cfg, preset, particle, energy, column):
    q = f"""
        SELECT e.{column} AS v
        FROM read_parquet('{store_glob(cfg, "cal_events")}') e
        JOIN read_parquet('{store_glob(cfg, "runs")}') r USING (run_id)
        WHERE r.preset = ? AND r.particle = ? AND r.e_min = r.e_max
          AND abs(r.e_min - ?) < 1e-9 AND r.run_id NOT LIKE 'calib-%' AND r.run_id NOT LIKE 'disp-%'
          AND e.n_clusters >= 1
    """
    return duckdb.execute(q, [preset, particle, energy]).df()["v"].to_numpy()


def scan(cfg, preset, particle, column):
    """Fit each energy point; returns dict with arrays E, res, res_err, mu, sig."""
    out = {"E": [], "res": [], "res_err": [], "mu": [], "sig": [], "n": []}
    for e in cfg.production.energies_gev:
        v = energy_points(cfg, preset, particle, float(e), column)
        if len(v) < cfg.analysis.min_events:
            print(f"  skip E={e}: only {len(v)} events")
            continue
        fit = fit_gauss_iterative(v, cfg.analysis.fit_nsigma, cfg.analysis.fit_iterations)
        if fit is None:
            print(f"  skip E={e}: fit failed")
            continue
        mu, sig, mu_err, sig_err = fit
        res = sig / mu
        res_err = res * np.hypot(sig_err / sig, mu_err / mu)
        for k, val in zip(("E", "res", "res_err", "mu", "sig", "n"),
                          (float(e), res, res_err, mu, sig, len(v))):
            out[k].append(val)
    return out


def fmt_abc(popt, perr, chi2, ndf):
    def term(v, e, unit):
        if v < 1e-6 and e > 1.0:   # parameter pinned at the zero bound, error meaningless
            return f"0 (at bound){unit}"
        return f"({100*v:.2f} ± {100*e:.2f})%{unit}"
    return (f"a = {term(popt[0], perr[0], '·√GeV')}\n"
            f"b = {term(popt[1], perr[1], '·GeV')}\n"
            f"c = {term(popt[2], perr[2], '')}\n"
            f"χ²/ndf = {chi2:.1f}/{ndf}")


def one_figure(cfg, preset, particle, results, out_png):
    meas = res_terms(cfg.response.smearing.measured_res)
    intr = res_terms(cfg.response.smearing.intrinsic_res)
    Ec = np.linspace(0.8, 21, 300)

    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for key, color, marker, label in (
            ("on", plotstyle.C_SIM, "o", "sim, full response chain"),
            ("off", plotstyle.C_OFF, "s", "sim, light-collection+smear OFF")):
        r, fitres = results[key]
        ax.errorbar(r["E"], 100 * np.array(r["res"]), yerr=100 * np.array(r["res_err"]),
                    fmt=marker, color=color, label=label, zorder=3)
        if fitres:
            popt, perr, chi2, ndf = fitres
            ax.plot(Ec, 100 * resolution_model(Ec, *popt), "-", color=color, alpha=0.85)
    ax.plot(Ec, 100 * resolution_model(Ec, *meas), "--", color=plotstyle.C_REF,
            label="measured PbWO4 model (halld)")
    ax.plot(Ec, 100 * resolution_model(Ec, *intr), ":", color=plotstyle.C_REF,
            label="intrinsic Geant part (halld)")

    popt, perr, chi2, ndf = results["on"][1]
    ax.text(0.97, 0.95, "fit (chain ON):\n" + fmt_abc(popt, perr, chi2, ndf),
            transform=ax.transAxes, ha="right", va="top", fontsize=8.5,
            bbox=dict(boxstyle="round", fc="white", ec="#cccccc"))
    ax.set_xlabel("kinetic energy [GeV]")
    ax.set_ylabel(r"$\sigma_E/E$ [%]")
    ax.set_title(f"{preset} {particle}, central-crystal illumination")
    ax.set_xlim(0, 21)
    ax.set_ylim(0, None)
    ax.legend(fontsize=8.5, loc="lower left")
    plotstyle.save_fig(fig, out_png, "phase1")
    plt.close(fig)


def overlay_figure(cfg, res3, res10, out_png):
    Ec = np.linspace(0.8, 21, 300)
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for r, fitres, color, marker, label in (
            (res10[0], res10[1], plotstyle.C_SIM, "o", "10x10 (reference conditions)"),
            (res3[0], res3[1], plotstyle.C_ALT, "^", "3x3 (beam-test size)")):
        ax.errorbar(r["E"], 100 * np.array(r["res"]), yerr=100 * np.array(r["res_err"]),
                    fmt=marker, color=color, label=label, zorder=3)
        popt, perr, chi2, ndf = fitres
        ax.plot(Ec, 100 * resolution_model(Ec, *popt), "-", color=color, alpha=0.85)
        ax.text(0.97, 0.95 if color == plotstyle.C_SIM else 0.68,
                label.split(" (")[0] + ":\n" + fmt_abc(popt, perr, chi2, ndf),
                transform=ax.transAxes, ha="right", va="top", fontsize=8,
                color=color, bbox=dict(boxstyle="round", fc="white", ec="#cccccc"))
    ax.plot(Ec, 100 * resolution_model(Ec, *res_terms(cfg.response.smearing.measured_res)), "--",
            color=plotstyle.C_REF, label="measured PbWO4 model (halld)")
    ax.set_xlabel("kinetic energy [GeV]")
    ax.set_ylabel(r"$\sigma_E/E$ [%]")
    ax.set_title("e-, transverse-leakage effect: 3x3 vs 10x10 (chain ON)")
    ax.set_xlim(0, 21)
    ax.set_ylim(0, None)
    ax.legend(fontsize=8.5, loc="lower left")
    plotstyle.save_fig(fig, out_png, "phase1")
    plt.close(fig)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--override", nargs="*", default=[],
                    help="OmegaConf dotlist, e.g. analysis.min_events=100")
    plotstyle.apply()
    cfg = load_config(ap.parse_args().override)
    plots = Path(cfg.paths.plots)
    reports = Path(cfg.paths.reports)
    plots.mkdir(exist_ok=True)
    reports.mkdir(exist_ok=True)

    all_fits = {}
    keep = {}
    for preset in cfg.production.presets:
        for particle in cfg.production.particles:
            print(f"== {preset} {particle}")
            results = {}
            for key, column in (("on", "reco_e"), ("off", "reco_e_true")):
                r = scan(cfg, preset, particle, column)
                if len(r["E"]) < 4:
                    print(f"  not enough points for {key}")
                    results[key] = (r, None)
                    continue
                popt, perr, chi2, ndf = fit_resolution(r["E"], r["res"], r["res_err"])
                results[key] = (r, (popt, perr, chi2, ndf))
                print(f"  chain {key:3s}: a={100*popt[0]:.2f}+-{100*perr[0]:.2f} "
                      f"b={100*popt[1]:.2f}+-{100*perr[1]:.2f} "
                      f"c={100*popt[2]:.2f}+-{100*perr[2]:.2f} [%] chi2/ndf={chi2:.1f}/{ndf}")
                all_fits[f"{preset}_{particle}_{key}"] = {
                    "a": popt[0], "b": popt[1], "c": popt[2],
                    "a_err": perr[0], "b_err": perr[1], "c_err": perr[2],
                    "chi2": chi2, "ndf": ndf,
                    "points": {k: list(map(float, v)) for k, v in r.items()},
                }
            if results["on"][1]:
                one_figure(cfg, preset, particle, results,
                           f"energy_resolution_{preset}_{particle}.png")
            keep[(preset, particle)] = results["on"]

    if ("p3x3", "e-") in keep and ("p10x10", "e-") in keep \
            and keep[("p3x3", "e-")][1] and keep[("p10x10", "e-")][1]:
        overlay_figure(cfg, keep[("p3x3", "e-")], keep[("p10x10", "e-")],
                       "energy_resolution_3x3_vs_10x10_e-.png")

    all_fits["reference"] = {
        "measured": res_terms(cfg.response.smearing.measured_res),
        "intrinsic": res_terms(cfg.response.smearing.intrinsic_res),
    }
    (reports / "resolution_fits.json").write_text(json.dumps(all_fits, indent=1))
    print(f"wrote {reports / 'resolution_fits.json'} and figures in {plots}")


if __name__ == "__main__":
    main()
