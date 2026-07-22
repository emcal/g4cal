#!/usr/bin/env python
"""Derive the single global calibration constant (response stage d).

Runs the calibration run from config (default: 5 GeV e- on the central crystal of the
10x10) with cal=1, Gaussian-fits the reconstructed energy peak, and writes
configs/calib.yaml with global = E_true / mu so reconstructed energy is unbiased.
"""
import multiprocessing as mp
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import duckdb

from g4cal.config import load_config, calib_path
from g4cal.store import Job, make_run_id, run_job, store_glob
from g4cal.fitting import fit_gauss_iterative


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--override", nargs="*", default=[],
                    help="OmegaConf dotlist, e.g. calibration_run.n_events=1000")
    cfg = load_config(ap.parse_args().override)
    cr = cfg.calibration_run
    det = cfg.detector.presets[cr.preset]
    chunk = int(cfg.production.chunk_size)
    seed0 = int(cfg.production.seed_base) - 1000  # outside the production seed range

    jobs = []
    nleft, k = int(cr.n_events), 0
    while nleft > 0:
        n = min(chunk, nleft)
        rid = "calib-" + make_run_id(cr.preset, cr.particle, cr.energy_gev, cr.energy_gev, k)
        jobs.append(Job(run_id=rid, preset=cr.preset, nx=det.nx, ny=det.ny,
                        particle=cr.particle, e_min=cr.energy_gev, e_max=cr.energy_gev,
                        gun_mode=cfg.production.gun_mode, n_events=n,
                        seed=seed0 + k, cal=1.0))
        nleft -= n
        k += 1

    print(f"calibration: {len(jobs)} jobs", flush=True)
    with mp.Pool(min(len(jobs), int(cfg.production.workers))) as pool:
        pool.starmap(run_job, [(j, cfg) for j in jobs])

    q = f"""
        SELECT reco_e FROM read_parquet('{store_glob(cfg, "cal_events")}')
        WHERE run_id LIKE 'calib-%' AND n_clusters >= 1
    """
    vals = duckdb.query(q).df()["reco_e"].to_numpy()
    fit = fit_gauss_iterative(vals, cfg.analysis.fit_nsigma, cfg.analysis.fit_iterations)
    if fit is None:
        print("calibration fit FAILED")
        return 1
    mu, sig, mu_err, _ = fit
    cal = float(cr.energy_gev) / mu
    print(f"events={len(vals)} mu={mu:.4f}+-{mu_err:.4f} GeV sigma={sig:.4f} "
          f"-> global calibration = {cal:.5f}", flush=True)

    calib_path().write_text(
        "# written by 150_calibrate.py: global calibration from "
        f"{cr.energy_gev} GeV {cr.particle} on {cr.preset} central crystal\n"
        f"# fit: mu={mu:.5f}+-{mu_err:.5f} GeV, sigma={sig:.5f}, n={len(vals)}\n"
        "response:\n  calibration:\n"
        f"    global: {cal:.6f}\n")
    print(f"wrote {calib_path()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
