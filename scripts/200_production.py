#!/usr/bin/env python
"""Run a production: matrix of (preset x particle x energy point) split into chunks,
executed as parallel single-threaded calsim+calreco jobs.

Examples:
    python 200_production.py                       # full validation matrix from config
    python 200_production.py --smoke               # tiny smoke matrix (store round-trip)
    python 200_production.py --presets p3x3 --particles e- --energies 5 --events 1000
"""
import argparse
import multiprocessing as mp
import sys
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from g4cal.config import load_config
from g4cal.store import Job, make_run_id, run_job

CFG = None  # set per worker via initializer (picklable OmegaConf is fine, keep global)


def _init(cfg):
    global CFG
    CFG = cfg


def _run(job):
    try:
        run_job(job, CFG)
        return (job.run_id, "ok")
    except Exception:
        return (job.run_id, traceback.format_exc(limit=3))


def build_jobs(cfg, presets, particles, energies, events, chunk_size, gun_mode):
    cal = float(cfg.response.calibration["global"])  # "global" is a python keyword
    jobs, seed = [], int(cfg.production.seed_base)
    for preset in presets:
        det = cfg.detector.presets[preset]
        for particle in particles:
            for e in energies:
                nleft, chunk = events, 0
                while nleft > 0:
                    n = min(chunk_size, nleft)
                    rid = make_run_id(preset, particle, e, e, chunk)
                    jobs.append(Job(run_id=rid, preset=preset, nx=det.nx, ny=det.ny,
                                    particle=particle, e_min=e, e_max=e,
                                    gun_mode=gun_mode, n_events=n, seed=seed, cal=cal))
                    seed += 1
                    nleft -= n
                    chunk += 1
    return jobs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--presets", nargs="*")
    ap.add_argument("--particles", nargs="*")
    ap.add_argument("--energies", nargs="*", type=float)
    ap.add_argument("--events", type=int)
    ap.add_argument("--workers", type=int)
    ap.add_argument("--gun-mode")
    ap.add_argument("--smoke", action="store_true",
                    help="tiny matrix: 3x3+100x100, e-/pi-/mu-, 5 GeV, few events")
    args = ap.parse_args()

    cfg = load_config()
    prod = cfg.production

    if args.smoke:
        jobs = (build_jobs(cfg, ["p3x3"], ["e-", "pi-", "mu-"], [5.0], 500, 500,
                           prod.gun_mode)
                + build_jobs(cfg, ["p100x100"], ["e-"], [5.0], 100, 100, prod.gun_mode))
    else:
        jobs = build_jobs(
            cfg,
            args.presets or list(prod.presets),
            args.particles or list(prod.particles),
            args.energies or list(prod.energies_gev),
            args.events or int(prod.events_per_point),
            int(prod.chunk_size),
            args.gun_mode or prod.gun_mode,
        )

    workers = args.workers or int(prod.workers)
    print(f"{len(jobs)} jobs on {workers} workers", flush=True)
    failed = 0
    with mp.Pool(workers, initializer=_init, initargs=(cfg,)) as pool:
        for i, (rid, status) in enumerate(pool.imap_unordered(_run, jobs), 1):
            if status != "ok":
                failed += 1
                print(f"[{i}/{len(jobs)}] FAILED {rid}\n{status}", flush=True)
            elif i % 10 == 0 or i == len(jobs):
                print(f"[{i}/{len(jobs)}] done (last: {rid})", flush=True)
    print(f"production finished: {len(jobs) - failed} ok, {failed} failed", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
