#!/usr/bin/env python
"""Event displays: per-module visible-energy yield with true impact, reconstructed
position, and recovered cluster energy overlaid.

For each particle (e-, gamma, pi-, mu-) x energy (default 2, 5, 20 GeV) a dedicated
small run (run_id 'disp-...', excluded from resolution analyses) is simulated on demand;
one figure per particle shows rows = energies, cols = events (default 4) -> with the
defaults 4 x 3 x 4 = 48 displayed events. Output: plots/event_displays/.
"""
import argparse
import multiprocessing as mp
import sys
from pathlib import Path

import duckdb
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

sys.path.insert(0, str(Path(__file__).parent))
from g4cal import plotstyle
from g4cal.config import load_config
from g4cal.store import Job, PARTICLE_TAGS, run_job, store_glob


def ensure_run(cfg, preset, particle, energy, n_events, existing_ids):
    rid = f"disp-{preset}-{PARTICLE_TAGS[particle]}-e{energy:g}"
    if rid in existing_ids:
        return rid
    det = cfg.detector.presets[preset]
    job = Job(run_id=rid, preset=preset, nx=det.nx, ny=det.ny, particle=particle,
              e_min=energy, e_max=energy, gun_mode=cfg.production.gun_mode,
              n_events=n_events, seed=910000 + hash(rid) % 10000,
              cal=float(cfg.response.calibration["global"]))
    return job  # caller runs pending jobs in parallel


def fetch(cfg, run_id, n_events):
    ev = duckdb.execute(f"""
        SELECT event, pdg, ekin, x, y, n_clusters, reco_e, reco_x, reco_y
        FROM read_parquet('{store_glob(cfg, "cal_events")}')
        WHERE run_id = ? ORDER BY event LIMIT ?""", [run_id, n_events]).df()
    hits = duckdb.execute(f"""
        SELECT event, ix, iy, e_vis FROM read_parquet('{store_glob(cfg, "cal_hits")}')
        WHERE run_id = ? AND event IN (SELECT event FROM (
            SELECT DISTINCT event FROM read_parquet('{store_glob(cfg, "cal_events")}')
            WHERE run_id = ? ORDER BY event LIMIT ?))""",
        [run_id, run_id, n_events]).df()
    return ev, hits


def draw_panel(ax, cfg, nx, ny, ev_row, ev_hits, vmax):
    pitch = float(cfg.detector.pitch_mm)
    half_x, half_y = 0.5 * nx * pitch, 0.5 * ny * pitch
    grid = np.zeros((ny, nx))
    for _, h in ev_hits.iterrows():
        grid[int(h.iy), int(h.ix)] = h.e_vis
    masked = np.ma.masked_less_equal(grid, 0.0)
    im = ax.imshow(masked, origin="lower", cmap="Blues",
                   norm=LogNorm(vmin=1e-3, vmax=max(vmax, 1e-2)),
                   extent=[-half_x, half_x, -half_y, half_y])
    # module boundaries
    for i in range(nx + 1):
        ax.axvline(-half_x + i * pitch, color="#cccccc", lw=0.4, zorder=2)
    for i in range(ny + 1):
        ax.axhline(-half_y + i * pitch, color="#cccccc", lw=0.4, zorder=2)
    ax.plot(ev_row.x, ev_row.y, "x", color=plotstyle.C_OFF, ms=9, mew=2.2,
            label="true impact", zorder=4)
    if ev_row.n_clusters >= 1 and np.isfinite(ev_row.reco_x):
        ax.plot(ev_row.reco_x, ev_row.reco_y, "o", mfc="none", mec=plotstyle.C_ALT,
                ms=11, mew=2.2, label="reco cluster", zorder=4)
    rec = f"E_rec={ev_row.reco_e:.2f} GeV, {int(ev_row.n_clusters)} cl" \
        if ev_row.n_clusters >= 1 else "no cluster"
    ax.set_title(f"evt {int(ev_row.event)}: {rec}", fontsize=8)
    ax.set_xticks([]); ax.set_yticks([])
    ax.grid(False)
    return im


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preset", default="p10x10")
    ap.add_argument("--particles", nargs="*", default=["e-", "gamma", "pi-", "mu-"])
    ap.add_argument("--energies", nargs="*", type=float, default=[2.0, 5.0, 20.0])
    ap.add_argument("--cols", type=int, default=4, help="events shown per row")
    ap.add_argument("--events", type=int, default=50, help="events simulated per run")
    ap.add_argument("-o", "--override", nargs="*", default=[])
    args = ap.parse_args()

    plotstyle.apply()
    cfg = load_config(args.override)
    det = cfg.detector.presets[args.preset]
    outdir = Path(cfg.paths.plots) / "event_displays"
    outdir.mkdir(parents=True, exist_ok=True)

    try:
        existing = set(duckdb.query(
            f"SELECT run_id FROM read_parquet('{store_glob(cfg, 'runs')}')"
        ).df()["run_id"])
    except Exception:
        existing = set()

    pending = []
    for particle in args.particles:
        for e in args.energies:
            j = ensure_run(cfg, args.preset, particle, e, args.events, existing)
            if isinstance(j, Job):
                pending.append(j)
    if pending:
        print(f"simulating {len(pending)} display runs of {args.events} events")
        with mp.Pool(min(len(pending), int(cfg.production.workers))) as pool:
            pool.starmap(run_job, [(j, cfg) for j in pending])

    for particle in args.particles:
        nrow, ncol = len(args.energies), args.cols
        fig, axes = plt.subplots(nrow, ncol, figsize=(2.6 * ncol + 1.2, 2.7 * nrow),
                                 squeeze=False)
        im = None
        for r, e in enumerate(args.energies):
            rid = f"disp-{args.preset}-{PARTICLE_TAGS[particle]}-e{e:g}"
            ev, hits = fetch(cfg, rid, ncol)
            for c in range(ncol):
                ax = axes[r][c]
                if c >= len(ev):
                    ax.axis("off")
                    continue
                row = ev.iloc[c]
                im = draw_panel(ax, cfg, det.nx, det.ny, row,
                                hits[hits.event == row.event], vmax=e)
                if c == 0:
                    ax.set_ylabel(f"E = {e:g} GeV", fontsize=10)
        h, l = axes[0][0].get_legend_handles_labels()
        fig.legend(h, l, loc="upper right", fontsize=9, ncols=2)
        fig.suptitle(f"{args.preset} {particle}: e_vis yield per module "
                     "(log color), true impact x, reco cluster o", fontsize=11)
        cbar = fig.colorbar(im, ax=[a for row in axes for a in row], shrink=0.85, pad=0.02)
        cbar.set_label("e_vis [GeV]")
        png = plotstyle.save_fig(fig, f"evtdisp_{args.preset}_{PARTICLE_TAGS[particle]}.png",
                                 "phase1/event_displays")
        plt.close(fig)
        print(f"wrote {png}")


if __name__ == "__main__":
    main()
