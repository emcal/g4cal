"""Log-likelihood profile tables: histograms of f = e_vis / E_incident per cell offset.

The LL method (ilreco PHASE II) needs, for a photon that hit the calorimeter at a known
point, the probability distribution of the energy fraction f = e_vis / E_incident that a
cell collects, as a function of where that cell sits relative to the impact point.

Table key. The displacement from the impact point to a cell center is measured in
0.01-cell units: (offset_x, offset_y). Crystals are square, so the distribution is
unchanged by flipping the sign of either component or swapping the axes; each
displacement is therefore stored once, as

    offset_major = max(|offset_x|, |offset_y|)
    offset_minor = min(|offset_x|, |offset_y|)

One table entry = one (offset_major, offset_minor) pair = one histogram of f. Entries
form a triangle 0 <= offset_minor <= offset_major <= max_cell_offset, flattened as

    entry_id = offset_major * (offset_major + 1) / 2 + offset_minor

(Larin's text profile prof[i][j] uses the same triangle with i = offset_major,
j = offset_minor.)

Cells below sparsification produce no hit row, so f = 0 samples are counted separately
as n_zero; nonzero f is histogrammed on log-spaced bins so the rapidly-falling
distributions of far cells keep shape resolution at small f.

Two subcommands, designed for farm production (simulation-pipeline
90-create-emcal-profile-jobs.py runs calsim + `fill` inside each job; `merge`
combines the partials after download):

  fill   one calsim point-gun run (hits.csv or hits.csv.gz) -> one sparse partial npz
         holding only the entries that run covers:
           python -m g4cal.llprof fill --hits hits.csv --run-json run.json \\
             --out partials/<tag>.npz --energy 5 --impact-x 13 --impact-y 7

  merge  sum the partials of one energy -> one dense table (the profile file):
           python -m g4cal.llprof merge --out profile_e5.npz partials/*.npz

Profile file (npz), rows in entry_id order:
  counts            uint32 [n_entries, n_bins]  nonzero-f histogram counts
  n_zero            int64  [n_entries]          f = 0 sample counts
  n_slots           int64  [n_entries]          total samples (counts.sum(1) + n_zero)
  edges             float64 [n_bins+1]          log-spaced f bin edges
  energy, max_cell_offset, n_events             scalars

Both commands write their npz atomically (tmp file + rename), so a file that exists is
always complete - killed jobs never leave half-written deliverables.
"""
import argparse
import json
import os

import numpy as np
import pandas as pd


def bin_edges(f_min, f_max, bins_per_decade):
    n_bins = int(round(np.log10(f_max / f_min) * bins_per_decade))
    return np.logspace(np.log10(f_min), np.log10(f_max), n_bins + 1)


def grid_entry_ids(crystals_nx, impact_x, impact_y):
    """Table entry id for every cell of the crystals_nx * crystals_nx grid, for a gun at
    the center of cell (crystals_nx//2, crystals_nx//2) plus (impact_x, impact_y)/100 cells.

    Returns (entry_id[crystals_nx, crystals_nx], offset_major[same]); offset_major is
    returned separately so callers can cut on max_cell_offset.
    """
    center = crystals_nx // 2
    offset_x = np.abs((np.arange(crystals_nx) - center) * 100 - impact_x)   # per column, 0.01-cell units
    offset_y = np.abs((np.arange(crystals_nx) - center) * 100 - impact_y)   # per row
    offset_major = np.maximum(offset_x[:, None], offset_y[None, :])
    offset_minor = np.minimum(offset_x[:, None], offset_y[None, :])
    return offset_major * (offset_major + 1) // 2 + offset_minor, offset_major


def entry_id_to_offsets(entry_id):
    """Inverse of the triangular flattening: entry_id -> (offset_major, offset_minor)."""
    offset_major = ((np.sqrt(8 * entry_id.astype(np.float64) + 1) - 1) // 2).astype(np.int64)
    return offset_major, entry_id - offset_major * (offset_major + 1) // 2


def save_npz_atomic(path, **arrays):
    """np.savez_compressed via tmp file + rename: an existing file is always complete."""
    tmp_path = path + ".part.npz"   # ends in .npz so numpy does not append another
    np.savez_compressed(tmp_path, **arrays)
    os.replace(tmp_path, path)


def fill(args):
    edges = bin_edges(args.f_min, args.f_max, args.bins_per_decade)
    n_bins = len(edges) - 1
    n_events = json.load(open(args.run_json))["n_processed"]

    entry_id_grid, offset_major_grid = grid_entry_ids(args.crystals_nx, args.impact_x, args.impact_y)
    in_range = offset_major_grid <= args.max_cell_offset
    # Several grid cells can map to the same table entry (mirror positions around the
    # impact point). cells_per_entry counts them - needed for the implicit-zero count.
    entry_ids, cells_per_entry = np.unique(entry_id_grid[in_range], return_counts=True)

    hits = pd.read_csv(args.hits, usecols=["ix", "iy", "e_vis"],
                       dtype={"ix": np.int32, "iy": np.int32, "e_vis": np.float64})
    f = hits.e_vis.to_numpy() / args.energy
    hit_ix, hit_iy = hits.ix.to_numpy(), hits.iy.to_numpy()
    keep = (f > 0) & in_range[hit_ix, hit_iy]
    row = np.searchsorted(entry_ids, entry_id_grid[hit_ix[keep], hit_iy[keep]])
    f_bin = np.clip(np.searchsorted(edges, f[keep], side="right") - 1, 0, n_bins - 1)
    counts = np.bincount(row * n_bins + f_bin, minlength=len(entry_ids) * n_bins)
    counts = counts.reshape(len(entry_ids), n_bins)
    # hits.csv only has cells above sparsification: every (event, cell) slot not counted
    # above is an f = 0 sample.
    n_zero = n_events * cells_per_entry - counts.sum(axis=1)

    save_npz_atomic(args.out,
                    entry_id=entry_ids.astype(np.int64),
                    counts=counts.astype(np.uint32),
                    n_zero=n_zero.astype(np.int64),
                    cells_per_entry=cells_per_entry.astype(np.int64),
                    edges=edges, n_events=n_events, energy=float(args.energy),
                    impact_x=args.impact_x, impact_y=args.impact_y,
                    max_cell_offset=args.max_cell_offset)
    print(f"llprof fill: {args.out}  n_events={n_events} entries={len(entry_ids)} "
          f"nonzero_samples={counts.sum()} zeros={n_zero.sum()}")


def merge(args):
    ref = np.load(args.partials[0])
    edges = ref["edges"]
    energy = float(ref["energy"])
    max_cell_offset = int(ref["max_cell_offset"])
    n_bins = len(edges) - 1
    n_entries = (max_cell_offset + 1) * (max_cell_offset + 2) // 2

    counts = np.zeros((n_entries, n_bins), np.int64)
    n_zero = np.zeros(n_entries, np.int64)
    n_slots = np.zeros(n_entries, np.int64)
    n_events = 0
    for path in args.partials:
        partial = np.load(path)
        if len(partial["edges"]) != len(edges) or float(partial["energy"]) != energy:
            raise ValueError(f"{path}: binning/energy mismatch with {args.partials[0]}")
        entry_id, partial_events = partial["entry_id"], int(partial["n_events"])
        counts[entry_id] += partial["counts"]
        n_zero[entry_id] += partial["n_zero"]
        n_slots[entry_id] += partial_events * partial["cells_per_entry"]
        n_events += partial_events

    save_npz_atomic(args.out, counts=counts.astype(np.uint32), n_zero=n_zero,
                    n_slots=n_slots, edges=edges, energy=energy,
                    max_cell_offset=max_cell_offset, n_events=n_events)
    covered = int((n_slots > 0).sum())
    print(f"llprof merge: {args.out}  partials={len(args.partials)} events={n_events} "
          f"entries covered {covered}/{n_entries}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    pf = sub.add_parser("fill", help="one calsim run -> sparse partial npz")
    pf.add_argument("--hits", required=True, help="calsim hits.csv (or .csv.gz)")
    pf.add_argument("--run-json", required=True, help="calsim run.json (event count)")
    pf.add_argument("--out", required=True)
    pf.add_argument("--energy", type=float, required=True, help="incident E, GeV")
    pf.add_argument("--impact-x", type=int, required=True,
                    help="gun offset from central-cell center, 0.01-cell units")
    pf.add_argument("--impact-y", type=int, required=True)
    pf.add_argument("--crystals-nx", type=int, default=20, help="grid size in crystals (nx * nx)")
    pf.add_argument("--max-cell-offset", type=int, default=500,
                    help="tabulate only cells within this offset of the impact point "
                         "(larger-axis distance, 0.01-cell units)")
    pf.add_argument("--f-min", type=float, default=1e-6)
    pf.add_argument("--f-max", type=float, default=1.5)
    pf.add_argument("--bins-per-decade", type=int, default=100)
    pf.set_defaults(func=fill)

    pm = sub.add_parser("merge", help="sum partials -> dense profile npz")
    pm.add_argument("--out", required=True)
    pm.add_argument("partials", nargs="+")
    pm.set_defaults(func=merge)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
