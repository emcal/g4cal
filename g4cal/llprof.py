"""Log-likelihood profile tables: per-node histograms of f = e_vis / E_incident.

The LL method (ilreco PHASE II) tabulates, for every node (i, j) = cell-center
distance from the photon impact point in 0.01-cell units (folded to
i >= j >= 0 by square symmetry, i <= i_max = 500), the distribution of the
cell energy fraction f. Zeros (cell below sparsification) are counted
separately as n_zero; nonzero f is histogrammed on log-spaced bins so the
rapidly-falling noncentral distributions keep shape resolution at small f.

Two subcommands, designed for farm production (simulation-pipeline
90-create-emcal-profile-jobs.py):

  fill   one calsim point-gun run (impact at sub-cell offset (fx, fy) in
         0.01-cell units from the central-cell center) -> one sparse partial:
         only the nodes that run covers.
           python -m g4cal.llprof fill --hits hits.csv --run-json run.json \
             --out partial.npz --energy 5 --fx 13 --fy 7

  merge  sum partials of one energy -> one dense table (the profile file):
           python -m g4cal.llprof merge --out profile_e5.npz partials/*.npz

Profile file (npz), triangular node order nid = i*(i+1)/2 + j:
  counts   uint32 [n_nodes, n_bins]  nonzero-f histogram counts
  n_zero   int64  [n_nodes]          zero-cell counts
  n_slots  int64  [n_nodes]          total samples (counts.sum(1) + n_zero)
  edges    float64 [n_bins+1]        log-spaced f bin edges
  energy, i_max, n_events            scalars
"""
import argparse
import json

import numpy as np
import pandas as pd


def bin_edges(f_min, f_max, bins_per_decade):
    n = int(round(np.log10(f_max / f_min) * bins_per_decade))
    return np.logspace(np.log10(f_min), np.log10(f_max), n + 1)


def grid_node_ids(nx, fx, fy):
    """Node ids for every cell of the nx*nx grid at impact offset (fx, fy).

    Returns (nid[nx, nx], max_offset[nx, nx]); the gun sits at the center of
    cell (nx//2, nx//2) plus (fx, fy)/100 cells.
    """
    c = nx // 2
    xoff = np.abs((np.arange(nx) - c) * 100 - fx)
    yoff = np.abs((np.arange(nx) - c) * 100 - fy)
    a = np.maximum(xoff[:, None], yoff[None, :])
    b = np.minimum(xoff[:, None], yoff[None, :])
    return a * (a + 1) // 2 + b, a


def nid_to_ij(nid):
    i = ((np.sqrt(8 * nid.astype(np.float64) + 1) - 1) // 2).astype(np.int64)
    return i, nid - i * (i + 1) // 2


def fill(a):
    edges = bin_edges(a.f_min, a.f_max, a.bins_per_decade)
    nb = len(edges) - 1
    n_ev = json.load(open(a.run_json))["n_processed"]

    nid_grid, max_off = grid_node_ids(a.nx, a.fx, a.fy)
    in_range = max_off <= a.i_max
    uniq, mult = np.unique(nid_grid[in_range], return_counts=True)

    hits = pd.read_csv(a.hits, usecols=["ix", "iy", "e_vis"],
                       dtype={"ix": np.int32, "iy": np.int32, "e_vis": np.float64})
    f = hits.e_vis.to_numpy() / a.energy
    hix, hiy = hits.ix.to_numpy(), hits.iy.to_numpy()
    sel = (f > 0) & in_range[hix, hiy]
    row = np.searchsorted(uniq, nid_grid[hix[sel], hiy[sel]])
    fbin = np.clip(np.searchsorted(edges, f[sel], side="right") - 1, 0, nb - 1)
    counts = np.bincount(row * nb + fbin,
                         minlength=len(uniq) * nb).reshape(len(uniq), nb)
    n_zero = n_ev * mult - counts.sum(axis=1)

    np.savez_compressed(a.out, nid=uniq.astype(np.int64),
                        counts=counts.astype(np.uint32),
                        n_zero=n_zero.astype(np.int64),
                        mult=mult.astype(np.int64),
                        edges=edges, n_events=n_ev, energy=float(a.energy),
                        fx=a.fx, fy=a.fy, i_max=a.i_max)
    print(f"llprof fill: {a.out}  n_events={n_ev} nodes={len(uniq)} "
          f"nonzero_entries={counts.sum()} zeros={n_zero.sum()}")


def merge(a):
    ref = np.load(a.partials[0])
    edges, energy, i_max = ref["edges"], float(ref["energy"]), int(ref["i_max"])
    nb = len(edges) - 1
    n_nodes = (i_max + 1) * (i_max + 2) // 2
    counts = np.zeros((n_nodes, nb), np.int64)
    n_zero = np.zeros(n_nodes, np.int64)
    n_slots = np.zeros(n_nodes, np.int64)
    n_events = 0
    for path in a.partials:
        d = np.load(path)
        if len(d["edges"]) != len(edges) or float(d["energy"]) != energy:
            raise ValueError(f"{path}: binning/energy mismatch with {a.partials[0]}")
        nid, n_ev = d["nid"], int(d["n_events"])
        counts[nid] += d["counts"]
        n_zero[nid] += d["n_zero"]
        n_slots[nid] += n_ev * d["mult"]
        n_events += n_ev
    np.savez_compressed(a.out, counts=counts.astype(np.uint32), n_zero=n_zero,
                        n_slots=n_slots, edges=edges, energy=energy,
                        i_max=i_max, n_events=n_events)
    covered = int((n_slots > 0).sum())
    print(f"llprof merge: {a.out}  partials={len(a.partials)} events={n_events} "
          f"nodes covered {covered}/{n_nodes}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    pf = sub.add_parser("fill", help="one calsim run -> sparse partial npz")
    pf.add_argument("--hits", required=True, help="calsim hits.csv")
    pf.add_argument("--run-json", required=True, help="calsim run.json (event count)")
    pf.add_argument("--out", required=True)
    pf.add_argument("--energy", type=float, required=True, help="incident E, GeV")
    pf.add_argument("--fx", type=int, required=True, help="impact x offset, 0.01-cell units")
    pf.add_argument("--fy", type=int, required=True)
    pf.add_argument("--nx", type=int, default=20)
    pf.add_argument("--i-max", type=int, default=500)
    pf.add_argument("--f-min", type=float, default=1e-6)
    pf.add_argument("--f-max", type=float, default=1.5)
    pf.add_argument("--bins-per-decade", type=int, default=100)
    pf.set_defaults(func=fill)

    pm = sub.add_parser("merge", help="sum partials -> dense profile npz")
    pm.add_argument("--out", required=True)
    pm.add_argument("partials", nargs="+")
    pm.set_defaults(func=merge)

    a = ap.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
