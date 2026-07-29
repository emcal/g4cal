# g4cal

Standalone Geant4 simulation of a PbWO4 crystal calorimeter, plus a parquet data
store and an ilreco-based reconstruction. Built to reproduce the measured HallD
ECAL (GlueX PbWO4 Compton calorimeter) energy resolution, and to generate training
and validation data for downstream reconstruction work (cal-fpga).

Two independent halves:

- **calsim** — the Geant4 application (C++). Parametric `n_x × n_y` grid of PbWO4
  crystals, a configurable response chain, CSV output. No dependency on ilreco.
- **calreco** — island-clustering reconstruction (C++), a thin adapter over the
  [ilreco](https://github.com/emcal/ilreco) library's context API.

The Python package `g4cal` drives productions (sim → store → reco), holds the
config, and produces the analysis plots.

## Layout

```
src/calsim/           Geant4 application (one executable: calsim)
src/calreco/          reconstruction driver (calreco_NxN executables)
src/ilreco_adapter/   C adapter over ilreco's context API
src/tests/            adapter round-trip test (ctest: adapter_3/10/100)
g4cal/                Python package: config, store, plotstyle, fitting
scripts/              build + production + analysis scripts (see below)
configs/              config.yaml (+ calib.yaml, written by 150_calibrate.py)
data/prof_pwo.dat     PbWO4 shower profile for ilreco
```

## Build

Needs Geant4 (with its datasets) and a C/C++ compiler. ilreco is resolved
automatically — `find_package(ilreco)` if installed, otherwise fetched from git.

```sh
cmake -B build -DCMAKE_PREFIX_PATH=/path/to/geant4
cmake --build build -j
ctest --test-dir build            # adapter round-trip, 3 tests
```

Produces `build/calsim`, `build/calreco_{3x3,10x10,100x100}`, and the test binaries.
The `calreco_NxN` binaries are identical — the array size is a runtime `--nx/--ny`
argument; the names exist only so production scripts have stable paths.

### ilreco resolution

`calreco` needs ilreco's **context API** (`ilreco_config_create` / `ilreco_reconstruct`),
present on github `emcal/ilreco` `master`. CMake tries, in order:

1. `find_package(ilreco)` — an installed package (`ilreco::ilreco`).
2. FetchContent from `G4CAL_ILRECO_GIT`@`G4CAL_ILRECO_TAG`. The default git is a
   sibling `../ilreco` checkout if present (offline-friendly), else the github URL;
   the default tag is `master`.

Override, e.g.:

```sh
cmake -B build -DG4CAL_ILRECO_GIT=https://github.com/emcal/ilreco -DG4CAL_ILRECO_TAG=master
```

To build only the simulation (no ilreco): `-DG4CAL_BUILD_RECO=OFF`.

## Python

```sh
uv venv .venv && uv pip install --python .venv -e .
```

Then source the Geant4 environment (`source <geant4>/bin/geant4.sh`) before running
anything that calls `calsim`.

## Detector and response

One parametric geometry: an `n_x × n_y` array of PbWO4 crystals, 20.55 mm square ×
200 mm long, 0.175 mm Tedlar wrap per lateral face (20.9 mm pitch), from the HallD
ECAL / hdds geometry. Beam along +z, origin at the array center; mm / GeV / ns.
Addressing is 0-based `(ix, iy)` with `cell_id = iy*n_x + ix`. Physics list is
`FTFP_BERT` with EM option 4 (`FTFP_BERT_EMZ`). `--hadronic off` disables the
gamma/electro/muon-nuclear processes — for gamma/e- primaries that makes the
shower purely electromagnetic (used for profile tables); recorded as `hadronic_on`
in `run.json`.

The gun fires e-/gamma/pi-/mu- at 1–20 GeV **kinetic** energy, normal incidence,
with four impact modes: uniform over the face, uniform over the central crystal, a
fixed scan grid, or a fixed point (`--gun-mode point --gun-x-mm <mm> --gun-y-mm <mm>`).

Geant4 simulates the physical shower itself (every track and energy deposit). On top
of that, a per-crystal **response chain** emulates what the real detector does to
that deposited energy — the effects Geant4 does not model: light collection,
readout sparsification, photostatistics/electronics fluctuations. Per crystal the
sim stores both `edep_true` (Geant truth, no chain) and `e_vis` (after the chain),
so chain-ON vs chain-OFF can be compared offline without re-simulating. The stages,
in order, each behind a flag, all ON by default:

- **attenuation** (`--attenuation`, `--atten-length-mm`) — light collection model,
  applied per Geant4 step in `CrystalSD`: `e_att = Σ edep · exp(−dist/λ)` where
  `dist` is the distance from the step to the back (readout) face, λ = 200 cm
  (HallD `hitECAL.c`). Light born deeper in the crystal attenuates less; this is a
  deterministic geometry-dependent weighting, not a random smear. With
  `--attenuation off`, `e_att = edep`.
- **threshold** (`--hit-threshold-gev`) — readout sparsification: a cell with
  attenuated energy below 5 MeV reads exactly zero (the `p0` zeros in the LL
  profile tables come from here).
- **energy scale** (`--energy-scale`) — multiply by 1.0962 (HallD `ECALSmearer`) to
  undo the *average* light loss from attenuation, so the energy scale is right on
  average while the event-by-event attenuation fluctuations survive.
- **smearing** (`--smearing`, `--measured-res-*`, `--intrinsic-res-*`) — the random
  part, in `EventAction`. A Geant4 shower fluctuates only with the *intrinsic*
  resolution (leakage, attenuation sampling); the real detector also has
  photostatistics, PMT gain and noise fluctuations, so its *measured* resolution is
  wider. Each cell's energy is multiplied by `1 + Gauss(0, σ)` with
  `σ² = measured_res²(E) − intrinsic_res²(E)`, evaluated term by term of
  σ/E = √(stochastic²/E + noise²/E² + constant²) at that cell's own energy.
  Injecting only the quadrature *difference* (independent variances add) lands the
  total simulated resolution on the measured HallD curve without double-counting
  what Geant4 already fluctuates. Note: `--smearing off` also disables the energy
  scale — it yields attenuated raw energy (`edep_true` is the column with no chain
  at all).
- **time** (`--timing`, `--time-sigma-ns`) — optional `t_vis`: energy-weighted mean
  arrival time incl. light propagation at `--light-speed-mm-ns`, plus a 0.4 ns
  Gaussian; off by default.

Constants and citations are in [REFERENCES.md](REFERENCES.md). A single global
**calibration constant** (response stage d) removes the attenuation bias; it lives in
`configs/calib.yaml`, written by `scripts/150_calibrate.py`.

### calsim options

All flags take a value (`--flag value`); booleans accept `on/off` (also `1/true/yes`).
Lengths in mm, energies in GeV, times in ns. Resolution terms follow
σ/E = √(stochastic²/E + noise²/E² + constant²).

| Flag | Default | Meaning |
|---|---|---|
| `--crystals-nx` / `--crystals-ny` | 3 / 3 | grid size in crystals |
| `--crystal-side-mm` | 20.55 | square crystal side |
| `--crystal-length-mm` | 200.0 | crystal length |
| `--wrap-thickness-mm` | 0.175 | Tedlar wrap per lateral face (0 = no wrap volume) |
| `--particle` | e- | gun particle (e-, gamma, pi-, mu-, ...) |
| `--energy-min-gev` / `--energy-max-gev` | 5.0 / 5.0 | kinetic energy range (uniform) |
| `--energy-gev` | — | shorthand: sets min and max to the same value |
| `--gun-mode` | central | `face` \| `central` \| `grid` \| `point` |
| `--gun-grid-n` | 5 | scan points per axis for gun-mode=grid |
| `--gun-x-mm` / `--gun-y-mm` | 0 / 0 | fixed impact for gun-mode=point |
| `--events` | 100 | number of events |
| `--seed` | 1 | random seed |
| `--hadronic` | on | gamma/electro/muon-nuclear on/off; `off` = pure EM shower for gamma/e- (profile tables) |
| `--out-dir` | . | output directory (hits.csv, events.csv, run.json) |
| `--run-id` | test | run tag written into run.json |
| `--attenuation` | on | light-attenuation stage |
| `--atten-length-mm` | 2000.0 | light attenuation length in the crystal |
| `--light-speed-mm-ns` | 130.0 | effective light speed, used for hit timing |
| `--smearing` | on | smearing stage (energy scale + extra Gaussian) |
| `--energy-scale` | 1.0962 | scale applied before smearing |
| `--measured-res-stochastic` | 0.0308 | target (measured halld ECAL) resolution, stochastic term |
| `--measured-res-noise` | 0.010 | ... noise term |
| `--measured-res-constant` | 0.007 | ... constant term |
| `--intrinsic-res-stochastic` | 0.0171216 | intrinsic Geant4 resolution, stochastic term |
| `--intrinsic-res-noise` | 0.0155070 | ... noise term |
| `--hit-threshold-gev` | 0.005 | sparsification threshold, on attenuated unscaled energy |
| `--store-min-gev` | 0.0001 | store hit if edep_true above this (or e_vis > 0) |
| `--timing` | off | timing simulation (adds t_vis column) |
| `--time-sigma-ns` | 0.4 | time smearing sigma |
| `--check-overlaps` | off | Geant4 geometry overlap check |

## Data store

`calsim` writes per-job CSV (`hits.csv`, `events.csv`, `run.json`); the Python store
layer (`g4cal.store`) runs reco and converts each job into **parquet shards** under
`${paths.data_root}/store/`, one writer per shard (no concurrent-write hazard).
**duckdb** is the query layer (`scripts/views.sql` defines convenience views).

| table        | columns |
|--------------|---------|
| `runs`       | run_id, preset, nx, ny, particle, gun_mode, e_min, e_max, seed, n_events, cal, config (full params JSON), timestamp |
| `cal_hits`   | run_id, event, ix, iy, edep_true, e_vis [, t_vis] |
| `cal_events` | run_id, event, truth (pdg, ekin, x, y, dx, dy, dz), reco (n_clusters, reco_e/x/y/chi2/size) and the chain-OFF `*_true` reco on `edep_true` |

## Scripts

Run from the repo root with the venv python, Geant4 sourced.

| script | does |
|--------|------|
| `scripts/build.sh` | configure + build the C++ (wraps cmake) |
| `scripts/150_calibrate.py` | derive the global calibration constant → `configs/calib.yaml` |
| `scripts/200_production.py` | run a `(preset × particle × energy)` matrix as parallel jobs into the store |
| `scripts/400_energy_resolution.py` | σ_E/E vs energy, fit `a/√E ⊕ b/E ⊕ c`, overlay the measured model |
| `scripts/410_position_resolution.py` | position resolution vs energy and the S-shape across a cell |
| `scripts/420_event_display.py` | per-module energy displays with truth + reco overlaid |

`200_production.py --smoke` runs a tiny matrix to check the round-trip. Example
production:

```sh
python scripts/200_production.py --presets p10x10 --particles e- gamma \
    --energies 1 2 3 5 8 12 16 20 --events 10000
python scripts/400_energy_resolution.py
```

## Config

`configs/config.yaml` is one self-contained OmegaConf file: per-stage blocks
(`paths`, `detector`, `response`, `reco`, `production`, `analysis`) chained by
`${...}` interpolation. `configs/calib.yaml` is merged on top when present.
`G4CAL_CONFIG` overrides the base file; `G4CAL_PLOTS_ROOT` overrides where figures go.
