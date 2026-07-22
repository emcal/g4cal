# References: halld_sim / hdds constant harvest

One-time harvest (2026-07-16) for the standalone PbWO4 simulation. Repos cloned to
`/data/cal-fpga/external/halld_sim` (HEAD 4db3352, 2026-07-15) and
`/data/cal-fpga/external/hdds` (HEAD 2026-04-21); paths below relative to those roots.
GlueX "ECAL" = their PbWO4 Compton calorimeter — same 2x2x20 cm class of block as ours.
No halld code is linked or run; everything below is hardcoded in our sources with citations.

## 1. hitECAL (HDGeant truth-hit formation)

File: `src/programs/Simulation/HDGeant/hitECAL.c`. Defaults overridable from CCDB
`ECAL/ecal_parms` (lines 62-99); repo defaults used.

| constant | value | units | how applied (formula) | citation |
|---|---|---|---|---|
| ATTEN_LENGTH | 200 | cm | `dEcorr = dEsum * exp(-dist/ATTEN_LENGTH)`, `dist = 0.5*LENGTH_OF_BLOCK - z_local` = distance from step midpoint to the **downstream/back end** (readout end) of the 20 cm crystal. NOT from the front face. | hitECAL.c:27,172-173 |
| C_EFFECTIVE | 13 | cm/ns | `tcorr = t + dist/C_EFFECTIVE`, `t` = step-midpoint global time | hitECAL.c:28,174 |
| LENGTH_OF_BLOCK | 20 | cm | enters `dist` only | hitECAL.c:30,172 |
| WIDTH_OF_BLOCK | 2 | cm | declared, **never used** | hitECAL.c:29 |
| TWO_HIT_RESOL | 75 | ns | same-block hits within 75 ns merged (energy-weighted time), else new hit | hitECAL.c:31,194-221 |
| MAX_HITS | 100 | — | max merged hits per block | hitECAL.c:32,217-227 |
| THRESH_MEV | 5 | MeV | (a) truth-shower record cut; (b) block hits with attenuated E <= 5 MeV dropped at pickup | hitECAL.c:33,135,275 |
| addressing | row/column, 0-based in HDDM | — | key `((row+1)<<16)+(column+1)` | hitECAL.c:169-183 |

## 2. ECALSmearer (mcsmear)

Files: `src/programs/Simulation/mcsmear/ECALSmearer.{cc,h}`. Defaults below (CCDB
`ECAL/mc_energy`, `ECAL/mc_time` override in production; not in repo).

| constant | value | units | how applied (formula) | citation |
|---|---|---|---|---|
| ECAL_EN_SCALE | 1.0962 | — | `E *= 1.0962` applied to truth-hit energy **first, unconditionally, before smearing** ("new calibration of the ECAL"); smearing sigma uses the scaled E | ECALSmearer.cc:13,80-83 |
| ECAL_EN_P0 (measured stochastic) | 3.08e-2 | GeV^1/2 | measured fractional variance `de_e_expect = (P0/sqrt(E))^2 + (P1/E)^2 + P2^2` | ECALSmearer.cc:16,88-89 |
| ECAL_EN_P1 (measured noise) | 1.0e-2 | GeV | as above | ECALSmearer.cc:17 |
| ECAL_EN_P2 (measured constant) | 0.7e-2 | — | as above | ECALSmearer.cc:18 |
| ECAL_EN_GP0 (Geant stochastic) | 1.71216e-2 | GeV^1/2 | intrinsic Geant variance `de_e_geant = (GP0/sqrt(E))^2 + (GP1/E)^2` | ECALSmearer.cc:21,92 |
| ECAL_EN_GP1 (Geant noise) | 1.55070e-2 | GeV | as above | ECALSmearer.cc:22 |
| ECAL_EN_GP2 | 0.0 | — | declared, never used | ECALSmearer.cc:23 |
| smearing sigma | — | fractional | `sig_res = sqrt(de_e_expect - de_e_geant)`; if > 0: `E *= (1 + Gaussian(sig_res))` — per-hit, multiplicative | ECALSmearer.cc:94-97 |
| ECAL_TSIGMA | 0.4 | ns | `t += Gaussian(0.4)` | ECALSmearer.cc:26,99 |
| ECAL_BLOCK_THRESHOLD | 15 | MeV | **dead code** — threshold check commented out; only the 5 MeV hitECAL cut is active | ECALSmearer.cc:30,105-111 |

Pileup/merging (later phases): ecal_max_hits = 3, ecal_min_delta_t_ns = 70,
ecal_integration_window_ns = 64 (hddm_s_merger.cc:78-80).

## 3. Freshness check (time-boxed)

halld_sim HEAD (2026-07-15) carries no fresher ECAL resolution constants; the
ECALSmearer defaults are byte-identical to the older CCAL (PbWO4 Compton calorimeter)
defaults (CCALSmearer.cc:15-32) — beam-test-era PbWO4 numbers, not a dedicated new-ECAL
calibration. Production values live in CCDB (not in the repo). Conclusion: use the
ECALSmearer defaults; anything fresher needs a CCDB query.

## 4. Geometry / material (hdds)

File: `CrystalECAL_HDDS.xml` (2023-12-17), materials in `Material_HDDS.xml`.

| constant | value | units | notes | citation |
|---|---|---|---|---|
| crystal (XTBL) | 2.055 x 2.055 x 20.0 | cm | GlueX ECAL crystal (ours is 2.0 x 2.0 x 20 per spec) | CrystalECAL_HDDS.xml:140-141 |
| module envelope (XTMD) | 2.09 x 2.09 x 20.0 | cm | Tedlar box -> **175 um Tedlar per lateral face, no extra air gap**; pitch = envelope | CrystalECAL_HDDS.xml:136-137 |
| lattice pitch | 2.09 | cm | dX = dY in all placements | CrystalECAL_HDDS.xml:50,83 |
| PWO density | 8.28 | g/cm3 | mass fractions W 0.4040, Pb 0.4554, O 0.1406 (stoichiometric PbWO4); G4 NIST `G4_PbWO4` matches (8.28) | Material_HDDS.xml:1359-1370 |
| Tedlar | 1.49 g/cm3, C2H3F (atoms 2/3/1) | | wrapping material | Material_HDDS.xml:697-702 |
| NIST/PDG cross-check | PbWO4 8.30 g/cm3, X0 0.8903 cm, R_M 1.959 cm | | NIST-sourced, for reference | PDG Atomic/Nuclear properties |

## Caveats

- hitECAL truth-hit energies are attenuation-corrected already; mcsmear's 1.0962
  compensates that plus sampling loss. We mirror the same order: attenuate ->
  5 MeV threshold -> scale -> smear.
- The 15 MeV block threshold is dead code in mcsmear; we keep it as an OFF-by-default
  reco knob only.
- Smearing silently skips when measured variance <= intrinsic variance (never happens
  with these defaults).
- `WIDTH_OF_BLOCK` (2 cm) disagrees with the 2.055 cm geometry but is unused.
