"""Job execution + parquet store (calsim -> calreco -> parquet).

One job = (preset, particle, energy point/range, chunk) -> calsim -> calreco (twice:
on e_vis with the global calibration, and on edep_true with cal=1 for chain-OFF
comparisons) -> parquet shards under <store>/{runs,cal_hits,cal_events}/<run_id>.parquet.
One writer per shard: no concurrent-write hazard by construction.

Store schema (see README.md):
  runs        run_id, preset, nx, ny, particle, gun_mode, e_min, e_max, seed, n_events,
              cal, config (full calsim params JSON), timestamp
  cal_hits    run_id, event, ix, iy, edep_true, e_vis [, t_vis]
  cal_events  run_id, event, pdg, ekin, x, y, dx, dy, dz,
              n_clusters, reco_e, reco_x, reco_y, reco_chi2, reco_size,
              n_clusters_true, reco_e_true, reco_x_true, reco_y_true
"""
import json
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd


@dataclass
class Job:
    run_id: str
    preset: str
    nx: int
    ny: int
    particle: str
    e_min: float
    e_max: float
    gun_mode: str
    n_events: int
    seed: int
    cal: float = 1.0
    extra_calsim_args: list = field(default_factory=list)  # e.g. ["--attenuation", "off"]


PARTICLE_TAGS = {"e-": "em", "gamma": "g", "pi-": "pim", "mu-": "mum"}


def make_run_id(preset, particle, e_min, e_max, chunk):
    tag = f"e{e_min:g}" if e_min == e_max else f"eU{e_min:g}-{e_max:g}"
    return f"{preset}-{PARTICLE_TAGS[particle]}-{tag}-c{chunk:02d}"


def calreco_binary(cfg, nx):
    return f"{cfg.paths.build}/calreco_{nx}x{nx}"


def run_job(job: Job, cfg) -> str:
    """Run one sim+reco+convert job. Returns run_id. Raises on any step failure."""
    tmp = Path(cfg.paths.tmp) / job.run_id
    tmp.mkdir(parents=True, exist_ok=True)
    Path(cfg.paths.logs).mkdir(parents=True, exist_ok=True)
    log = Path(cfg.paths.logs) / f"{job.run_id}.log"
    pitch = float(cfg.detector.pitch_mm)
    r = cfg.response

    calsim_cmd = [
        f"{cfg.paths.build}/calsim",
        "--crystals-nx", str(job.nx), "--crystals-ny", str(job.ny),
        "--crystal-side-mm", str(cfg.detector.crystal_side_mm),
        "--crystal-length-mm", str(cfg.detector.crystal_length_mm),
        "--wrap-thickness-mm", str(cfg.detector.wrap_thickness_mm),
        "--particle", job.particle,
        "--energy-min-gev", str(job.e_min), "--energy-max-gev", str(job.e_max),
        "--gun-mode", job.gun_mode,
        "--events", str(job.n_events), "--seed", str(job.seed),
        "--out-dir", str(tmp), "--run-id", job.run_id,
        "--attenuation", "on" if r.attenuation.enabled else "off",
        "--atten-length-mm", str(r.attenuation.length_mm),
        "--light-speed-mm-ns", str(r.attenuation.light_speed_mm_ns),
        "--smearing", "on" if r.smearing.enabled else "off",
        "--energy-scale", str(r.smearing.energy_scale),
        "--measured-res-stochastic", str(r.smearing.measured_res.stochastic),
        "--measured-res-noise", str(r.smearing.measured_res.noise),
        "--measured-res-constant", str(r.smearing.measured_res.constant),
        "--intrinsic-res-stochastic", str(r.smearing.intrinsic_res.stochastic),
        "--intrinsic-res-noise", str(r.smearing.intrinsic_res.noise),
        "--hit-threshold-gev", str(r.smearing.hit_threshold_gev),
        "--store-min-gev", str(r.store_min_gev),
        "--timing", "on" if r.time.enabled else "off",
        "--time-sigma-ns", str(r.time.sigma_ns),
    ] + list(job.extra_calsim_args)

    reco_base = [
        calreco_binary(cfg, job.nx),
        "--hits", str(tmp / "hits.csv"), "--profile", str(cfg.paths.ilreco_profile),
        "--nx", str(job.nx), "--ny", str(job.ny), "--pitch", str(pitch),
        "--min-e", str(cfg.reco.block_threshold_gev),
    ]
    with open(log, "w") as lf:
        subprocess.run(calsim_cmd, stdout=lf, stderr=subprocess.STDOUT, check=True)
        subprocess.run(reco_base + ["--out", str(tmp / "reco.csv"),
                                    "--energy-col", "e_vis", "--cal", str(job.cal)],
                       stdout=lf, stderr=subprocess.STDOUT, check=True)
        subprocess.run(reco_base + ["--out", str(tmp / "reco_true.csv"),
                                    "--energy-col", "edep_true", "--cal", "1.0"],
                       stdout=lf, stderr=subprocess.STDOUT, check=True)

    convert_to_store(job, cfg, tmp)
    shutil.rmtree(tmp)
    return job.run_id


def convert_to_store(job: Job, cfg, tmp: Path):
    store = Path(cfg.paths.store)
    for sub in ("runs", "cal_hits", "cal_events"):
        (store / sub).mkdir(parents=True, exist_ok=True)

    hits = pd.read_csv(tmp / "hits.csv")
    events = pd.read_csv(tmp / "events.csv")
    reco = pd.read_csv(tmp / "reco.csv")
    reco_true = pd.read_csv(tmp / "reco_true.csv").rename(columns={
        "n_clusters": "n_clusters_true", "e": "reco_e_true",
        "x": "reco_x_true", "y": "reco_y_true"})[
        ["event", "n_clusters_true", "reco_e_true", "reco_x_true", "reco_y_true"]]
    reco = reco.rename(columns={"e": "reco_e", "x": "reco_x", "y": "reco_y",
                                "chi2": "reco_chi2", "size": "reco_size"})

    ev = events.merge(reco, on="event", how="left").merge(reco_true, on="event", how="left")
    ev["n_clusters"] = ev["n_clusters"].fillna(0).astype("int32")
    ev["n_clusters_true"] = ev["n_clusters_true"].fillna(0).astype("int32")
    ev.insert(0, "run_id", job.run_id)
    hits.insert(0, "run_id", job.run_id)

    run_json = json.loads((tmp / "run.json").read_text())
    manifest = pd.DataFrame([{
        "run_id": job.run_id, "preset": job.preset, "nx": job.nx, "ny": job.ny,
        "particle": job.particle, "gun_mode": job.gun_mode,
        "e_min": job.e_min, "e_max": job.e_max, "seed": job.seed,
        "n_events": run_json["n_processed"], "cal": job.cal,
        "config": json.dumps(run_json), "timestamp": run_json["timestamp"],
    }])

    hits.to_parquet(store / "cal_hits" / f"{job.run_id}.parquet", index=False)
    ev.to_parquet(store / "cal_events" / f"{job.run_id}.parquet", index=False)
    manifest.to_parquet(store / "runs" / f"{job.run_id}.parquet", index=False)


def store_glob(cfg, table):
    return f"{cfg.paths.store}/{table}/*.parquet"
