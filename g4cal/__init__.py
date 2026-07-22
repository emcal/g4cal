"""g4cal: PbWO4 crystal-calorimeter Geant4 simulation + data store + reconstruction.

The C++ half (calsim, calreco) builds via CMake; this package is the Python half:
config loading, running sim+reco jobs into a parquet store, and the analysis/plots.

    from g4cal import load_config, run_job, store_glob
"""
from g4cal.config import load_config
from g4cal.store import Job, make_run_id, run_job, store_glob, PARTICLE_TAGS

__all__ = ["load_config", "Job", "make_run_id", "run_job", "store_glob", "PARTICLE_TAGS"]
