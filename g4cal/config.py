"""Config loader (OmegaConf idiom from submodules/simulation-pipeline).

Usage:
    from g4cal import load_config
    cfg = load_config()                       # configs/config.yaml (+ calib.yaml merge)
    cfg = load_config(overrides=["production.workers=8"])

G4CAL_CONFIG env var selects an alternative base YAML.
"""
import os
from pathlib import Path

from omegaconf import OmegaConf

REPO = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = REPO / "configs" / "config.yaml"


def config_path():
    """The base YAML that will be loaded (G4CAL_CONFIG env var or the default)."""
    return Path(os.environ.get("G4CAL_CONFIG", str(DEFAULT_CONFIG)))


def calib_path():
    """calib.yaml sits next to the base config, so a foreign config (e.g.
    cal-fpga's) merges its own calibration, not g4cal's."""
    return config_path().parent / "calib.yaml"


def load_config(overrides=None):
    cfg = OmegaConf.load(config_path())
    calib = calib_path()
    if calib.exists():
        cfg = OmegaConf.merge(cfg, OmegaConf.load(calib))
    if overrides:
        cfg = OmegaConf.merge(cfg, OmegaConf.from_dotlist(list(overrides)))
    return cfg
