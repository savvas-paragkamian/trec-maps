"""Configuration: YAML parameters, transect definitions, small shared helpers."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
PARAMS_FILE = ROOT / "config" / "params.yml"
TRANSECTS_FILE = ROOT / "config" / "transects.yml"


def log(step: str, msg: str) -> None:
    print(f"[{step}] {msg}", flush=True)


def load_params(path: Path | str | None = None) -> dict:
    """params.yml with every `paths.*` entry made absolute against the repo root."""
    with open(Path(path) if path else PARAMS_FILE) as f:
        params = yaml.safe_load(f)
    for key, value in params["paths"].items():
        p = Path(value)
        params["paths"][key] = p if p.is_absolute() else ROOT / p
    return params


def load_transect(name: str, path: Path | str | None = None) -> dict:
    with open(Path(path) if path else TRANSECTS_FILE) as f:
        transects = yaml.safe_load(f)["transects"]
    if name not in transects:
        raise SystemExit(f"unknown transect {name!r}; known: {', '.join(transects)}")
    cfg = dict(transects[name])
    cfg["name"] = name
    return cfg


def transect_dirs(params: dict, name: str) -> tuple[Path, Path]:
    """(derived_dir/name, output_dir/name), created."""
    derived = params["paths"]["derived_dir"] / name
    out = params["paths"]["output_dir"] / name
    derived.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    return derived, out


def hex_to_rgb(color: str) -> np.ndarray:
    h = color.lstrip("#")
    return np.array([int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4)], dtype=np.float32)


def apply_overrides(cfg: dict, items: list[str] | None) -> dict:
    """`--set camera.azimuth_deg=200 exaggeration=4` style overrides, parsed as YAML scalars."""
    for item in items or []:
        key, sep, value = item.partition("=")
        if not sep:
            raise SystemExit(f"--set expects key=value, got {item!r}")
        node = cfg
        parts = key.split(".")
        for part in parts[:-1]:
            node = node[part]
        node[parts[-1]] = yaml.safe_load(value)
    return cfg
