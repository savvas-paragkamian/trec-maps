"""Configuration: YAML parameters, transect definitions, small shared helpers."""
from __future__ import annotations

import gzip
import json
import re
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


def _transects_doc(path: Path | str | None = None) -> dict:
    with open(Path(path) if path else TRANSECTS_FILE) as f:
        return yaml.safe_load(f)


def all_transects(params: dict, path: Path | str | None = None) -> dict[str, dict]:
    """The explicit entries of transects.yml plus one auto entry per remaining LSI site of the exports."""
    doc = _transects_doc(path)
    transects = {k: dict(v) for k, v in (doc.get("transects") or {}).items()}
    auto = doc.get("auto") or {}
    if auto.get("enabled", True):
        from . import sites
        taken = {t.get("land_site") for t in transects.values()}
        for name, cfg in sites.auto_transects(params, auto, taken).items():
            deep_update(cfg, (auto.get("overrides") or {}).get(name) or {})
            transects.setdefault(name, cfg)
    for name, cfg in transects.items():
        cfg["name"] = name
    return transects


def output_stem(tcfg: dict) -> str:
    """The name of a transect's figures and of their folder: the land site's ID first, then the transect
    name without its `_lsi_<n>` (the ID carries that number): thrace_lsi_1 -> 042-01_thrace."""
    return f"{tcfg['land_site']}_{re.sub(r'_lsi_[0-9]+$', '', tcfg['name'])}"


def load_transect(name: str, params: dict | None = None, path: Path | str | None = None) -> dict:
    """A transect by its name (thrace_lsi_1) or by the name of its figures (042-01_thrace)."""
    transects = all_transects(params or load_params(), path)
    by_stem = {output_stem(t): t for t in transects.values()}
    if name not in transects and name not in by_stem:
        raise SystemExit(f"unknown transect {name!r}; `make list` prints the known ones")
    return transects.get(name) or by_stem[name]


def transect_dirs(params: dict, tcfg: dict) -> tuple[Path, Path, str]:
    """(derived_dir/<name>, output_dir/<stem>, stem), created; a figure is `out / f"{stem}_profile.png"`."""
    stem = output_stem(tcfg)
    derived = params["paths"]["derived_dir"] / tcfg["name"]
    out = params["paths"]["output_dir"] / stem
    derived.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    return derived, out, stem


def read_json(path: Path | str):
    """A JSON file, gzipped when its name ends in .gz."""
    path = Path(path)
    with (gzip.open if path.suffix == ".gz" else open)(path, "rt", encoding="utf-8") as f:
        return json.load(f)


def write_json(obj, dest: Path | str) -> None:
    """Write JSON (gzipped when the name ends in .gz) as `<dest>.part`, then rename."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    with (gzip.open if dest.suffix == ".gz" else open)(part, "wt", encoding="utf-8") as f:
        json.dump(obj, f, separators=(",", ":"))
    part.rename(dest)


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


def deep_update(base: dict, overrides: dict) -> dict:
    """Merge `overrides` into `base` key by key (nested dicts merge, anything else replaces)."""
    for key, value in (overrides or {}).items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            deep_update(base[key], value)
        else:
            base[key] = value
    return base


def apply_transect_overrides(params: dict, tcfg: dict) -> dict:
    """A transect may override the grid, profile, schematic and render sections of params.yml."""
    for section in ("grid", "profile", "schematic", "render"):
        if isinstance(tcfg.get(section), dict):
            deep_update(params[section], tcfg[section])
    return params
