"""Shared argument parsing for the numbered scripts (run them from anywhere)."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from trecmaps import config as C  # noqa: E402


def parser(doc: str) -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=doc)
    ap.add_argument("--transect", default="anavyssos", help="name in config/transects.yml")
    ap.add_argument("--params", default=None, help="alternative params.yml")
    return ap


def setup(args) -> tuple[dict, dict]:
    return C.load_params(args.params), C.load_transect(args.transect)
