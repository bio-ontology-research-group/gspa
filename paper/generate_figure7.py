#!/usr/bin/env python3
"""Thin wrapper: the figure scripts live in paper-reproducibility/figures/ (see README.md there)."""
import runpy
import sys
from pathlib import Path

script = Path(__file__).resolve().parents[1] / "paper-reproducibility" / "figures" / "generate_figures.py"
sys.argv = [str(script)] + ['7'] + sys.argv[1:]
runpy.run_path(str(script), run_name="__main__")
