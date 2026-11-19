"""Shared pytest fixtures/path setup.

src/python is not an installed package in this repo layout, so make the
pipeline modules importable from there.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_PYTHON = PROJECT_ROOT / "src" / "python"

if str(SRC_PYTHON) not in sys.path:
    sys.path.insert(0, str(SRC_PYTHON))
