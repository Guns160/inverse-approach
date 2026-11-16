#!/usr/bin/env python3
"""Run from any working directory using only this package's parent."""
from importlib import import_module
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
main = import_module("eta_c_pseudoscalar+BW.__main__").main

if __name__ == "__main__":
    main()
