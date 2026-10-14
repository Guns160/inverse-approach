#!/usr/bin/env python3
"""Run from any working directory using only this package's parent."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ηc_赝标流.__main__ import main

if __name__ == "__main__":
    main()
