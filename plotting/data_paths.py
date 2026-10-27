"""Explicit paths shared by the standalone figure entry points."""
from __future__ import annotations

import argparse
import os
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[1]


def resolve_data_root(value: str | Path | None = None) -> Path | None:
    """Use the CLI root, then COPPER_DATA_ROOT, then the bundled data directory."""
    value = value or os.environ.get("COPPER_DATA_ROOT")
    if not value:
        value = PACKAGE_ROOT / "data"
        if not value.is_dir():
            return None
    root = Path(value).expanduser().resolve()
    for folder in ("8Final_Results", "9Figure_Inputs"):
        if not (root / folder).is_dir():
            raise ValueError(f"Data root must contain {folder}: {root}")
    return root


def data_root_from_arguments() -> Path | None:
    """Resolve main-figure constants before their drawing modules import them."""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--data-root", type=Path)
    args, _ = parser.parse_known_args()
    return resolve_data_root(args.data_root)


def external_output_root(value: Path, data_root: Path) -> Path:
    """Keep generated files outside the downloaded code and data packages."""
    output = value.expanduser().resolve()
    code_root = Path(__file__).resolve().parents[1]
    for protected in (code_root, data_root):
        if output == protected or protected in output.parents:
            raise ValueError(f"Choose an output directory outside the code/data packages: {output}")
    return output
