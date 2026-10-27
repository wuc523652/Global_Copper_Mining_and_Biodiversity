"""Build or verify a SHA-256 inventory of the complete local dataset."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

MANIFEST = "deposit_manifest.json"
REQUIRED = ("1Copper_Mine_Data", "2Species_Inventory", "3Grid_and_Biodiversity",
            "4Mine_Buffers", "5Spatial_Overlap", "6Randomization",
            "7Ownership_Attribution", "8Final_Results", "9Figure_Inputs", "10Raw_Inputs",
            "11Reanalysis_Checkpoints")


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def inventory(root: Path) -> list[dict]:
    records = []
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"Deposit must contain regular files, not links: {path}")
        if path.is_file() and path != root / MANIFEST:
            records.append({"path": path.relative_to(root).as_posix(),
                            "bytes": path.stat().st_size, "sha256": digest(path)})
    return records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("build", "verify"))
    parser.add_argument("--data-root", type=Path, default=Path(__file__).resolve().parents[1] / "data",
                        help="Complete dataset directory (default: data in the code package)")
    args = parser.parse_args()
    root = args.data_root.expanduser().resolve()
    if not root.is_dir():
        parser.error(f"Not a dataset directory: {root}")
    missing = [name for name in REQUIRED if not (root / name).is_dir()]
    if missing:
        parser.error("Missing dataset folders: " + ", ".join(missing))
    target = root / MANIFEST
    if args.action == "build":
        records = inventory(root)
        target.write_text(json.dumps({"schema_version": 1, "algorithm": "sha256",
                                     "files": records}, indent=2) + "\n", encoding="utf-8")
        print(f"Inventoried {len(records)} files ({sum(r['bytes'] for r in records):,} bytes).")
        return 0
    if not target.is_file():
        parser.error(f"Missing deposit manifest: {target}")
    payload = json.loads(target.read_text(encoding="utf-8"))
    if payload.get("schema_version") != 1 or payload.get("algorithm") != "sha256":
        raise ValueError("Unsupported deposit manifest format")
    seen = set()
    failures = []
    for record in payload["files"]:
        relative = record["path"]
        path = (root / relative).resolve()
        if not path.is_relative_to(root) or path == root / MANIFEST or relative in seen:
            raise ValueError(f"Unsafe or duplicate manifest path: {relative}")
        seen.add(relative)
        if not path.is_file():
            failures.append(f"Missing: {relative}")
        elif path.stat().st_size != record["bytes"] or digest(path) != record["sha256"]:
            failures.append(f"Changed: {relative}")
    actual = {p.relative_to(root).as_posix() for p in root.rglob("*")
              if p.is_file() and p != target}
    failures.extend(f"Unlisted: {p}" for p in sorted(actual - seen))
    if failures:
        print("Dataset verification failed:\n" + "\n".join(failures), file=sys.stderr)
        return 2
    print(f"Verified {len(seen)} files against {MANIFEST}.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, KeyError, TypeError, OSError) as exc:
        print(f"Cannot verify dataset: {exc}", file=sys.stderr)
        raise SystemExit(2)
