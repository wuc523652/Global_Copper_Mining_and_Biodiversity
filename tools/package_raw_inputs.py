"""Copy study-matched local inputs without changing their bytes or access terms.

This prepares a LOCAL reproduction package. It does not upload files, grant a
redistribution licence, or fetch replacement data from current provider releases.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil


SHAPEFILES = [
    "Global_GADM/GC_GADM.shp",
    "Global_Copper_Mine/Global_Copper_Mine/Global_Copper_Mine-allinfo.shp",
    "Global_Copper_Mine/Global_Copper_Mine/Buffer-30km/Global_Copper_Mine-30km.shp",
    *[f"IUCNredlist/202502/{taxon}-{status}/data_0.shp"
      for taxon in ("Amphibians", "Birds", "Fish", "Mammals", "Reptiles")
      for status in ("CR", "EN", "VU")],
]
WORKBOOKS = [
    f"Global_Copper_Mine/Excel/{name}" for name in (
        "Global_Copper_Mine-regionrenew.xlsx", "Copper_mine_all_info.xlsx",
        "Copper_mine_ownerlist-countryrenew.xlsx",
    )
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def input_files(source: Path) -> list[str]:
    files = list(WORKBOOKS)
    for relative in SHAPEFILES:
        stem = Path(relative)
        files.extend(stem.with_suffix(ext).as_posix() for ext in (".shp", ".shx", ".dbf", ".prj"))
        if (source / stem.with_suffix(".cpg")).is_file():
            files.append(stem.with_suffix(".cpg").as_posix())
    return sorted(files)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True, help="Companion dataset root")
    args = parser.parse_args()
    source = args.source_root.expanduser().resolve()
    dataset = args.data_root.expanduser().resolve()
    target = dataset / "10Raw_Inputs"
    files = input_files(source)
    missing = [relative for relative in files if not (source / relative).is_file()]
    if missing:
        parser.error("Missing study inputs: " + ", ".join(missing))
    if source == target or source.is_relative_to(target):
        parser.error("The source cannot be the destination or a directory inside it")
    records = []
    for relative in files:
        src, dst = source / relative, target / relative
        expected = sha256(src)
        if dst.exists() and sha256(dst) != expected:
            raise FileExistsError(f"Refusing to overwrite a different existing input: {dst}")
        if not dst.exists():
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
        actual = sha256(dst)
        if actual != expected:
            raise RuntimeError(f"Copy checksum mismatch: {relative}")
        records.append({"path": relative, "bytes": dst.stat().st_size, "sha256": actual})
    manifest = {
        "purpose": "Local, study-matched inputs for the complete spatial workflow",
        "access": "Original provider conditions apply; not certified for public redistribution",
        "provider_versions": {
            "iucn": "February 2025 local snapshot",
            "gadm": "Historical release not established; exact supplied bytes are fingerprinted",
            "ownership": "Historical snapshot date not established; exact supplied bytes are fingerprinted",
        },
        "files": records,
    }
    (target / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"files": len(records), "bytes": sum(r["bytes"] for r in records), "destination": str(target)}))


if __name__ == "__main__":
    main()
