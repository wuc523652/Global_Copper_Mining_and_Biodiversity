"""List, check, or run the essential manuscript analyses in a separate workspace."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import runpy
import sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parent
MANIFEST = json.loads((ROOT / "analysis_manifest.json").read_text(encoding="utf-8"))
sys.path.insert(0, str(ROOT / "analysis"))


def configure(args):
    from _paths import validate_roots
    data_root, work_root = validate_roots(args.data_root, args.work_root)
    os.environ["NC_DATA_ROOT"] = str(data_root)
    os.environ["NC_WORK_ROOT"] = str(work_root)
    os.environ["NC_CAP_CHECK"] = "1" if args.capped_check else "0"
    os.environ["MPLBACKEND"] = "Agg"
    os.environ["MPLCONFIGDIR"] = str(work_root / ".matplotlib")


def raw_requirements():
    from _paths import data_path
    return {
        "gadm": [Path(data_path("Global_GADM/GC_GADM.shp"))],
        "footprint": [Path(data_path("Global_Copper_Mine/Global_Copper_Mine/Global_Copper_Mine-allinfo.shp"))],
        "buffer": [Path(data_path("Global_Copper_Mine/Global_Copper_Mine/Buffer-30km/Global_Copper_Mine-30km.shp"))],
        "iucn": [Path(data_path(f"IUCNredlist/202502/{taxon}-{status}/data_0.shp"))
                 for taxon in ("Amphibians", "Birds", "Fish", "Mammals", "Reptiles") for status in ("CR", "EN", "VU")],
        "ownership": [Path(data_path("Global_Copper_Mine/Excel/" + name)) for name in
                      ("Global_Copper_Mine-regionrenew.xlsx", "Copper_mine_all_info.xlsx", "Copper_mine_ownerlist-countryrenew.xlsx")],
    }


def load_module(file, name):
    spec = importlib.util.spec_from_file_location(name, file)
    if spec is None or spec.loader is None:
        raise ImportError(file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def prepare_buffer():
    """Reuse the original strict local-AEQD function; write no raw input files."""
    import geopandas as gpd
    from _paths import analysis_script, work_root
    destination = work_root() / "prepared_inputs/primary_buffer.gpkg"
    if destination.exists():
        raise FileExistsError(f"Buffer already exists: {destination}")
    source = raw_requirements()["footprint"][0]
    if not source.is_file():
        raise FileNotFoundError(source)
    common = load_module(analysis_script("15_buffer_distance_sensitivity.py"), "nc_buffer_functions")
    footprints = gpd.read_file(source, engine="pyogrio")
    if footprints.crs is None:
        raise ValueError("Mine footprints must have a declared CRS")
    footprints = footprints.to_crs("EPSG:4326")
    buffers = []
    for number, geometry in enumerate(footprints.geometry, 1):
        buffered, _ = common.strict_buffer_wgs84(geometry, 30.0)
        if buffered is None or buffered.is_empty:
            raise ValueError(f"Invalid buffer for mine row {number}")
        buffers.append(buffered)
    footprints.geometry = buffers
    destination.parent.mkdir(parents=True, exist_ok=True)
    footprints.to_file(destination, layer="primary_buffer", driver="GPKG", engine="pyogrio")
    print(f"Prepared {len(footprints)} 30-km buffers: {destination}")


def missing_intermediates(selected, work_root):
    """Allow prerequisites produced by earlier selected steps, but check resumed inputs."""
    planned_directories = []
    missing = []
    for record in selected:
        for relative in record["intermediate_inputs"]:
            path = work_root / relative
            will_be_produced = any(path.is_relative_to(directory) and step in record["dependencies"]
                                   for step, directory in planned_directories)
            if not path.is_file() and not will_be_produced:
                missing.append(path)
        planned_directories.append((record["step"], work_root / record["output_directory"]))
    return missing


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true", help="List the scientific steps without executing them")
    parser.add_argument("--check", action="store_true", help="Check external inputs and prerequisites for selected steps; run no analysis")
    parser.add_argument("--steps", nargs="+", choices=[r["step"] for r in MANIFEST])
    parser.add_argument("--all", action="store_true", help="Run all steps including spatial sensitivities")
    parser.add_argument("--data-root", type=Path, default=ROOT / "inputs")
    parser.add_argument("--work-root", type=Path, default=ROOT / "work")
    parser.add_argument("--prepare-buffer", action="store_true", help="Build the primary buffer from licensed mine footprints")
    parser.add_argument("--capped-check", action="store_true", help="Run the previously reported 1,000,000-km2 cap check in step 08")
    args = parser.parse_args()
    if args.all and args.steps:
        parser.error("Choose --all or --steps, not both")
    if args.check and args.prepare_buffer:
        parser.error("--check does not execute preparation; run --prepare-buffer separately")
    if args.list:
        for record in MANIFEST:
            print(record["step"], record["purpose"], "[dependencies:", ",".join(record["dependencies"]) or "none", "]")
        return 0
    if not (args.check or args.steps or args.all or args.prepare_buffer):
        parser.print_help()
        return 0
    configure(args)
    if args.prepare_buffer:
        prepare_buffer()
    chosen = set(args.steps or [r["step"] for r in MANIFEST])
    selected = [r for r in MANIFEST if r["step"] in chosen]
    if args.check or args.all or args.steps:
        groups = {g for r in MANIFEST if r["step"] in chosen for g in r["raw_inputs"]}
        missing = [p for group, paths in raw_requirements().items() if group in groups for p in paths if not p.is_file()]
        # A shapefile needs its attributes, index and CRS sidecars as well.
        missing.extend(p.with_suffix(ext) for group, paths in raw_requirements().items() if group in groups
                       for p in paths if p.suffix.lower() == ".shp" and p.is_file()
                       for ext in (".shx", ".dbf", ".prj") if not p.with_suffix(ext).is_file())
        missing.extend(missing_intermediates(selected, args.work_root.expanduser().resolve()))
        if missing:
            print("Missing required input files (including intermediate outputs for resumed steps):")
            for path in missing:
                print(path)
            return 2
        print("Input preflight passed. Earlier selected steps must produce their listed outputs before dependent steps run.")
    if args.check:
        return 0
    if args.all or args.steps:
        for record in MANIFEST:
            if record["step"] in chosen:
                missing = [args.work_root.expanduser().resolve() / p for p in record["intermediate_inputs"]
                           if not (args.work_root.expanduser().resolve() / p).is_file()]
                if missing:
                    raise FileNotFoundError("Missing intermediate inputs: " + ", ".join(map(str, missing)))
                print("\nRUN", record["step"], record["purpose"], flush=True)
                runpy.run_path(str(ROOT / record["file"]), run_name="__main__")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, FileNotFoundError, FileExistsError) as exc:
        print(f"Cannot run: {exc}", file=sys.stderr)
        raise SystemExit(2)
