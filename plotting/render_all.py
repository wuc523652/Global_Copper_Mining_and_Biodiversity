"""Reproduce manuscript figures using the code package and complete dataset."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.dont_write_bytecode = True
from data_paths import external_output_root, resolve_data_root

CODE = Path(__file__).resolve().parent
ALL_FIGURES = [str(number) for number in range(1, 8)] + [f"S{number}" for number in range(1, 10)]


def selected_figures(values: list[str]) -> list[str]:
    selected: list[str] = []
    groups = {"all": ALL_FIGURES, "main": ALL_FIGURES[:7], "supplementary": ALL_FIGURES[7:]}
    for value in values:
        additions = groups.get(value.lower(), [value.upper()])
        for item in additions:
            if item not in ALL_FIGURES:
                raise ValueError(f"Unknown figure {value!r}; use 1–7, S1–S9, all, main, or supplementary")
            if item not in selected:
                selected.append(item)
    return selected


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, help="Dataset root; overrides COPPER_DATA_ROOT and bundled data/")
    parser.add_argument("--output-root", type=Path, required=True, help="Directory outside both downloaded packages")
    parser.add_argument("--figures", nargs="+", default=["all"], help="all, main, supplementary, or selected IDs, e.g. 1 7 S2")
    parser.add_argument("--formats", nargs="+", choices=["png", "pdf", "svg"], default=["png"])
    args = parser.parse_args()
    try:
        data_root = resolve_data_root(args.data_root)
        if data_root is None:
            raise ValueError("Extract the dataset into data/, supply --data-root, or set COPPER_DATA_ROOT")
        output_root = external_output_root(args.output_root, data_root)
        figures = selected_figures(args.figures)
    except ValueError as error:
        parser.error(str(error))
    output_root.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment.update({"COPPER_DATA_ROOT": str(data_root), "PYTHONDONTWRITEBYTECODE": "1", "MPLBACKEND": "Agg"})
    report = {"data_root": str(data_root), "output_root": str(output_root), "figures": [], "complete": False}
    report_path = output_root / "figure_render_report.json"
    for identifier in figures:
        is_supplementary = identifier.startswith("S")
        number = int(identifier[1:] if is_supplementary else identifier)
        output = output_root / ("Supplementary" if is_supplementary else "Main")
        if is_supplementary:
            commands = [[sys.executable, "-B", str(CODE / "supplementary" / "render_supplementary.py"),
                         "--data-root", str(data_root), "--output", str(output), "--figures", str(number),
                         "--formats", *dict.fromkeys(args.formats)]]
        else:
            commands = [[sys.executable, "-B", str(CODE / "main" / f"figure_{number:02d}.py"),
                         "--data-root", str(data_root), "--output-dir", str(output), "--format", extension]
                        for extension in dict.fromkeys(args.formats)]
        print(f"Rendering {identifier} ...", flush=True)
        started = time.perf_counter()
        returncode = 0
        for command in commands:
            result = subprocess.run(command, env=environment, check=False)
            if result.returncode:
                returncode = result.returncode
                break
        report["figures"].append({"figure": identifier, "elapsed_seconds": round(time.perf_counter() - started, 3),
                                  "returncode": returncode, "formats": list(dict.fromkeys(args.formats))})
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        if returncode:
            print(f"Figure {identifier} failed; see its error above. Report: {report_path}", file=sys.stderr)
            return returncode
    report["complete"] = True
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Rendered {len(figures)} figures. Report: {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
