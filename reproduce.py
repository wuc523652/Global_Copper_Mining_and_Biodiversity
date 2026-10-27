"""Verify the deposit, recompute statistics, redraw figures, or run raw analyses."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def run(script: str, *arguments: str) -> None:
    env = dict(os.environ, PYTHONUTF8="1", MPLBACKEND="Agg")
    subprocess.run([sys.executable, "-B", str(ROOT / script), *arguments],
                   check=True, env=env, cwd=ROOT)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=ROOT / "data",
                        help="Complete dataset directory (default: data beside this script)")
    parser.add_argument("--output-root", type=Path,
                        help="A separate directory for generated outputs (required except check)")
    parser.add_argument("--mode", choices=("check", "calculations", "figures", "all", "raw"),
                        default="check")
    args = parser.parse_args()
    data = args.data_root.expanduser().resolve()
    if not data.is_dir():
        parser.error(f"Dataset directory does not exist: {data}. Extract the data release into data/ or supply --data-root.")
    if args.mode != "check":
        if args.output_root is None:
            parser.error("--output-root is required for calculation or rendering")
        output = args.output_root.expanduser().resolve()
        if (output.is_relative_to(data) or data.is_relative_to(output)
                or output.is_relative_to(ROOT) or ROOT.is_relative_to(output)):
            parser.error("Output must be separate from both the dataset and code directories")
        if (output / "0code").exists() or (output / ".git").exists():
            parser.error("Use a dedicated output directory, not an existing research project")
    # The deposit manifest is checked before any computation or output creation.
    run("tools/data_manifest.py", "verify", "--data-root", str(data))
    if args.mode == "check":
        print("Dataset integrity verified. This does not execute any scientific analysis.")
        return 0
    output.mkdir(parents=True, exist_ok=True)
    if args.mode in ("calculations", "all"):
        run("reproduction/recompute.py", "--data-root", str(data),
            "--output-root", str(output / "calculations"))
    if args.mode in ("figures", "all"):
        run("plotting/render_all.py", "--data-root", str(data),
            "--output-root", str(output / "figures"), "--figures", "all")
    if args.mode == "raw":
        run("run_analysis.py", "--all", "--data-root", str(data),
            "--work-root", str(output / "raw_analysis"))
    print(f"Requested mode completed: {args.mode}. Outputs: {output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as exc:
        print(f"Reproduction stopped: a component failed (exit {exc.returncode}).", file=sys.stderr)
        raise SystemExit(exc.returncode)
