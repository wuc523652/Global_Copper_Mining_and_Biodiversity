"""Resolve separate input and output roots for the analysis scripts."""
import os
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[1]


def _contained(root: Path, relative: str) -> Path:
    result = (root / relative).resolve()
    if not result.is_relative_to(root):
        raise ValueError("A configured relative path escapes its root")
    return result


def work_root() -> Path:
    _, work = validate_roots(
        Path(os.environ.get("NC_DATA_ROOT", str(PACKAGE_ROOT / "data" / "10Raw_Inputs"))),
        Path(os.environ.get("NC_WORK_ROOT", str(PACKAGE_ROOT / "work"))),
    )
    return work


def validate_roots(data: Path, work: Path) -> tuple[Path, Path]:
    data, work = data.expanduser().resolve(), work.expanduser().resolve()
    if work.is_relative_to(data) or data.is_relative_to(work):
        raise ValueError("Use separate, non-overlapping input and output directories")
    if PACKAGE_ROOT.is_relative_to(work) or any(
        work.is_relative_to(PACKAGE_ROOT / name)
        for name in ("analysis", "demo", "data", "plotting", "reproduction", "tools", "validation")
    ):
        raise ValueError("Do not use the code package or its source/data directories as the output directory")
    if (work / "0code").is_dir() or (work / ".git").exists():
        raise ValueError("Use a separate output directory, not an existing research project")
    return data, work


def project_path(relative: str = "") -> str:
    return str(_contained(work_root(), relative))


def data_path(relative: str = "") -> str:
    root = Path(os.environ.get("NC_DATA_ROOT", str(PACKAGE_ROOT / "data" / "10Raw_Inputs"))).expanduser().resolve()
    if relative == "Global_Copper_Mine/Global_Copper_Mine/Buffer-30km/Global_Copper_Mine-30km.shp":
        generated = work_root() / "prepared_inputs/primary_buffer.gpkg"
        if generated.is_file():
            return str(generated)
    return str(_contained(root, relative))


def analysis_script(name: str) -> str:
    if Path(name).name != name:
        raise ValueError("Expected a single analysis filename")
    path = PACKAGE_ROOT / "analysis" / name
    if not path.is_file():
        raise FileNotFoundError(path)
    return str(path)
