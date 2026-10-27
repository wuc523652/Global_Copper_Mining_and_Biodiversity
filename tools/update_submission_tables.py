"""Apply the approved ownership presentation and Spearman corrections to a release.

Input: an existing materials folder containing Table/Supplementary and the
companion Zenodo package. Default: validate and report without writing. --apply
updates the paired Data 6-9 workbooks and the two ownership-weight CSVs in place.
No spatial analysis is run. Repeated application is byte-for-byte idempotent.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import posixpath
import re
import sys
from typing import Any
from xml.etree import ElementTree as ET
from xml.sax.saxutils import escape, unescape
from zipfile import ZipFile

import openpyxl
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
from _ownership_labels import REGISTRATION_DISPLAY_NAMES

ZENODO = Path("Zenodo/Global_Copper_Mining_and_Biodiversity")
AREA_HEADERS = ("Main cross border area (km²)", "Complete case cross border area (km²)")
OLD_RHO = {"D2": 0.9924419090176239, "G2": 0.9774228646048222}
NEW_RHO = {"D2": "0.990050709206434", "G2": "0.976950656711279"}
TEXT_REPLACEMENTS = {
    "Hong Kong": REGISTRATION_DISPLAY_NAMES["HKG"],
    "All calculation countries and territories are retained, including HKG. Area is equity-attributed mapped-range exposure, not ecological damage or legal liability.":
        "All countries and regions in the calculations are retained. Area is equity-attributed mapped-range exposure, not ecological damage or legal liability.",
    "Complete cases require 100% recorded equity and resolution of every recorded share to a registration country; 959 mines. No imputation or renormalization.":
        "Complete cases require 100% recorded equity and resolution of every recorded share to a registration country or region; 959 mines. No imputation or renormalization.",
    "Registration country resolved for the project equity holder; not necessarily ultimate control.":
        "Registration country or region resolved for the project equity holder; not necessarily ultimate control.",
    "Spearman correlation of cross-border country-pair area rankings.":
        "Spearman correlation of cross-border area rankings for pairs of countries and regions.",
    "Spearman correlation of owner-registration-country cross-border-area rankings.":
        "Spearman correlation of cross-border-area rankings for owner-registration countries and regions.",
    "Only owner-registration and mine-host countries that differ. Equity percentages summed across mines can exceed 100%; this quantity is not a single-mine ownership percentage.":
        "Only records whose owner-registration country or region differs from the mine-host country. Equity percentages summed across mines can exceed 100%; this quantity is not a single-mine ownership percentage.",
    "Sum of recorded resolved equity percentages across all mines in the country pair; can exceed 100%.":
        "Sum of recorded resolved equity percentages across all mines in the pair of countries and regions; can exceed 100%.",
    "Traceable country resolved mine equivalent": "Traceable country or region resolved mine equivalent",
    "Traceable country resolved mine equivalent.": "Traceable country or region resolved mine equivalent.",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def approved_text(value: str) -> str:
    """Transform only approved display text, never machine keys or host labels."""
    if value in TEXT_REPLACEMENTS:
        return TEXT_REPLACEMENTS[value]
    prefix = "Owner registration country"
    if value.startswith(prefix) and not value.startswith(prefix + " or region"):
        return value.replace(prefix, prefix + " or region", 1)
    return value


def snapshot(payload: bytes) -> dict[str, dict[str, object]]:
    book = openpyxl.load_workbook(io.BytesIO(payload), read_only=True, data_only=False)
    try:
        return {
            sheet.title: {cell.coordinate: cell.value for row in sheet for cell in row if cell.value is not None}
            for sheet in book
        }
    finally:
        book.close()


def area_correlation(payload: bytes, sheet: str, expected_rows: int) -> float:
    book = openpyxl.load_workbook(io.BytesIO(payload), read_only=True, data_only=True)
    try:
        rows = iter(book[sheet].values)
        headers = next(rows)
        indices = [headers.index(header) for header in AREA_HEADERS]
        pairs = [[row[index] for index in indices] for row in rows]
    finally:
        book.close()
    require(len(pairs) == expected_rows, f"Unexpected sample size in {sheet}: {len(pairs)}")
    require(all(isinstance(x, (int, float)) and math.isfinite(x) for row in pairs for x in row),
            f"Non-finite or missing area in {sheet}")
    return float(spearmanr([row[0] for row in pairs], [row[1] for row in pairs]).statistic)


def worksheet_member(archive: ZipFile, sheet_name: str) -> str:
    ns = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    workbook = ET.fromstring(archive.read("xl/workbook.xml"))
    sheets = workbook.findall("s:sheets/s:sheet", ns)
    sheet = next(item for item in sheets if item.get("name") == sheet_name)
    relation_id = sheet.get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id")
    relations = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    target = next(item.get("Target") for item in relations if item.get("Id") == relation_id)
    if target is None:
        raise ValueError(f"Missing worksheet relationship target: {sheet_name}")
    return target.lstrip("/") if target.startswith("/") else posixpath.normpath(posixpath.join("xl", target))


def rewrite_text_nodes(payload: bytes) -> bytes:
    """Patch string nodes without reserializing numbers, styles or other XML."""
    def replace(match: re.Match[bytes]) -> bytes:
        old = unescape(match[2].decode("utf-8"))
        new = approved_text(old)
        if new == old:
            return match[0]
        return match[1] + escape(new).encode("utf-8") + match[3]

    # Excel producers use shared/inline <t> strings or t="str" cells with <v>.
    prefix = rb"(?:[A-Za-z_][\w.-]*:)?"
    pattern = rb"(<" + prefix + rb"(?:t|v)(?:\s[^>]*)?>)(.*?)(</" + prefix + rb"(?:t|v)>)"
    return re.sub(pattern, replace, payload, flags=re.DOTALL)


def rewrite_table_headers(payload: bytes) -> bytes:
    """Keep Excel table-column names consistent with the display header cells."""
    def replace(match: re.Match[bytes]) -> bytes:
        old = unescape(match[2].decode("utf-8"), {"&quot;": '"', "&apos;": "'"})
        new = approved_text(old)
        if new == old:
            return match[0]
        return match[1] + escape(new, {'"': "&quot;"}).encode("utf-8") + match[3]

    return re.sub(rb'(<(?:[A-Za-z_][\w.-]*:)?tableColumn\b[^>]*\bname=")([^"]+)(")', replace, payload)


def patch_workbook(payload: bytes, number: int) -> tuple[bytes, dict]:
    before = snapshot(payload)
    if number == 6:
        for coordinate, new in NEW_RHO.items():
            old = before["Completeness_Stability"][coordinate]
            require(old in (OLD_RHO[coordinate], float(new)), f"Unexpected starting rho at {coordinate}: {old}")
    output = io.BytesIO()
    changed_members = []
    with ZipFile(io.BytesIO(payload)) as source:
        rho_member = worksheet_member(source, "Completeness_Stability") if number == 6 else None
        with ZipFile(output, "w") as destination:
            for member in source.infolist():
                original = source.read(member.filename)
                updated = original
                if member.filename == "xl/sharedStrings.xml" or member.filename.startswith("xl/worksheets/") and member.filename.endswith(".xml"):
                    updated = rewrite_text_nodes(updated)
                if member.filename.startswith("xl/tables/") and member.filename.endswith(".xml"):
                    updated = rewrite_table_headers(updated)
                if member.filename == rho_member:
                    for coordinate, new in NEW_RHO.items():
                        prefix = rb"(?:[A-Za-z_][\w.-]*:)?"
                        pattern = (rb'(<'+ prefix + rb'c\b[^>]*\br="' + coordinate.encode()
                                   + rb'"[^>]*>.*?<' + prefix + rb'v>)[^<]+(</' + prefix + rb'v>)')
                        updated, count = re.subn(pattern, lambda match: match[1] + new.encode() + match[2], updated, flags=re.DOTALL)
                        require(count == 1, f"Expected one numeric cell {coordinate}, found {count}")
                if updated != original:
                    changed_members.append(member.filename)
                destination.writestr(member, updated)
            destination.comment = source.comment
    revised = output.getvalue() if changed_members else payload
    after = snapshot(revised)
    require(before.keys() == after.keys(), "Worksheet names changed")
    changes = []
    numeric_count = 0
    for sheet, cells in before.items():
        require(cells.keys() == after[sheet].keys(), f"Cell coordinates changed in {sheet}")
        for coordinate, old in cells.items():
            updated_value = after[sheet][coordinate]
            expected = approved_text(old) if isinstance(old, str) else old
            if number == 6 and sheet == "Completeness_Stability" and coordinate in NEW_RHO:
                expected = float(NEW_RHO[coordinate])
            require(updated_value == expected, f"Unapproved change in Data {number} {sheet}!{coordinate}")
            if isinstance(old, (int, float)):
                numeric_count += 1
            if old != updated_value:
                changes.append({"sheet": sheet, "cell": coordinate, "old": old, "new": updated_value})
    return revised, {"numeric_cells_checked": numeric_count, "changes": changes,
                     "changed_xml_members": changed_members}


def patch_csv(payload: bytes) -> tuple[bytes, dict]:
    bom = payload.startswith(b"\xef\xbb\xbf")
    original_text = payload.decode("utf-8-sig")
    rows = list(csv.reader(io.StringIO(original_text, newline="")))
    require(bool(rows), "Empty ownership CSV")
    headers = rows[0]
    code_index = headers.index("investor_iso3")
    name_index = headers.index("investor_country")
    revised = [list(row) for row in rows]
    changes = []
    for line_number, row in enumerate(revised[1:], start=2):
        require(len(row) == len(headers), f"CSV field count differs on line {line_number}")
        expected = REGISTRATION_DISPLAY_NAMES.get(row[code_index])
        if expected is None:
            continue
        require(row[name_index] in ("Hong Kong", expected), f"Unexpected registration name on line {line_number}")
        if row[name_index] != expected:
            changes.append({"line": line_number, "field": "investor_country", "old": row[name_index], "new": expected})
            row[name_index] = expected
    if not changes:
        return payload, {"rows_checked": len(rows) - 1, "changes": []}
    output = io.StringIO(newline="")
    writer = csv.writer(output, lineterminator="\r\n" if "\r\n" in original_text else "\n")
    writer.writerows(revised)
    new_payload = output.getvalue().encode("utf-8-sig" if bom else "utf-8")
    parsed = list(csv.reader(io.StringIO(new_payload.decode("utf-8-sig"), newline="")))
    require(len(parsed) == len(rows) and parsed[0] == headers, "CSV rows or machine headers changed")
    for old, new in zip(rows[1:], parsed[1:]):
        expected_row = list(old)
        if old[code_index] in REGISTRATION_DISPLAY_NAMES:
            expected_row[name_index] = REGISTRATION_DISPLAY_NAMES[old[code_index]]
        require(new == expected_row, "Unapproved CSV field change")
    return new_payload, {"rows_checked": len(rows) - 1, "changes": changes}


def sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def prepare_updates(root: Path) -> tuple[list[tuple[Path, bytes, bytes]], dict[str, Any]]:
    table_dir = root / "Table/Supplementary"
    zenodo = root / ZENODO
    inputs = {number: (table_dir / f"Supplementary_Data_{number}.xlsx").read_bytes() for number in range(6, 10)}
    owner_rho = area_correlation(inputs[6], "Owner_Rank_Comparison", 78)
    pair_rho = area_correlation(inputs[7], "Pair_Rank_Comparison", 183)
    for coordinate, rho in (("D2", owner_rho), ("G2", pair_rho)):
        require(math.isclose(rho, float(NEW_RHO[coordinate]), rel_tol=0, abs_tol=1e-14),
                f"Recomputed coefficient differs from the approved result at {coordinate}: {rho}")
    updates = []
    audit: dict[str, Any] = {"materials_root": str(root), "standard_spearman": {"owner_n": 78, "owner_rho": owner_rho,
             "bilateral_n": 183, "bilateral_rho": pair_rho}, "files": []}
    for number, before in inputs.items():
        name = f"Supplementary_Data_{number}.xlsx"
        paths = [table_dir / name, zenodo / "8Final_Results" / name]
        require(paths[1].read_bytes() == before, f"Table and Zenodo copies differ before update: {name}")
        after, details = patch_workbook(before, number)
        for path in paths:
            updates.append((path, before, after))
        audit["files"].append({"paths": [str(path.relative_to(root)) for path in paths],
                               "before_sha256": sha256(before), "after_sha256": sha256(after), **details})
    for name in ("Mine_OwnerCountry_DerivedWeights_30km_Primary_HKGIncluded.csv",
                 "Mine_OwnerCountry_DerivedWeights_30km_CompleteCase.csv"):
        path = zenodo / "7Ownership_Attribution" / name
        before = path.read_bytes()
        after, details = patch_csv(before)
        updates.append((path, before, after))
        audit["files"].append({"paths": [str(path.relative_to(root))], "before_sha256": sha256(before),
                               "after_sha256": sha256(after), **details})
    return updates, audit


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--materials-root", required=True, type=Path)
    parser.add_argument("--apply", action="store_true", help="Write validated changes (default: dry run)")
    parser.add_argument("--audit-json", type=Path, help="Optional detailed verification report")
    args = parser.parse_args()
    root = args.materials_root.expanduser().resolve()
    updates, audit = prepare_updates(root)
    changed = [(path, before, after) for path, before, after in updates if before != after]
    audit["mode"] = "apply" if args.apply else "dry_run"
    audit["changed_files"] = len(changed)
    if args.audit_json:
        report_path = args.audit_json.expanduser().resolve()
        require(report_path not in {path.resolve() for path, _, _ in updates}, "Audit output overlaps a data file")
    if args.apply:
        # Complete all checks before the first write and detect concurrent edits.
        for path, before, _ in updates:
            require(path.read_bytes() == before, f"File changed during validation: {path}")
        for path, _, after in changed:
            temporary = path.with_name(path.name + ".submission-update.tmp")
            with temporary.open("xb") as stream:
                stream.write(after)
            temporary.replace(path)
        for path, _, after in updates:
            require(path.read_bytes() == after, f"Post-write verification failed: {path}")
    if args.audit_json:
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"mode": audit["mode"], "changed_files": len(changed),
                      "standard_spearman": audit["standard_spearman"]}, indent=2))


if __name__ == "__main__":
    main()
