#!/usr/bin/env python3
from __future__ import annotations

import argparse
import html
import re
import subprocess
from pathlib import Path


LEVEL_ROW_RE = re.compile(
    r"<tr><td class=['\"]cell elvl['\"]>(.*?)</td>.*?<td class=['\"]cellc t12['\"]>(.*?)</td>",
    re.S,
)
NUMBER_RE = re.compile(r"[-+]?\d+(?:\.\d+)?")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare PDG2025 nuclei entries with NNDC NuDat levels.")
    parser.add_argument(
        "--list-file",
        default="input/list/PDG2025_codex/list-withexcitednuclei.dat",
        help="Path to Thermal-FIST particle list to validate.",
    )
    parser.add_argument(
        "--cache-dir",
        default="input/list/PDG2025_codex/pdglisting/nndc-cache",
        help="Directory for cached NNDC HTML pages.",
    )
    parser.add_argument(
        "--fetch",
        action="store_true",
        help="Fetch/update NNDC pages with curl before validating.",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=8,
        help="Fetch retries per isotope when --fetch is used.",
    )
    parser.add_argument(
        "--curl-timeout",
        type=int,
        default=25,
        help="curl max-time in seconds for each attempt.",
    )
    parser.add_argument(
        "--energy-tol-kev",
        type=float,
        default=50.0,
        help="Energy mismatch threshold in keV.",
    )
    parser.add_argument(
        "--width-tol-mev",
        type=float,
        default=0.05,
        help="Width mismatch threshold in MeV.",
    )
    return parser.parse_args()


def parse_list_rows(path: Path) -> dict[str, list[dict[str, float | int | str]]]:
    symbol_by_z = {
        1: "H",
        2: "He",
        3: "Li",
    }
    rows_by_iso: dict[str, list[dict[str, float | int | str]]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            pdgid = int(parts[0])
            if abs(pdgid) < 1000000000:
                continue
            baryon = int(parts[6])
            charge = int(parts[7])
            strange = int(parts[8])
            charm = int(parts[9])
            if strange != 0 or charm != 0:
                continue
            if charge not in symbol_by_z:
                continue
            isotope = f"{baryon}{symbol_by_z[charge]}"
            row = {
                "pdgid": pdgid,
                "name": parts[1],
                "mass": float(parts[3]),
                "width": float(parts[12]),
            }
            rows_by_iso.setdefault(isotope, []).append(row)
    for isotope in rows_by_iso:
        rows_by_iso[isotope].sort(key=lambda item: float(item["mass"]))
    return rows_by_iso


def fetch_page(isotope: str, out_path: Path, retries: int, timeout: int) -> bool:
    url = f"https://www.nndc.bnl.gov/nudat3/getdataset.jsp?nucleus={isotope}&unc=NDS"
    for _ in range(max(retries, 1)):
        cmd = [
            "curl",
            "-sS",
            "-L",
            "--max-time",
            str(timeout),
            url,
            "-o",
            str(out_path),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0 and out_path.exists() and out_path.stat().st_size > 0:
            return True
    return False


def clean_html_text(fragment: str) -> str:
    text = re.sub(r"<br\s*/?>", " ", fragment, flags=re.I)
    text = re.sub(r"<[^>]+>", " ", text)
    return " ".join(html.unescape(text).split())


def parse_nndc_levels(html_text: str) -> list[tuple[float, float | None]]:
    if "Empty Data Set" in html_text:
        return []
    levels: list[tuple[float, float | None]] = []
    for energy_cell, width_cell in LEVEL_ROW_RE.findall(html_text):
        energy_text = clean_html_text(energy_cell)
        if "+X" in energy_text:
            continue
        energy_match = NUMBER_RE.search(energy_text)
        if energy_match is None:
            continue
        energy_kev = float(energy_match.group())

        width_text = clean_html_text(width_cell).upper()
        if "STABLE" in width_text:
            width_mev: float | None = 0.0
        else:
            width_match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*MEV", width_text)
            width_mev = float(width_match.group(1)) if width_match else None
        levels.append((energy_kev, width_mev))

    levels.sort(key=lambda value: value[0])
    deduped: list[tuple[float, float | None]] = []
    seen_energies: set[float] = set()
    for entry in levels:
        if entry[0] in seen_energies:
            continue
        seen_energies.add(entry[0])
        deduped.append(entry)
    return deduped


def compare_rows(
    rows_by_iso: dict[str, list[dict[str, float | int | str]]],
    levels_by_iso: dict[str, list[tuple[float, float | None]]],
    energy_tol_kev: float,
    width_tol_mev: float,
) -> tuple[int, list[str], list[str]]:
    checked = 0
    missing_isotopes: list[str] = []
    mismatch_lines: list[str] = []

    for isotope in sorted(rows_by_iso):
        if isotope not in levels_by_iso:
            missing_isotopes.append(isotope)
            continue
        nndc_levels = levels_by_iso[isotope]
        if not nndc_levels:
            continue
        rows = rows_by_iso[isotope]
        ground_mass = float(rows[0]["mass"])
        for row in rows:
            excitation_kev = (float(row["mass"]) - ground_mass) * 1.0e6
            width_mev = float(row["width"]) * 1.0e3
            best = min(nndc_levels, key=lambda item: abs(item[0] - excitation_kev))
            delta_energy = excitation_kev - best[0]
            delta_width = None if best[1] is None else width_mev - best[1]
            checked += 1

            bad_energy = abs(delta_energy) > energy_tol_kev
            bad_width = delta_width is not None and abs(delta_width) > width_tol_mev
            if bad_energy or bad_width:
                width_ref = "NA" if best[1] is None else f"{best[1]:.3f}"
                width_delta = "NA" if delta_width is None else f"{delta_width:+.3f}"
                mismatch_lines.append(
                    f"{isotope:>4s}  {int(row['pdgid']):>12d}  {str(row['name']):<12s}  "
                    f"dE={delta_energy:+8.1f} keV  dW={width_delta:>8s} MeV  "
                    f"(NNDC E={best[0]:.1f} keV, W={width_ref} MeV)"
                )

    return checked, missing_isotopes, mismatch_lines


def main() -> None:
    args = parse_args()
    list_path = Path(args.list_file).resolve()
    cache_dir = Path(args.cache_dir).resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)

    rows_by_iso = parse_list_rows(list_path)
    if not rows_by_iso:
        raise RuntimeError(f"No non-strange nuclei rows found in {list_path}")

    fetched_ok: list[str] = []
    fetch_failed: list[str] = []
    if args.fetch:
        for isotope in sorted(rows_by_iso):
            out_path = cache_dir / f"{isotope}.html"
            if fetch_page(isotope, out_path, retries=args.retries, timeout=args.curl_timeout):
                fetched_ok.append(isotope)
            else:
                fetch_failed.append(isotope)

    levels_by_iso: dict[str, list[tuple[float, float | None]]] = {}
    for isotope in sorted(rows_by_iso):
        html_path = cache_dir / f"{isotope}.html"
        if not html_path.exists():
            continue
        text = html_path.read_text(encoding="utf-8", errors="ignore")
        levels_by_iso[isotope] = parse_nndc_levels(text)

    checked, missing_isotopes, mismatches = compare_rows(
        rows_by_iso=rows_by_iso,
        levels_by_iso=levels_by_iso,
        energy_tol_kev=args.energy_tol_kev,
        width_tol_mev=args.width_tol_mev,
    )

    print("NNDC nuclei check")
    print(f"- Input list: {list_path}")
    print(f"- Cache dir:  {cache_dir}")
    if args.fetch:
        print(f"- Fetch success: {len(fetched_ok)} ({', '.join(fetched_ok) if fetched_ok else 'none'})")
        print(f"- Fetch failed:  {len(fetch_failed)} ({', '.join(fetch_failed) if fetch_failed else 'none'})")
    print(f"- Isotopes in list: {', '.join(sorted(rows_by_iso))}")
    print(f"- Isotopes with cached NNDC pages: {', '.join(sorted(levels_by_iso)) if levels_by_iso else 'none'}")
    print(f"- States checked: {checked}")
    print(f"- Mismatches: {len(mismatches)}")

    if missing_isotopes:
        print("Missing isotopes (no cached NNDC page):")
        for isotope in missing_isotopes:
            print(f"  - {isotope}")

    if mismatches:
        print("Mismatched states:")
        for line in mismatches:
            print(f"  - {line}")


if __name__ == "__main__":
    main()
