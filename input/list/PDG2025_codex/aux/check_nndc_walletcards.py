#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
from dataclasses import dataclass
from pathlib import Path


ISOTOPE_BY_CHARGE = {
    1: "H",
    2: "He",
    3: "Li",
}


@dataclass
class LocalLevel:
    pdgid: int
    name: str
    stable_flag: int
    mass_gev: float
    width_mev: float
    excitation_kev: float


@dataclass
class WalletLevel:
    level_index: int
    excitation_kev: float
    stable_flag: bool
    width_mev: float | None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare light nuclei in PDG2025_codex with NNDC Wallet Cards."
    )
    parser.add_argument(
        "--list-file",
        default="input/list/PDG2025_codex/list-withexcitednuclei.dat",
        help="Thermal-FIST list file to validate.",
    )
    parser.add_argument(
        "--walletcards-url",
        default="https://www.nndc.bnl.gov/walletcards/StandardSearchServlet",
        help="NNDC Wallet Cards search endpoint.",
    )
    parser.add_argument(
        "--cache-dir",
        default="input/list/PDG2025_codex/pdglisting/nndc-walletcards-cache",
        help="Directory for cached JSON responses.",
    )
    parser.add_argument(
        "--energy-tol-kev",
        type=float,
        default=150.0,
        help="Energy matching tolerance in keV.",
    )
    parser.add_argument(
        "--width-tol-mev",
        type=float,
        default=0.12,
        help="Width mismatch tolerance in MeV.",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=5,
        help="Number of network retries per isotope.",
    )
    parser.add_argument(
        "--curl-timeout",
        type=int,
        default=25,
        help="curl timeout in seconds.",
    )
    parser.add_argument(
        "--summary-out",
        default="input/list/PDG2025_codex/pdglisting/nndc-walletcards/check-summary.tsv",
        help="Summary TSV output path.",
    )
    parser.add_argument(
        "--details-out",
        default="input/list/PDG2025_codex/pdglisting/nndc-walletcards/check-details.tsv",
        help="Detailed TSV output path.",
    )
    parser.add_argument(
        "--cache-only",
        action="store_true",
        help="Skip network fetch and use cached JSON responses only.",
    )
    return parser.parse_args()


def parse_local_rows(path: Path) -> dict[str, list[LocalLevel]]:
    rows_by_iso: dict[str, list[dict[str, float | int | str]]] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
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
        if charge not in ISOTOPE_BY_CHARGE:
            continue

        isotope = f"{baryon}{ISOTOPE_BY_CHARGE[charge]}"
        rows_by_iso.setdefault(isotope, []).append(
            {
                "pdgid": pdgid,
                "name": parts[1],
                "stable_flag": int(parts[2]),
                "mass_gev": float(parts[3]),
                "width_mev": float(parts[12]) * 1.0e3,
            }
        )

    result: dict[str, list[LocalLevel]] = {}
    for isotope, rows in rows_by_iso.items():
        ordered = sorted(rows, key=lambda item: float(item["mass_gev"]))
        ground_mass = float(ordered[0]["mass_gev"])
        result[isotope] = [
            LocalLevel(
                pdgid=int(row["pdgid"]),
                name=str(row["name"]),
                stable_flag=int(row["stable_flag"]),
                mass_gev=float(row["mass_gev"]),
                width_mev=float(row["width_mev"]),
                excitation_kev=(float(row["mass_gev"]) - ground_mass) * 1.0e6,
            )
            for row in ordered
        ]
    return result


def get_width_mev(entry: dict) -> float | None:
    if entry.get("stable") is True:
        return 0.0
    decay_width = entry.get("decayWidth")
    if not isinstance(decay_width, dict):
        return None
    value = decay_width.get("value")
    unit = decay_width.get("unit")
    if value is None or unit is None:
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if unit == "MeV":
        return numeric
    if unit == "keV":
        return numeric / 1.0e3
    if unit == "eV":
        return numeric / 1.0e6
    return None


def fetch_walletcards(
    isotope: str,
    endpoint: str,
    cache_file: Path,
    retries: int,
    timeout: int,
) -> list[dict]:
    payload = json.dumps(
        {
            "nuclide": isotope,
            "element": None,
            "zMin": None,
            "zMax": None,
            "nMin": None,
            "nMax": None,
            "aMin": None,
            "aMax": None,
        },
        separators=(",", ":"),
    )

    last_error = None
    for _ in range(max(retries, 1)):
        cmd = [
            "curl",
            "-sS",
            "-L",
            "--max-time",
            str(timeout),
            endpoint,
            "-H",
            "Content-Type: application/json",
            "--data",
            payload,
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            last_error = proc.stderr.strip() or f"curl return code {proc.returncode}"
            continue
        try:
            parsed = json.loads(proc.stdout)
            if not isinstance(parsed, list):
                last_error = "response is not a JSON list"
                continue
            cache_file.write_text(proc.stdout, encoding="utf-8")
            return parsed
        except json.JSONDecodeError as exc:
            last_error = f"invalid JSON response: {exc}"
            continue

    raise RuntimeError(f"Failed to fetch Wallet Cards for {isotope}: {last_error}")


def parse_wallet_levels(entries: list[dict]) -> list[WalletLevel]:
    levels: list[WalletLevel] = []
    for entry in entries:
        level = entry.get("levelEnergy", {})
        value = level.get("value")
        if value is None:
            continue
        try:
            excitation_kev = float(value)
        except (TypeError, ValueError):
            continue
        levels.append(
            WalletLevel(
                level_index=int(entry.get("levelIndex", 0)),
                excitation_kev=excitation_kev,
                stable_flag=bool(entry.get("stable", False)),
                width_mev=get_width_mev(entry),
            )
        )
    return sorted(levels, key=lambda item: (item.excitation_kev, item.level_index))


def match_levels(
    local_levels: list[LocalLevel],
    wallet_levels: list[WalletLevel],
    energy_tol_kev: float,
) -> tuple[list[tuple[LocalLevel, WalletLevel, float]], list[LocalLevel], list[WalletLevel]]:
    matched: list[tuple[LocalLevel, WalletLevel, float]] = []
    used_local_indices: set[int] = set()
    wallet_only: list[WalletLevel] = []

    for wallet in wallet_levels:
        best_idx = None
        best_delta = None
        for index, local in enumerate(local_levels):
            if index in used_local_indices:
                continue
            delta = abs(local.excitation_kev - wallet.excitation_kev)
            if best_delta is None or delta < best_delta:
                best_delta = delta
                best_idx = index
        if best_idx is not None and best_delta is not None and best_delta <= energy_tol_kev:
            used_local_indices.add(best_idx)
            matched.append((local_levels[best_idx], wallet, best_delta))
        else:
            wallet_only.append(wallet)

    local_unmatched = [row for index, row in enumerate(local_levels) if index not in used_local_indices]
    return matched, local_unmatched, wallet_only


def main() -> None:
    args = parse_args()
    list_file = Path(args.list_file).resolve()
    cache_dir = Path(args.cache_dir).resolve()
    summary_out = Path(args.summary_out).resolve()
    details_out = Path(args.details_out).resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)

    local_by_iso = parse_local_rows(list_file)
    isotopes = sorted(local_by_iso)
    if not isotopes:
        raise RuntimeError(f"No non-strange nuclei found in {list_file}")

    wallet_by_iso: dict[str, list[WalletLevel]] = {}
    failed: list[str] = []
    source_by_iso: dict[str, str] = {}
    fetched_live = 0
    loaded_cache = 0
    for isotope in isotopes:
        cache_file = cache_dir / f"{isotope}.json"
        payload = None
        if not args.cache_only:
            try:
                payload = fetch_walletcards(
                    isotope=isotope,
                    endpoint=args.walletcards_url,
                    cache_file=cache_file,
                    retries=args.retries,
                    timeout=args.curl_timeout,
                )
                fetched_live += 1
                source_by_iso[isotope] = "live"
            except RuntimeError:
                payload = None

        if payload is None and cache_file.exists():
            try:
                payload = json.loads(cache_file.read_text(encoding="utf-8"))
                if not isinstance(payload, list):
                    payload = None
                else:
                    loaded_cache += 1
                    source_by_iso[isotope] = "cache"
            except json.JSONDecodeError:
                payload = None

        if payload is None:
            failed.append(isotope)
            source_by_iso[isotope] = "missing"
            wallet_by_iso[isotope] = []
        else:
            wallet_by_iso[isotope] = parse_wallet_levels(payload)

    if fetched_live == 0 and loaded_cache == 0:
        missing = ", ".join(failed) if failed else "none"
        raise RuntimeError(
            "No Wallet Cards data available (network and cache both unavailable). "
            f"Missing isotopes: {missing}"
        )

    summary_lines = [
        "isotope\tsource\tlocal_levels\twallet_levels\tmatched\twallet_only\tlocal_unmatched\twidth_mismatch"
    ]
    detail_lines = [
        "isotope\tcategory\tname_or_level\tpdgid\texcitation_kev_local\texcitation_kev_wallet\twidth_mev_local\twidth_mev_wallet\tdelta_energy_kev\tdelta_width_mev"
    ]

    for isotope in isotopes:
        local_levels = local_by_iso.get(isotope, [])
        wallet_levels = wallet_by_iso.get(isotope, [])
        matched, local_unmatched, wallet_only = match_levels(
            local_levels=local_levels,
            wallet_levels=wallet_levels,
            energy_tol_kev=args.energy_tol_kev,
        )

        width_mismatch = 0
        for local, wallet, delta_e in matched:
            if wallet.width_mev is None:
                continue
            delta_w = local.width_mev - wallet.width_mev
            if abs(delta_w) > args.width_tol_mev:
                width_mismatch += 1
                detail_lines.append(
                    f"{isotope}\twidth-mismatch\t{local.name}\t{local.pdgid}\t"
                    f"{local.excitation_kev:.1f}\t{wallet.excitation_kev:.1f}\t"
                    f"{local.width_mev:.3f}\t{wallet.width_mev:.3f}\t"
                    f"{delta_e:+.1f}\t{delta_w:+.3f}"
                )

        for wallet in wallet_only:
            wallet_name = f"idx{wallet.level_index}@{wallet.excitation_kev:.1f}keV"
            wallet_width = "NA" if wallet.width_mev is None else f"{wallet.width_mev:.3f}"
            detail_lines.append(
                f"{isotope}\twallet-only\t{wallet_name}\tNA\tNA\t{wallet.excitation_kev:.1f}\tNA\t{wallet_width}\tNA\tNA"
            )

        for local in local_unmatched:
            detail_lines.append(
                f"{isotope}\tlocal-unmatched\t{local.name}\t{local.pdgid}\t{local.excitation_kev:.1f}\tNA\t"
                f"{local.width_mev:.3f}\tNA\tNA\tNA"
            )

        summary_lines.append(
            f"{isotope}\t{source_by_iso[isotope]}\t{len(local_levels)}\t{len(wallet_levels)}\t{len(matched)}\t"
            f"{len(wallet_only)}\t{len(local_unmatched)}\t{width_mismatch}"
        )

    summary_out.parent.mkdir(parents=True, exist_ok=True)
    summary_out.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    details_out.parent.mkdir(parents=True, exist_ok=True)
    details_out.write_text("\n".join(detail_lines) + "\n", encoding="utf-8")

    print("NNDC Wallet Cards nuclei check")
    print(f"- Input list: {list_file}")
    print(f"- Endpoint:   {args.walletcards_url}")
    print(f"- Cache dir:  {cache_dir}")
    print(f"- Isotopes:   {', '.join(isotopes)}")
    print(f"- Live fetch: {fetched_live}")
    print(f"- Cache load: {loaded_cache}")
    print(f"- Failed:     {', '.join(failed) if failed else 'none'}")
    for line in summary_lines:
        print(line)
    print(f"- Summary written to: {summary_out}")
    print(f"- Details written to: {details_out}")


if __name__ == "__main__":
    main()
