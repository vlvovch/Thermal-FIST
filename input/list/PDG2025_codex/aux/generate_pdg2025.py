#!/usr/bin/env python3
"""
Generate Thermal-FIST PDG2025 hadron list files.

Workflow:
1) Start from PDG2020 charm-inclusive hadron membership and conventions.
2) Update masses/widths and known BRs from PDG API 2025 SQLite.
3) Expand with additional light/strange/charmed states present in PDG2025 (MCID-based).
4) Fill missing decay information via nearest-state templates with matching
   conserved quantum numbers and renormalize BR to 100%.
5) Write a single canonical decays.dat containing all generated decays.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
import re
import sqlite3


@dataclass
class ParticleEntry:
    pdgid: int
    name: str
    stable: int
    mass: float
    degeneracy: float
    statistics: int
    baryon: int
    charge: int
    strange: int
    charm: int
    abs_s: float
    abs_c: float
    width: float
    threshold: float


@dataclass
class DecayChannel:
    br: float
    daughters: List[int]


EXCLUDED_EXPANSION_MCIDS = {130, 310}
EXCLUDED_ESTABLISHED_WITH_CHARM_PDGS = {4422}
LIGHT_DECAY_REFERENCE_OVERRIDE_IDS = {
    9010225,  # f(2)(1565)
    9050221,  # f(0)(2020)
    9070225,  # f(2)(2150)
    9010315,  # K(2)*(1980)0
    9010325,  # K(2)*(1980)+
}
MANUAL_DECAY_OVERRIDE_CHANNELS: Dict[int, List[DecayChannel]] = {
    # Enforce Viewer-style representation for K(1)(1650):
    # K pi pi + K phi.
    9000313: [
        DecayChannel(br=0.9, daughters=[311, 211, -211]),
        DecayChannel(br=0.1, daughters=[311, 333]),
    ],
    9000323: [
        DecayChannel(br=0.9, daughters=[321, 211, -211]),
        DecayChannel(br=0.1, daughters=[321, 333]),
    ],
    # D(s1)(2536): keep dominant measured D* K modes explicit (incl. D*(2007)0 K+),
    # preserve small measured D pi K mode, keep upper-limit channels conservative,
    # and assign the residual to D_s* gamma.
    10433: [
        DecayChannel(br=0.31, daughters=[413, 311]),
        DecayChannel(br=0.36, daughters=[423, 321]),
        DecayChannel(br=0.01, daughters=[411, -211, 321]),
        DecayChannel(br=0.08, daughters=[411, 311]),
        DecayChannel(br=0.04, daughters=[421, 321]),
        DecayChannel(br=0.20, daughters=[433, 22]),
    ],
}
DISPLAY_NAME_OVERRIDES = {
    9000313: "K(1)(1650)0",
    9000323: "K(1)(1650)+",
    9010225: "f(2)(1565)",
    9050221: "f(0)(2020)",
    9070225: "f(2)(2150)",
}
EXTERNAL_MASSES_GEV = {
    11: 0.000510999,
    12: 0.0,
    13: 0.10565837,
    14: 0.0,
    15: 1.77686,
    16: 0.0,
    22: 0.0,
}
MANUAL_ESTABLISHED_STATE_OVERRIDES = (
    # D and Ds established states without DB MCIDs.
    {
        "mcid": 417,
        "parent_pdgid": "M203",
        "name": "D3*(2750)+",
        "flags": "M",
        "quantum_j": "3",
        "charge": 1.0,
    },
    {
        "mcid": 427,
        "parent_pdgid": "M203",
        "name": "D3*(2750)0",
        "flags": "M",
        "quantum_j": "3",
        "charge": 0.0,
    },
    {
        "mcid": 30433,
        "parent_pdgid": "M182",
        "name": "D(s1)*(2700)+",
        "flags": "M",
        "quantum_j": "1",
        "charge": 1.0,
    },
    {
        "mcid": 437,
        "parent_pdgid": "M226",
        "name": "D(s3)*(2860)+",
        "flags": "M",
        "quantum_j": "3",
        "charge": 1.0,
    },
    # Charmonium established states without DB MCIDs.
    {
        "mcid": 9060443,
        "parent_pdgid": "M074",
        "name": "psi(4230)",
        "flags": "M",
        "quantum_j": "1",
        "charge": 0.0,
    },
    {
        "mcid": 9080443,
        "parent_pdgid": "M181",
        "name": "psi(4360)",
        "flags": "M",
        "quantum_j": "1",
        "charge": 0.0,
    },
    {
        "mcid": 9090443,
        "parent_pdgid": "M189",
        "name": "psi(4660)",
        "flags": "M",
        "quantum_j": "1",
        "charge": 0.0,
    },
    {
        "mcid": 9040441,
        "parent_pdgid": "M159",
        "name": "chi(c0)(3915)",
        "flags": "M",
        "quantum_j": "0",
        "charge": 0.0,
    },
    {
        "mcid": 9030443,
        "parent_pdgid": "M176",
        "name": "chi(c1)(3872)",
        "flags": "M",
        "quantum_j": "1",
        "charge": 0.0,
    },
    {
        "mcid": 9050443,
        "parent_pdgid": "M193",
        "name": "chi(c1)(4140)",
        "flags": "M",
        "quantum_j": "1",
        "charge": 0.0,
    },
    {
        "mcid": 9070443,
        "parent_pdgid": "M233",
        "name": "chi(c1)(4274)",
        "flags": "M",
        "quantum_j": "1",
        "charge": 0.0,
    },
    {
        "mcid": 30445,
        "parent_pdgid": "M212",
        "name": "psi(2)(3823)",
        "flags": "M",
        "quantum_j": "2",
        "charge": 0.0,
    },
    {
        "mcid": 447,
        "parent_pdgid": "M241",
        "name": "psi(3)(3842)",
        "flags": "M",
        "quantum_j": "3",
        "charge": 0.0,
    },
    # Charmed-baryon established states without DB MCIDs.
    {
        "mcid": 204122,
        "parent_pdgid": "B178",
        "name": "Lambda(c)(2860)+",
        "flags": "B",
        "quantum_j": "3/2",
        "charge": 1.0,
    },
    {
        "mcid": 304122,
        "parent_pdgid": "B122",
        "name": "Lambda(c)(2940)+",
        "flags": "B",
        "quantum_j": "3/2",
        "charge": 1.0,
    },
    {
        "mcid": 204314,
        "parent_pdgid": "B130",
        "name": "Xi(c)(2970)0",
        "flags": "B",
        "quantum_j": "1/2",
        "charge": 0.0,
    },
    {
        "mcid": 204324,
        "parent_pdgid": "B130",
        "name": "Xi(c)(2970)+",
        "flags": "B",
        "quantum_j": "1/2",
        "charge": 1.0,
    },
    {
        "mcid": 9104324,
        "parent_pdgid": "B157",
        "name": "Xi(c)(3055)+",
        "flags": "B",
        "quantum_j": "3/2",
        "charge": 1.0,
    },
    {
        "mcid": 9204314,
        "parent_pdgid": "B147",
        "name": "Xi(c)(3080)0",
        "flags": "B",
        "quantum_j": "5/2",
        "charge": 0.0,
    },
    {
        "mcid": 9204324,
        "parent_pdgid": "B147",
        "name": "Xi(c)(3080)+",
        "flags": "B",
        "quantum_j": "5/2",
        "charge": 1.0,
    },
    {
        "mcid": 9004112,
        "parent_pdgid": "B155",
        "name": "Sigma(c)(2800)0",
        "flags": "B",
        "quantum_j": "3/2",
        "charge": 0.0,
    },
    {
        "mcid": 9004212,
        "parent_pdgid": "B155",
        "name": "Sigma(c)(2800)+",
        "flags": "B",
        "quantum_j": "3/2",
        "charge": 1.0,
    },
    {
        "mcid": 9004222,
        "parent_pdgid": "B155",
        "name": "Sigma(c)(2800)++",
        "flags": "B",
        "quantum_j": "3/2",
        "charge": 2.0,
    },
    {
        "mcid": 9104332,
        "parent_pdgid": "B173",
        "name": "Omega(c)(3000)",
        "flags": "B",
        "quantum_j": "1/2",
        "charge": 0.0,
    },
    {
        "mcid": 9204332,
        "parent_pdgid": "B174",
        "name": "Omega(c)(3050)",
        "flags": "B",
        "quantum_j": "1/2",
        "charge": 0.0,
    },
    {
        "mcid": 9104334,
        "parent_pdgid": "B175",
        "name": "Omega(c)(3065)",
        "flags": "B",
        "quantum_j": "3/2",
        "charge": 0.0,
    },
    {
        "mcid": 9204334,
        "parent_pdgid": "B176",
        "name": "Omega(c)(3090)",
        "flags": "B",
        "quantum_j": "3/2",
        "charge": 0.0,
    },
    {
        "mcid": 9104336,
        "parent_pdgid": "B177",
        "name": "Omega(c)(3120)",
        "flags": "B",
        "quantum_j": "5/2",
        "charge": 0.0,
    },
    {
        "mcid": 9304332,
        "parent_pdgid": "B209",
        "name": "Omega(c)(3185)",
        "flags": "B",
        "quantum_j": "1/2",
        "charge": 0.0,
    },
    {
        "mcid": 9304334,
        "parent_pdgid": "B210",
        "name": "Omega(c)(3327)",
        "flags": "B",
        "quantum_j": "3/2",
        "charge": 0.0,
    },
    {
        "mcid": 4422,
        "parent_pdgid": "S068",
        "name": "Xi(cc)++",
        "flags": "B",
        "quantum_j": "1/2",
        "charge": 2.0,
    },
)
EXTERNAL_COMMENT_NAMES = {
    22: "gamma",
    11: "e-",
    -11: "e+",
    13: "mu-",
    -13: "mu+",
    15: "tau-",
    -15: "tau+",
    12: "nu(e)",
    -12: "anti-nu(e)",
    14: "nu(mu)",
    -14: "anti-nu(mu)",
    16: "nu(tau)",
    -16: "anti-nu(tau)",
}


def parse_particle_list(path: Path) -> List[ParticleEntry]:
    entries: List[ParticleEntry] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 14:
                continue
            entries.append(
                ParticleEntry(
                    pdgid=int(parts[0]),
                    name=parts[1],
                    stable=int(parts[2]),
                    mass=float(parts[3]),
                    degeneracy=float(parts[4]),
                    statistics=int(parts[5]),
                    baryon=int(parts[6]),
                    charge=int(parts[7]),
                    strange=int(parts[8]),
                    charm=int(parts[9]),
                    abs_s=float(parts[10]),
                    abs_c=float(parts[11]),
                    width=float(parts[12]),
                    threshold=float(parts[13]),
                )
            )
    return entries


def parse_decay_file(path: Path) -> Dict[int, List[DecayChannel]]:
    stripped: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.split("#", 1)[0].strip()
            if line:
                stripped.append(line)

    decays: Dict[int, List[DecayChannel]] = {}
    pos = 0
    while pos < len(stripped):
        parent = int(stripped[pos])
        pos += 1
        count = int(stripped[pos])
        pos += 1
        channels: List[DecayChannel] = []
        for _ in range(count):
            fields = stripped[pos].split()
            pos += 1
            channels.append(
                DecayChannel(
                    br=float(fields[0]),
                    daughters=[int(item) for item in fields[1:]],
                )
            )
        decays[parent] = channels
    return decays


def parse_mcd_established_abs_ids(path: Path) -> set[int]:
    """Return abs(PDGID) values present in the PDG MC mass-width table."""
    ids: set[int] = set()
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw or raw.startswith("*"):
                continue
            id_field = raw[:32]
            for token in id_field.split():
                try:
                    value = int(token)
                except ValueError:
                    continue
                if value != 0:
                    ids.add(abs(value))
    return ids


def parse_mcd_mass_width_by_abs_id(path: Path) -> Dict[int, Tuple[Optional[float], Optional[float]]]:
    """Return abs(PDGID) -> (mass_GeV, width_GeV) parsed from mass_width_2025.mcd."""
    data: Dict[int, Tuple[Optional[float], Optional[float]]] = {}
    num_pattern = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)")

    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw or raw.startswith("*"):
                continue

            ids: List[int] = []
            id_field = raw[:32]
            for token in id_field.split():
                try:
                    value = int(token)
                except ValueError:
                    continue
                if value != 0:
                    ids.append(abs(value))
            if not ids:
                continue

            tail = raw[32:]
            tokens = num_pattern.findall(tail)
            if not tokens:
                continue

            mass: Optional[float] = None
            width: Optional[float] = None
            try:
                mass = float(tokens[0].replace("D", "E").replace("d", "E"))
            except ValueError:
                mass = None
            if len(tokens) >= 4:
                try:
                    width = float(tokens[3].replace("D", "E").replace("d", "E"))
                except ValueError:
                    width = None

            for pid in ids:
                old_mass, old_width = data.get(pid, (None, None))
                use_mass = old_mass if old_mass is not None else mass
                use_width = old_width if old_width is not None else width
                data[pid] = (use_mass, use_width)

    return data


def parse_established_pdg_nodes(path: Path) -> set[str]:
    nodes: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            nodes.add(line)
    return nodes


def format_float(value: float, precision: int = 8) -> str:
    return f"{value:.{precision}g}"


def is_hidden_charm_name(name: str) -> bool:
    lname = name.lower()
    return (
        "(c)" in lname
        or lname.startswith("j/psi")
        or lname.startswith("psi(")
        or "chi(c" in lname
        or "eta(c" in lname
        or "h(c" in lname
    )


def is_charm_sector_entry(entry: ParticleEntry) -> bool:
    return entry.charm != 0 or is_hidden_charm_name(entry.name)


def pdg2020_style_sector_index(entry: ParticleEntry) -> int:
    # Keep non-charm hadrons first (mesons then baryons), and charm sector at end.
    if not is_charm_sector_entry(entry):
        if entry.baryon == 0:
            abs_s = abs(entry.strange)
            if abs_s == 0:
                return 0
            if abs_s == 1:
                return 1
            return 2
        abs_s = abs(entry.strange)
        if abs_s == 0:
            return 3
        if abs_s == 1:
            return 4
        if abs_s == 2:
            return 5
        return 6

    if entry.baryon == 0:
        # Hidden charm (charmonium-like) before open charm mesons.
        if entry.charm == 0:
            return 7
        if abs(entry.strange) == 0:
            return 8
        return 9

    # Charmed baryons, grouped by strangeness; double-charm at the very end.
    if abs(entry.charm) >= 2:
        return 13
    abs_s = abs(entry.strange)
    if abs_s == 0:
        return 10
    if abs_s == 1:
        return 11
    return 12


def charge_sort_rank(charge: int) -> int:
    rank_map = {
        -2: 0,
        -1: 1,
        0: 2,
        1: 3,
        2: 4,
    }
    return rank_map.get(charge, 10 + charge)


def family_sort_key(name: str) -> str:
    return normalize_name(strip_charge_suffix(name)).lower()


def sorted_entries_and_ids(
    entries: Sequence[ParticleEntry],
    pdg2020_rank_map: Optional[Dict[int, int]] = None,
) -> Tuple[List[ParticleEntry], List[int]]:
    rank_map = pdg2020_rank_map or {}
    sectors: Dict[int, int] = {}
    families: Dict[int, str] = {}
    for entry in entries:
        sectors[entry.pdgid] = pdg2020_style_sector_index(entry)
        families[entry.pdgid] = family_sort_key(entry.name)

    family_anchor_mass: Dict[Tuple[int, str], float] = {}
    family_known_rank: Dict[Tuple[int, str], int] = {}
    for entry in entries:
        key = (sectors[entry.pdgid], families[entry.pdgid])
        family_anchor_mass[key] = min(family_anchor_mass.get(key, 1.0e9), entry.mass)
        if entry.pdgid in rank_map:
            family_known_rank[key] = min(family_known_rank.get(key, 10**9), rank_map[entry.pdgid])

    def sort_key(entry: ParticleEntry) -> Tuple[int, int, int, float, str, float, int, int]:
        sector = sectors[entry.pdgid]
        family = families[entry.pdgid]
        fam_key = (sector, family)
        reference_rank = rank_map.get(entry.pdgid, family_known_rank.get(fam_key, 10**9))
        exact_rank_flag = 0 if entry.pdgid in rank_map else 1
        return (
            sector,
            reference_rank,
            exact_rank_flag,
            family_anchor_mass.get(fam_key, entry.mass),
            family,
            entry.mass,
            charge_sort_rank(entry.charge),
            entry.pdgid,
        )

    ordered = sorted(entries, key=sort_key)
    return ordered, [entry.pdgid for entry in ordered]


def antiparticle_comment_name(name: str) -> str:
    if name.startswith("anti-"):
        return name[5:]
    if name.endswith("++"):
        return f"{name[:-2]}--"
    if name.endswith("--"):
        return f"{name[:-2]}++"
    if name.endswith("+"):
        return f"{name[:-1]}-"
    if name.endswith("-"):
        return f"{name[:-1]}+"
    return f"anti-{name}"


def decay_comment_name(pid: int, names: Dict[int, str]) -> str:
    if pid in names:
        return names[pid]
    if pid in EXTERNAL_COMMENT_NAMES:
        return EXTERNAL_COMMENT_NAMES[pid]
    if -pid in names:
        return antiparticle_comment_name(names[-pid])
    if -pid in EXTERNAL_COMMENT_NAMES:
        return antiparticle_comment_name(EXTERNAL_COMMENT_NAMES[-pid])
    return str(pid)


def write_with_comment(handle, lhs: str, comment: str, comment_col: int = 36) -> None:
    if len(lhs) < comment_col:
        lhs = lhs.ljust(comment_col)
    handle.write(f"{lhs} # {comment}\n")


def write_particle_list(
    path: Path,
    entries: Iterable[ParticleEntry],
    source_text: str,
    description_line: str = "PDG2025 list containing strange and non-strange hadrons only",
) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"# {description_line}\n")
        handle.write(f"# Source: {source_text}\n")
        handle.write(
            "#         pdgid                name         stable      mass[GeV]     degeneracy"
            "     statistics              B              Q              S              C"
            "            |S|            |C|     width[GeV] threshold[GeV]\n"
        )
        for entry in entries:
            handle.write(
                f"{entry.pdgid:15d} {entry.name:>20s} {entry.stable:14d}"
                f" {format_float(entry.mass):>14s} {format_float(entry.degeneracy):>14s}"
                f" {entry.statistics:14d} {entry.baryon:14d} {entry.charge:14d}"
                f" {entry.strange:14d} {entry.charm:14d} {format_float(entry.abs_s):>14s}"
                f" {format_float(entry.abs_c):>14s} {format_float(entry.width):>14s}"
                f" {format_float(entry.threshold):>14s}\n"
            )


def write_decay_file(
    path: Path,
    decays: Dict[int, List[DecayChannel]],
    ordered_ids: List[int],
    names: Dict[int, str],
    source_text: str,
) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# the list of decays\n")
        handle.write(f"# Source: {source_text}\n")
        handle.write("# each entry consists of the following:\n")
        handle.write("# a line with the pdgid of decaying particle\n")
        handle.write("# a line with the number of decay channels\n")
        handle.write("# for each channel a line containing whitespace-separated values of the channel branching ratio and pdg ids of the daughter products\n")
        handle.write("# everything after the # symbol is treated as a comment and ignored\n")
        handle.write("# decays of antiparticles are not listed but generated from the listed decays of particles\n\n")

        for parent in ordered_ids:
            channels = decays.get(parent)
            if not channels:
                continue
            pname = names.get(parent, str(parent))
            write_with_comment(handle, str(parent), pname)
            write_with_comment(handle, str(len(channels)), f"{len(channels)} decay channels")
            for channel in channels:
                daughters = " ".join(str(pid) for pid in channel.daughters)
                dnames = " + ".join(decay_comment_name(pid, names) for pid in channel.daughters)
                lhs = f"{format_float(channel.br, precision=10):<15} {daughters}"
                write_with_comment(handle, lhs, f"{pname} -> {dnames}")
            handle.write("\n")


def unit_scale_to_gev(unit_text: str) -> Optional[float]:
    unit = unit_text.strip()
    if unit in ("", "GeV", "GeV/c^2", "GeV/c", "GeV/c^{2}", "GeV/c^{2} "):
        return 1.0
    if "MeV" in unit:
        return 1.0e-3
    if "keV" in unit:
        return 1.0e-6
    if "eV" in unit:
        return 1.0e-9
    return None


def pick_summary_value(cur: sqlite3.Cursor, pdgid: str) -> Optional[Tuple[float, str]]:
    rows = cur.execute(
        """
        SELECT value, unit_text, value_type, sort
        FROM pdgdata
        WHERE pdgid = ?
          AND in_summary_table = 1
          AND limit_type IS NULL
          AND value IS NOT NULL
        ORDER BY
          CASE value_type
            WHEN 'AC' THEN 0
            WHEN 'FC' THEN 1
            WHEN 'V'  THEN 2
            WHEN 'D'  THEN 3
            WHEN 'E'  THEN 4
            WHEN 'O'  THEN 5
            ELSE 9
          END,
          sort
        """,
        (pdgid,),
    ).fetchall()
    if not rows:
        return None
    return (float(rows[0]["value"]), str(rows[0]["unit_text"]))


def pick_mass_or_width(cur: sqlite3.Cursor, parent_pdgid: str, kind: str) -> Optional[float]:
    assert kind in ("M", "G")
    candidates = cur.execute(
        """
        SELECT pdgid
        FROM pdgid
        WHERE parent_pdgid = ?
          AND data_type = ?
        ORDER BY CASE WHEN instr(flags, 'D') > 0 THEN 0 ELSE 1 END, sort
        """,
        (parent_pdgid, kind),
    ).fetchall()
    for row in candidates:
        picked = pick_summary_value(cur, row["pdgid"])
        if picked is None:
            continue
        value, unit = picked
        scale = unit_scale_to_gev(unit)
        if scale is not None:
            return value * scale

    keyword = "%MASS%" if kind == "M" else "%WIDTH%"
    sec_candidates = cur.execute(
        """
        SELECT pdgid
        FROM pdgid
        WHERE parent_pdgid = ?
          AND data_type = 'SEC'
          AND UPPER(description) LIKE ?
        ORDER BY sort
        """,
        (parent_pdgid, keyword),
    ).fetchall()
    for row in sec_candidates:
        picked = pick_summary_value(cur, row["pdgid"])
        if picked is None:
            continue
        value, unit = picked
        scale = unit_scale_to_gev(unit)
        if scale is not None:
            return value * scale
    return None


def pick_mass_or_width_upper_limit(cur: sqlite3.Cursor, parent_pdgid: str, kind: str) -> Optional[float]:
    assert kind in ("M", "G")
    candidates = cur.execute(
        """
        SELECT pdgid
        FROM pdgid
        WHERE parent_pdgid = ?
          AND data_type = ?
        ORDER BY CASE WHEN instr(flags, 'D') > 0 THEN 0 ELSE 1 END, sort
        """,
        (parent_pdgid, kind),
    ).fetchall()

    keyword = "%MASS%" if kind == "M" else "%WIDTH%"
    sec_candidates = cur.execute(
        """
        SELECT pdgid
        FROM pdgid
        WHERE parent_pdgid = ?
          AND data_type = 'SEC'
          AND UPPER(description) LIKE ?
        ORDER BY sort
        """,
        (parent_pdgid, keyword),
    ).fetchall()
    candidates = list(candidates) + list(sec_candidates)

    limits: List[float] = []
    for row in candidates:
        rows = cur.execute(
            """
            SELECT value, unit_text
            FROM pdgdata
            WHERE pdgid = ?
              AND in_summary_table = 1
              AND limit_type = 'U'
              AND value IS NOT NULL
            ORDER BY sort
            """,
            (str(row["pdgid"]),),
        ).fetchall()
        for item in rows:
            scale = unit_scale_to_gev(str(item["unit_text"]))
            if scale is None:
                continue
            value = float(item["value"]) * scale
            if value >= 0.0:
                limits.append(value)

    if not limits:
        return None
    return min(limits)


def fetch_parent_map(cur: sqlite3.Cursor) -> Dict[int, str]:
    rows = cur.execute(
        """
        SELECT ABS(p.mcid) AS amcid, p.mcid, p.pdgid, p.cc_type
        FROM pdgparticle AS p
        JOIN pdgid AS g ON g.pdgid = p.pdgid
        WHERE p.mcid IS NOT NULL
          AND g.data_type = 'PART'
        ORDER BY
          ABS(p.mcid),
          CASE
            WHEN p.mcid > 0 THEN 0
            WHEN p.cc_type = 'S' THEN 1
            ELSE 2
          END
        """
    ).fetchall()
    ret: Dict[int, str] = {}
    for row in rows:
        amcid = int(row["amcid"])
        if amcid not in ret:
            ret[amcid] = str(row["pdgid"])
    return ret


def resolve_decay_daughter_mcid(cur: sqlite3.Cursor, pdgitem_id: int, name: str) -> Optional[int]:
    direct = cur.execute(
        """
        SELECT mcid
        FROM pdgparticle
        WHERE pdgitem_id = ?
          AND name = ?
          AND mcid IS NOT NULL
        """,
        (pdgitem_id, name),
    ).fetchall()
    if len(direct) == 1:
        return int(direct[0]["mcid"])

    fallback = cur.execute(
        """
        SELECT DISTINCT mcid
        FROM pdgparticle
        WHERE pdgitem_id = ?
          AND mcid IS NOT NULL
        """,
        (pdgitem_id,),
    ).fetchall()
    if len(fallback) == 1:
        return int(fallback[0]["mcid"])
    return None


def fetch_explicit_bfx_modes(cur: sqlite3.Cursor, parent_pdgid: str) -> List[DecayChannel]:
    modes = cur.execute(
        """
        SELECT pdgid
        FROM pdgid
        WHERE parent_pdgid = ?
          AND data_type LIKE 'BFX%'
        ORDER BY sort
        """,
        (parent_pdgid,),
    ).fetchall()

    channels: List[DecayChannel] = []
    for mode_row in modes:
        mode_pdgid = str(mode_row["pdgid"])
        picked = pick_summary_value(cur, mode_pdgid)
        if picked is None:
            continue
        br_value, _ = picked
        if br_value <= 0.0:
            continue

        daughters_rows = cur.execute(
            """
            SELECT pdgitem_id, name, multiplier
            FROM pdgdecay
            WHERE pdgid = ?
              AND is_outgoing = 1
            ORDER BY sort
            """,
            (mode_pdgid,),
        ).fetchall()

        daughters: List[int] = []
        unresolved = False
        for daughter_row in daughters_rows:
            daughter_mcid = resolve_decay_daughter_mcid(
                cur,
                int(daughter_row["pdgitem_id"]),
                str(daughter_row["name"]),
            )
            if daughter_mcid is None:
                unresolved = True
                break
            daughters.extend([daughter_mcid] * max(int(daughter_row["multiplier"]), 1))

        if unresolved or not daughters:
            continue
        channels.append(DecayChannel(br=br_value, daughters=daughters))

    return channels


def channel_key(channel: DecayChannel) -> Tuple[int, ...]:
    return tuple(sorted(channel.daughters))


def aggregate_channels(channels: Sequence[DecayChannel]) -> List[DecayChannel]:
    merged: Dict[Tuple[int, ...], float] = {}
    for ch in channels:
        key = channel_key(ch)
        merged[key] = merged.get(key, 0.0) + ch.br
    return [DecayChannel(br=br, daughters=list(key)) for key, br in sorted(merged.items())]


def normalize_channels(channels: List[DecayChannel]) -> bool:
    total = sum(max(ch.br, 0.0) for ch in channels)
    if total <= 0.0:
        return False
    for ch in channels:
        ch.br = max(ch.br, 0.0) / total
    return True


def channels_with_known_daughters(channels: Sequence[DecayChannel], known_ids: Sequence[int]) -> List[DecayChannel]:
    known = {abs(pid) for pid in known_ids}
    out: List[DecayChannel] = []
    for channel in channels:
        if all(abs(pid) in known or abs(pid) in EXTERNAL_MASSES_GEV for pid in channel.daughters):
            out.append(DecayChannel(br=channel.br, daughters=list(channel.daughters)))
    return out


def signed_qnums(pid: int, by_absid: Dict[int, ParticleEntry]) -> Optional[Tuple[int, int, int, int]]:
    ad = abs(pid)
    if ad in EXTERNAL_MASSES_GEV:
        if ad in (11, 13, 15):
            q = -1
        else:
            q = 0
        if pid < 0:
            q = -q
        return (0, q, 0, 0)
    if ad not in by_absid:
        return None
    base = by_absid[ad]
    if pid < 0:
        return (-base.baryon, -base.charge, -base.strange, -base.charm)
    return (base.baryon, base.charge, base.strange, base.charm)


def filter_channels_for_parent(
    parent: ParticleEntry,
    channels: Sequence[DecayChannel],
    by_absid: Dict[int, ParticleEntry],
) -> List[DecayChannel]:
    out: List[DecayChannel] = []
    pchg = (parent.baryon, parent.charge, parent.strange, parent.charm)
    for channel in channels:
        total = [0, 0, 0, 0]
        valid = True
        for daughter in channel.daughters:
            dchg = signed_qnums(daughter, by_absid)
            if dchg is None:
                valid = False
                break
            for idx in range(4):
                total[idx] += dchg[idx]
        if not valid:
            continue
        # For unstable resonances enforce exact B,Q,S,C conservation.
        if parent.stable == 0 and tuple(total) != pchg:
            continue
        # For weakly decaying stable hadrons allow S/C changes.
        # Still require B and Q conservation for sanity.
        if parent.stable == 1:
            if total[0] != pchg[0] or total[1] != pchg[1]:
                continue
        out.append(DecayChannel(br=channel.br, daughters=list(channel.daughters)))
    return out


def update_channels_from_pdg(base_channels: List[DecayChannel], pdg_channels: List[DecayChannel]) -> Tuple[List[DecayChannel], bool]:
    base_agg = aggregate_channels(base_channels)
    pdg_agg = aggregate_channels(pdg_channels)
    base_by_key: Dict[Tuple[int, ...], float] = {channel_key(ch): ch.br for ch in base_agg}
    pdg_by_key: Dict[Tuple[int, ...], float] = {channel_key(ch): ch.br for ch in pdg_agg}

    matched_keys = set(base_by_key.keys()) & set(pdg_by_key.keys())
    if not matched_keys:
        return base_channels, False
    new_keys = set(pdg_by_key.keys()) - set(base_by_key.keys())

    old_known_sum = sum(base_by_key[key] for key in matched_keys)
    new_known_sum = sum(pdg_by_key[key] for key in matched_keys)
    new_only_sum = sum(pdg_by_key[key] for key in new_keys)
    new_total_known = new_known_sum + new_only_sum

    # Guard against large inflation from sparse/ratio-only partial BR updates.
    # If there are no newly mapped modes and explicit BR coverage is poor,
    # keep the existing calibrated block.
    if not new_keys and old_known_sum > 0.0 and (new_known_sum / old_known_sum) < 0.6:
        return base_channels, False

    if new_total_known >= 1.0 + 1.0e-6:
        return base_channels, False

    old_rest_sum = 1.0 - old_known_sum
    new_rest_target = max(0.0, 1.0 - new_total_known)
    scale = 0.0
    if old_rest_sum > 1.0e-12:
        scale = new_rest_target / old_rest_sum

    updated: List[DecayChannel] = []
    for ch in base_agg:
        key = channel_key(ch)
        if key in matched_keys:
            updated.append(DecayChannel(br=pdg_by_key[key], daughters=list(key)))
        else:
            updated.append(DecayChannel(br=(ch.br * scale if old_rest_sum > 1.0e-12 else 0.0), daughters=list(key)))
    for key in sorted(new_keys):
        updated.append(DecayChannel(br=pdg_by_key[key], daughters=list(key)))

    if not normalize_channels(updated):
        return base_channels, False
    return updated, True


def merge_pdg_and_template(pdg_channels: List[DecayChannel], template_channels: List[DecayChannel]) -> List[DecayChannel]:
    pdg = aggregate_channels(pdg_channels)
    if not pdg and not template_channels:
        return []
    if pdg and not template_channels:
        normalize_channels(pdg)
        return pdg
    if not pdg and template_channels:
        ret = [DecayChannel(br=ch.br, daughters=list(ch.daughters)) for ch in template_channels]
        normalize_channels(ret)
        return ret

    pdg_sum = sum(ch.br for ch in pdg)
    if pdg_sum >= 1.0 - 1.0e-9:
        normalize_channels(pdg)
        return pdg

    pdg_keys = {channel_key(ch) for ch in pdg}
    extras = [DecayChannel(br=ch.br, daughters=list(ch.daughters)) for ch in template_channels if channel_key(ch) not in pdg_keys]
    if extras:
        normalize_channels(extras)
        for ch in extras:
            ch.br *= max(0.0, 1.0 - pdg_sum)
        ret = [DecayChannel(br=ch.br, daughters=list(ch.daughters)) for ch in pdg] + extras
        normalize_channels(ret)
        return ret

    normalize_channels(pdg)
    return pdg


def abs_charm_content_from_mcid(mcid: int, flags: str) -> int:
    code = abs(int(mcid))
    if flags == "M":
        q1 = (code // 100) % 10
        q2 = (code // 10) % 10
        return int((1 if q1 == 4 else 0) + (1 if q2 == 4 else 0))
    if flags == "B":
        q1 = (code // 1000) % 10
        q2 = (code // 100) % 10
        q3 = (code // 10) % 10
        return int((1 if q1 == 4 else 0) + (1 if q2 == 4 else 0) + (1 if q3 == 4 else 0))
    return 0


def read_pdg_source_text(cur: sqlite3.Cursor) -> str:
    release = cur.execute("SELECT value FROM pdginfo WHERE name='data_release_timestamp'").fetchone()
    citation = cur.execute("SELECT value FROM pdginfo WHERE name='citation'").fetchone()
    release_value = str(release["value"]) if release is not None else "unknown release timestamp"
    citation_value = str(citation["value"]) if citation is not None else "PDG API citation unavailable"
    return f"PDG API 2025 SQLite ({release_value}); {citation_value}"


def parse_j_value(quantum_j: Optional[str]) -> Optional[float]:
    if quantum_j is None:
        return None
    text = quantum_j.strip()
    if not text:
        return None

    frac = re.search(r"(-?\d+)\s*/\s*(\d+)", text)
    if frac:
        num = int(frac.group(1))
        den = int(frac.group(2))
        if den != 0:
            return num / den

    num = re.search(r"(-?\d+(\.\d+)?)", text)
    if num:
        return float(num.group(1))
    return None


def normalize_name(name: str) -> str:
    ret = name.strip().replace(" ", "")
    ret = ret.replace("^", "").replace("_", "")
    ret = ret.replace("{", "").replace("}", "")
    ret = ret.replace("\\", "")
    # Keep D-meson spin notation consistent with PDG2020 style.
    ret = re.sub(r"^D([0-9])\*\(", r"D(\1)*(", ret)
    ret = re.sub(r"^D([0-9])\(", r"D(\1)(", ret)
    # Keep kaon spin notation consistent with PDG2020/mass_width tables:
    # K0*(700) -> K(0)*(700), K2*(1980) -> K(2)*(1980), etc.
    ret = re.sub(r"^K([0-9]+)\*", r"K(\1)*", ret)
    # Keep isoscalar f-meson labels consistent with PDG2020 style.
    ret = re.sub(r"^f([0-9])\(([^)]+)\)0$", r"f(\1)(\2)", ret)
    return ret


def normalize_display_name(name: str, mcid: int) -> str:
    if mcid in DISPLAY_NAME_OVERRIDES:
        return DISPLAY_NAME_OVERRIDES[mcid]
    amcid = abs(mcid)
    if amcid in DISPLAY_NAME_OVERRIDES:
        return DISPLAY_NAME_OVERRIDES[amcid]
    return normalize_name(name)


def infer_mass_from_name_gev(name: str) -> Optional[float]:
    matches = re.findall(r"\((\d+(?:\.\d+)?)\)", name)
    if not matches:
        return None
    value = float(matches[-1])
    if value > 20.0:
        return value * 1.0e-3
    return value


def strip_charge_suffix(name: str) -> str:
    text = name.strip()
    for suffix in ("++", "--", "+", "-", "0"):
        if text.endswith(suffix) and len(text) > len(suffix):
            return text[: -len(suffix)]
    return text


def name_match_keys(name: str) -> set[str]:
    keys: set[str] = set()
    normalized = normalize_name(name).lower()
    if normalized:
        keys.add(normalized)
    stripped = normalize_name(strip_charge_suffix(name)).lower()
    if stripped:
        keys.add(stripped)
    return keys


def fetch_part_name_keys(cur: sqlite3.Cursor) -> set[str]:
    keys: set[str] = set()

    part_descriptions = cur.execute(
        """
        SELECT DISTINCT description
        FROM pdgid
        WHERE data_type = 'PART'
          AND description IS NOT NULL
        """
    ).fetchall()
    for row in part_descriptions:
        keys.update(name_match_keys(str(row["description"])))

    part_particle_names = cur.execute(
        """
        SELECT DISTINCT p.name
        FROM pdgparticle AS p
        JOIN pdgid AS g ON g.pdgid = p.pdgid
        WHERE g.data_type = 'PART'
          AND p.name IS NOT NULL
        """
    ).fetchall()
    for row in part_particle_names:
        keys.update(name_match_keys(str(row["name"])))

    return keys


def fetch_part_parent_map_by_name(cur: sqlite3.Cursor) -> Dict[str, str]:
    """Map normalized name keys to unique PART parent PDG IDs."""
    by_key: Dict[str, set[str]] = {}

    rows_desc = cur.execute(
        """
        SELECT DISTINCT pdgid, description
        FROM pdgid
        WHERE data_type = 'PART'
          AND description IS NOT NULL
        """
    ).fetchall()
    for row in rows_desc:
        parent = str(row["pdgid"])
        for key in name_match_keys(str(row["description"])):
            by_key.setdefault(key, set()).add(parent)

    rows_name = cur.execute(
        """
        SELECT DISTINCT p.pdgid, p.name
        FROM pdgparticle AS p
        JOIN pdgid AS g ON g.pdgid = p.pdgid
        WHERE g.data_type = 'PART'
          AND p.name IS NOT NULL
        """
    ).fetchall()
    for row in rows_name:
        parent = str(row["pdgid"])
        for key in name_match_keys(str(row["name"])):
            by_key.setdefault(key, set()).add(parent)

    ret: Dict[str, str] = {}
    for key, parents in by_key.items():
        if len(parents) == 1:
            ret[key] = next(iter(parents))
    return ret


def fetch_node_part_name_keys(cur: sqlite3.Cursor, nodes: Sequence[str]) -> set[str]:
    if not nodes:
        return set()
    placeholders = ",".join("?" for _ in nodes)
    params = tuple(nodes)
    keys: set[str] = set()

    rows_desc = cur.execute(
        f"""
        SELECT DISTINCT description
        FROM pdgid
        WHERE data_type = 'PART'
          AND pdgid IN ({placeholders})
          AND description IS NOT NULL
        """,
        params,
    ).fetchall()
    for row in rows_desc:
        keys.update(name_match_keys(str(row["description"])))

    rows_name = cur.execute(
        f"""
        SELECT DISTINCT name
        FROM pdgparticle
        WHERE pdgid IN ({placeholders})
          AND name IS NOT NULL
        """,
        params,
    ).fetchall()
    for row in rows_name:
        keys.update(name_match_keys(str(row["name"])))

    return keys


def fetch_node_abs_mcids(cur: sqlite3.Cursor, nodes: Sequence[str]) -> set[int]:
    if not nodes:
        return set()
    placeholders = ",".join("?" for _ in nodes)
    params = tuple(nodes)
    rows = cur.execute(
        f"""
        SELECT DISTINCT ABS(mcid) AS amcid
        FROM pdgparticle
        WHERE pdgid IN ({placeholders})
          AND mcid IS NOT NULL
        """,
        params,
    ).fetchall()
    return {int(row["amcid"]) for row in rows}


def quark_charge(qdigit: int) -> float:
    if qdigit in (1, 3, 5):
        return -1.0 / 3.0
    if qdigit in (2, 4, 6):
        return 2.0 / 3.0
    return 0.0


def meson_quantum_assignment(quark: int, antiquark: int) -> Tuple[float, int, int]:
    charge = quark_charge(quark) - quark_charge(antiquark)
    strange = (-1 if quark == 3 else 0) + (1 if antiquark == 3 else 0)
    charm = (1 if quark == 4 else 0) + (-1 if antiquark == 4 else 0)
    return charge, strange, charm


def quantum_numbers_from_mcid(
    mcid: int,
    flags: str,
    charge_hint: Optional[float] = None,
    name_hint: str = "",
) -> Optional[Tuple[int, int, int, int]]:
    sign = 1 if mcid >= 0 else -1
    code = abs(mcid)

    if flags == "M":
        q1 = (code // 100) % 10
        q2 = (code // 10) % 10
        if q1 == 0:
            return (0, 0, 0, 0)

        # Two possible conventions exist for q1/q2 ordering in meson MCIDs.
        # Use measured electric charge when available to resolve orientation.
        c1, s1, ch1 = meson_quantum_assignment(q1, q2)
        c2, s2, ch2 = meson_quantum_assignment(q2, q1)
        options = [(c1, s1, ch1), (c2, s2, ch2)]

        if charge_hint is not None:
            matched = [item for item in options if abs(item[0] - charge_hint) < 1.0e-6]
            if len(matched) == 1:
                charge, strange, charm = matched[0]
            elif len(matched) == 2:
                lname = name_hint.lower()
                if "bar" in lname:
                    charge, strange, charm = sorted(matched, key=lambda item: (item[1], item[2]))[0]
                else:
                    charge, strange, charm = sorted(matched, key=lambda item: (-item[1], -item[2]))[0]
            else:
                charge, strange, charm = min(options, key=lambda item: abs(item[0] - charge_hint))
        else:
            charge, strange, charm = options[0]

        if sign < 0:
            charge, strange, charm = -charge, -strange, -charm
        return (0, int(round(charge)), int(round(strange)), int(round(charm)))

    if flags == "B":
        q1 = (code // 1000) % 10
        q2 = (code // 100) % 10
        q3 = (code // 10) % 10
        charge = quark_charge(q1) + quark_charge(q2) + quark_charge(q3)
        strange = sum(-1 for q in (q1, q2, q3) if q == 3)
        charm = sum(1 for q in (q1, q2, q3) if q == 4)
        baryon = 1
        if sign < 0:
            baryon, charge, strange, charm = -baryon, -charge, -strange, -charm
        return (baryon, int(round(charge)), int(round(strange)), int(round(charm)))

    return None


def infer_degeneracy(mcid: int, quantum_j: Optional[str]) -> float:
    j = parse_j_value(quantum_j)
    if j is not None and j >= 0.0:
        val = 2.0 * j + 1.0
        if abs(val - round(val)) < 1.0e-7:
            return float(int(round(val)))
        return val
    j2p1 = abs(mcid) % 10
    if j2p1 > 0:
        return float(j2p1)
    return 1.0


def fetch_expansion_rows(cur: sqlite3.Cursor, existing_abs_ids: Sequence[int]) -> List[sqlite3.Row]:
    rows = cur.execute(
        """
        SELECT p.mcid, p.pdgid, p.name, p.cc_type, p.charge, p.quantum_j, g.flags, g.description
        FROM pdgparticle AS p
        JOIN pdgid AS g ON g.pdgid = p.pdgid
        WHERE p.mcid IS NOT NULL
          AND g.data_type = 'PART'
          AND g.flags IN ('M', 'B')
          AND (
            (g.flags='M' AND ((ABS(p.mcid)/100)%10)<=4 AND ((ABS(p.mcid)/10)%10)<=4)
            OR
            (g.flags='B' AND ((ABS(p.mcid)/1000)%10)<=4 AND ((ABS(p.mcid)/100)%10)<=4 AND ((ABS(p.mcid)/10)%10)<=4)
          )
        ORDER BY ABS(p.mcid), p.cc_type, p.charge DESC, p.mcid DESC
        """
    ).fetchall()

    by_abs: Dict[int, List[sqlite3.Row]] = {}
    for row in rows:
        amcid = abs(int(row["mcid"]))
        by_abs.setdefault(amcid, []).append(row)

    existing = {abs(int(pid)) for pid in existing_abs_ids}
    selected: List[sqlite3.Row] = []
    for amcid, group in sorted(by_abs.items()):
        if amcid in existing or amcid in EXCLUDED_EXPANSION_MCIDS:
            continue

        def rank(row: sqlite3.Row) -> Tuple[int, int, float]:
            flags = str(row["flags"])
            cc = str(row["cc_type"] or "")
            charge = float(row["charge"] or 0.0)
            if flags == "M":
                # prefer particle/self-conjugate and non-negative charge for mesons
                return (
                    0 if cc in ("P", "S") else 1,
                    0 if charge >= -1.0e-9 else 1,
                    -charge,
                )
            # baryons: keep particle states (including negative electric charge baryons)
            return (
                0 if cc in ("P", "S") else 1,
                0,
                -charge,
            )

        selected.append(sorted(group, key=rank)[0])

    return selected


def build_entry_from_row(
    cur: sqlite3.Cursor,
    row: sqlite3.Row,
    mcd_mass_width_by_absid: Optional[Dict[int, Tuple[Optional[float], Optional[float]]]] = None,
) -> Optional[ParticleEntry]:
    mcid = int(row["mcid"])
    parent = str(row["pdgid"])
    flags = str(row["flags"])
    qnums = quantum_numbers_from_mcid(
        mcid,
        flags,
        charge_hint=float(row["charge"]) if row["charge"] is not None else None,
        name_hint=str(row["name"]),
    )
    if qnums is None:
        return None
    baryon, charge, strange, charm = qnums

    mass = None
    if mcd_mass_width_by_absid is not None:
        mcd_mass = mcd_mass_width_by_absid.get(abs(mcid), (None, None))[0]
        if mcd_mass is not None and mcd_mass > 0.0:
            mass = mcd_mass
    if mass is None or mass <= 0.0:
        mass = pick_mass_or_width(cur, parent, "M")
    if mass is None or mass <= 0.0:
        return None

    width = None
    if mcd_mass_width_by_absid is not None:
        mcd_width = mcd_mass_width_by_absid.get(abs(mcid), (None, None))[1]
        if mcd_width is not None and mcd_width >= 0.0:
            width = mcd_width
    if width is None or width < 0.0:
        width = pick_mass_or_width(cur, parent, "G")
    if width is None or width < 0.0:
        width = 0.0
    width_ul = pick_mass_or_width_upper_limit(cur, parent, "G")
    if width_ul is not None and width > width_ul:
        width = width_ul

    name = normalize_display_name(str(row["name"]), mcid)
    degeneracy = infer_degeneracy(mcid, row["quantum_j"])
    statistics = -1 if flags == "M" else 1
    abs_c = abs_charm_content_from_mcid(mcid, flags)

    return ParticleEntry(
        pdgid=int(mcid),
        name=name,
        stable=0,
        mass=float(mass),
        degeneracy=float(degeneracy),
        statistics=statistics,
        baryon=baryon,
        charge=charge,
        strange=strange,
        charm=charm,
        abs_s=float(abs(strange)),
        abs_c=float(abs_c if abs_c > 0 else abs(charm)),
        width=float(width),
        threshold=0.0,
    )


def build_entry_from_manual_override(
    cur: sqlite3.Cursor,
    mcid: int,
    parent_pdgid: str,
    name: str,
    flags: str,
    quantum_j: Optional[str],
    charge: Optional[float],
) -> Optional[ParticleEntry]:
    qnums = quantum_numbers_from_mcid(
        mcid,
        flags,
        charge_hint=charge,
        name_hint=name,
    )
    if qnums is None:
        return None
    baryon, qcharge, strange, charm = qnums

    mass = pick_mass_or_width(cur, parent_pdgid, "M")
    if mass is None or mass <= 0.0:
        mass = infer_mass_from_name_gev(name)
        if mass is None or mass <= 0.0:
            return None
    width = pick_mass_or_width(cur, parent_pdgid, "G")
    if width is None or width < 0.0:
        width = 0.0
    width_ul = pick_mass_or_width_upper_limit(cur, parent_pdgid, "G")
    if width_ul is not None and width > width_ul:
        width = width_ul
    abs_c = abs_charm_content_from_mcid(mcid, flags)

    return ParticleEntry(
        pdgid=int(mcid),
        name=normalize_display_name(name, mcid),
        stable=0,
        mass=float(mass),
        degeneracy=infer_degeneracy(mcid, quantum_j),
        statistics=-1 if flags == "M" else 1,
        baryon=baryon,
        charge=qcharge,
        strange=strange,
        charm=charm,
        abs_s=float(abs(strange)),
        abs_c=float(abs_c if abs_c > 0 else abs(charm)),
        width=float(width),
        threshold=0.0,
    )


def select_template_decay(
    entry: ParticleEntry,
    entries_by_id: Dict[int, ParticleEntry],
    decays: Dict[int, List[DecayChannel]],
) -> Optional[List[DecayChannel]]:
    candidates: List[Tuple[float, int]] = []
    for pid, other in entries_by_id.items():
        if pid == entry.pdgid:
            continue
        channels = decays.get(pid)
        if not channels:
            continue
        if other.stable == 1:
            continue
        if other.statistics != entry.statistics:
            continue
        if (other.baryon, other.charge, other.strange, other.charm) != (
            entry.baryon,
            entry.charge,
            entry.strange,
            entry.charm,
        ):
            continue
        mass_diff = abs(other.mass - entry.mass)
        deg_penalty = abs(other.degeneracy - entry.degeneracy)
        score = mass_diff + 0.15 * deg_penalty
        candidates.append((score, pid))

    if not candidates:
        return None
    candidates.sort()
    template_pid = candidates[0][1]
    return [DecayChannel(br=ch.br, daughters=list(ch.daughters)) for ch in decays[template_pid]]


def fallback_generic_decay(entry: ParticleEntry, entries_by_id: Dict[int, ParticleEntry]) -> List[DecayChannel]:
    # Last-resort placeholder preserving charges using nearest lower-mass state.
    matches: List[Tuple[float, int]] = []
    for pid, other in entries_by_id.items():
        if pid == entry.pdgid:
            continue
        if (other.baryon, other.charge, other.strange, other.charm) == (
            entry.baryon,
            entry.charge,
            entry.strange,
            entry.charm,
        ) and other.mass < entry.mass:
            matches.append((entry.mass - other.mass, pid))
    if not matches:
        return []
    matches.sort()
    daughter = matches[0][1]
    pi0_mass = entries_by_id.get(111, ParticleEntry(111, "pi0", 1, 0.1349768, 1.0, -1, 0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0)).mass
    daughters = [daughter]
    if entries_by_id[daughter].mass + pi0_mass < entry.mass:
        daughters.append(111)
    return [DecayChannel(br=1.0, daughters=daughters)]


def compute_hadronic_strong_threshold(
    parent: ParticleEntry,
    channels: Sequence[DecayChannel],
    by_absid: Dict[int, ParticleEntry],
    mass_by_absid: Dict[int, float],
) -> float:
    parent_qnums = (parent.baryon, parent.charge, parent.strange, parent.charm)
    best: Optional[float] = None
    for channel in channels:
        total = 0.0
        total_qnums = [0, 0, 0, 0]
        valid = True
        for daughter in channel.daughters:
            ad = abs(daughter)
            # hadronic/strong thresholds: only hadronic daughters in-list
            if ad not in by_absid or ad not in mass_by_absid:
                valid = False
                break
            total += mass_by_absid[ad]
            dchg = signed_qnums(daughter, by_absid)
            if dchg is None:
                valid = False
                break
            for idx in range(4):
                total_qnums[idx] += dchg[idx]
        if valid and tuple(total_qnums) != parent_qnums:
            valid = False
        if valid:
            best = total if best is None else min(best, total)
    return 0.0 if best is None else best


def clone_particle_entries(entries: Sequence[ParticleEntry]) -> List[ParticleEntry]:
    return [
        ParticleEntry(
            pdgid=item.pdgid,
            name=item.name,
            stable=item.stable,
            mass=item.mass,
            degeneracy=item.degeneracy,
            statistics=item.statistics,
            baryon=item.baryon,
            charge=item.charge,
            strange=item.strange,
            charm=item.charm,
            abs_s=item.abs_s,
            abs_c=item.abs_c,
            width=item.width,
            threshold=item.threshold,
        )
        for item in entries
    ]


def clone_decay_map(decays: Dict[int, List[DecayChannel]]) -> Dict[int, List[DecayChannel]]:
    return {
        pid: [DecayChannel(br=channel.br, daughters=list(channel.daughters)) for channel in channels]
        for pid, channels in decays.items()
    }


def select_extra_entries(
    base_entries: Sequence[ParticleEntry],
    extended_entries: Sequence[ParticleEntry],
) -> List[ParticleEntry]:
    base_ids = {entry.pdgid for entry in base_entries}
    return [entry for entry in extended_entries if entry.pdgid not in base_ids]


def sanitize_decays_for_entries(
    entries: Sequence[ParticleEntry],
    ids_order: Sequence[int],
    decays: Dict[int, List[DecayChannel]],
    recompute_threshold_ids: Optional[set[int]] = None,
) -> Dict[int, List[DecayChannel]]:
    entries_by_id = {entry.pdgid: entry for entry in entries}
    known_ids = [entry.pdgid for entry in entries]
    by_absid = {abs(entry.pdgid): entry for entry in entries}
    mass_by_absid = {abs(entry.pdgid): entry.mass for entry in entries}

    sanitized: Dict[int, List[DecayChannel]] = {}
    for pid in ids_order:
        parent = entries_by_id.get(pid)
        if parent is None:
            continue
        recompute_threshold = parent.stable == 0 and (
            recompute_threshold_ids is None or pid in recompute_threshold_ids
        )

        channels = decays.get(pid)
        if not channels:
            if recompute_threshold:
                parent.threshold = 0.0
            continue

        kept = channels_with_known_daughters(channels, known_ids)
        kept = filter_channels_for_parent(parent, kept, by_absid)
        if not kept or not normalize_channels(kept):
            if recompute_threshold:
                parent.threshold = 0.0
            continue

        if recompute_threshold:
            parent.threshold = compute_hadronic_strong_threshold(parent, kept, by_absid, mass_by_absid)
        sanitized[pid] = kept

    return sanitized


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    base_dir = script_dir.parent
    listing_dir = base_dir / "pdglisting"
    pdg2020_dir = base_dir.parent / "PDG2020"
    db_path = listing_dir / "pdg-2025-v0.2.2.sqlite"
    mcd_path = listing_dir / "mass_width_2025.mcd"
    pdglive_nodes_path = listing_dir / "pdglive-established-nodes-2025.txt"
    out_list = base_dir / "list.dat"
    out_decays = base_dir / "decays.dat"
    out_list_with_charm = base_dir / "list-withcharm.dat"
    out_list_extra = base_dir / "list-extra.dat"
    out_list_with_nuclei = base_dir / "list-withnuclei.dat"
    out_list_with_excited_nuclei = base_dir / "list-withexcitednuclei.dat"

    pdg2020_base_entries = parse_particle_list(pdg2020_dir / "list-withcharm.dat")
    pdg2020_rank_map = {entry.pdgid: index for index, entry in enumerate(pdg2020_base_entries)}
    pdg2020_with_nuclei_entries = parse_particle_list(pdg2020_dir / "list-withnuclei-withcharm.dat")
    pdg2020_with_excited_nuclei_entries = parse_particle_list(pdg2020_dir / "list-withexcitednuclei.dat")
    pdg2025_ref_with_charm_path = base_dir.parent / "PDG2025" / "list-withcharm.dat"
    pdg2025_ref_with_charm_entries = (
        parse_particle_list(pdg2025_ref_with_charm_path)
        if pdg2025_ref_with_charm_path.exists()
        else []
    )
    pdg2025_ref_with_charm_by_id = {entry.pdgid: entry for entry in pdg2025_ref_with_charm_entries}
    pdg2025_ref_decays_path = base_dir.parent / "PDG2025" / "decays.dat"
    pdg2025_ref_decays = (
        parse_decay_file(pdg2025_ref_decays_path)
        if pdg2025_ref_decays_path.exists()
        else {}
    )
    pdg2020_nuclei_thresholds = {
        entry.pdgid: entry.threshold
        for entry in pdg2020_with_excited_nuclei_entries
        if abs(entry.pdgid) >= 1000000000
    }
    nuclei_entries = select_extra_entries(pdg2020_base_entries, pdg2020_with_nuclei_entries)
    excited_nuclei_entries = select_extra_entries(pdg2020_with_nuclei_entries, pdg2020_with_excited_nuclei_entries)

    entries = clone_particle_entries(pdg2020_base_entries)
    decays_pdg2020 = parse_decay_file(pdg2020_dir / "decays.dat")
    decays = clone_decay_map(decays_pdg2020)
    established_abs_ids = parse_mcd_established_abs_ids(mcd_path)
    mcd_mass_width_by_absid = parse_mcd_mass_width_by_abs_id(mcd_path)
    if not established_abs_ids:
        raise RuntimeError(f"Failed to load established-state IDs from {mcd_path}")
    pdglive_established_nodes = parse_established_pdg_nodes(pdglive_nodes_path)
    if not pdglive_established_nodes:
        raise RuntimeError(f"Failed to load established pdgLive nodes from {pdglive_nodes_path}")

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    source_text = read_pdg_source_text(cur)
    parent_map = fetch_parent_map(cur)
    parent_map_by_name = fetch_part_parent_map_by_name(cur)
    part_name_keys = fetch_part_name_keys(cur)
    established_node_part_keys = fetch_node_part_name_keys(cur, sorted(pdglive_established_nodes))
    established_node_abs_mcids = fetch_node_abs_mcids(cur, sorted(pdglive_established_nodes))
    manual_established_abs_ids = {
        abs(int(spec["mcid"]))
        for spec in MANUAL_ESTABLISHED_STATE_OVERRIDES
        if str(spec["parent_pdgid"]) in pdglive_established_nodes
    }
    established_abs_ids_augmented = set(established_abs_ids) | established_node_abs_mcids | manual_established_abs_ids
    db_abs_ids = set(parent_map.keys())
    allowed_abs_ids = established_abs_ids | db_abs_ids

    baseline_total = len(entries)
    kept_by_name = 0
    filtered_entries: List[ParticleEntry] = []
    for entry in entries:
        if abs(entry.pdgid) in allowed_abs_ids:
            filtered_entries.append(entry)
            continue
        keys = name_match_keys(entry.name)
        if any(key in part_name_keys for key in keys):
            filtered_entries.append(entry)
            kept_by_name += 1
    entries = filtered_entries
    baseline_removed_unestablished = baseline_total - len(entries)
    ids_order = [entry.pdgid for entry in entries]
    entry_ids_set = set(ids_order)
    decays = {pid: channels for pid, channels in decays.items() if pid in entry_ids_set}

    mass_updates = 0
    width_updates = 0
    decay_updates = 0
    added_states = 0
    added_states_with_template = 0
    added_states_with_pdg_modes = 0
    added_states_with_generic_fallback = 0

    # Update existing masses and widths.
    for entry in entries:
        parent = parent_map.get(abs(entry.pdgid))
        if not parent:
            parent_candidates = {
                parent_map_by_name[key]
                for key in name_match_keys(entry.name)
                if key in parent_map_by_name
            }
            if len(parent_candidates) == 1:
                parent = next(iter(parent_candidates))
        updated_mass = None
        if mcd_mass_width_by_absid is not None:
            mcd_mass = mcd_mass_width_by_absid.get(abs(entry.pdgid), (None, None))[0]
            if mcd_mass is not None and mcd_mass > 0.0:
                updated_mass = mcd_mass
        if (updated_mass is None or updated_mass <= 0.0) and parent:
            updated_mass = pick_mass_or_width(cur, parent, "M")
        if updated_mass is not None and updated_mass > 0.0:
            entry.mass = updated_mass
            mass_updates += 1

        updated_width = None
        if mcd_mass_width_by_absid is not None:
            mcd_width = mcd_mass_width_by_absid.get(abs(entry.pdgid), (None, None))[1]
            if mcd_width is not None and mcd_width >= 0.0:
                updated_width = mcd_width
        if (updated_width is None or updated_width < 0.0) and parent:
            updated_width = pick_mass_or_width(cur, parent, "G")
        if updated_width is not None and updated_width >= 0.0:
            if parent:
                width_ul = pick_mass_or_width_upper_limit(cur, parent, "G")
                if width_ul is not None and updated_width > width_ul:
                    updated_width = width_ul
            entry.width = updated_width
            width_updates += 1
        if parent:
            width_ul = pick_mass_or_width_upper_limit(cur, parent, "G")
            if width_ul is not None and entry.width > width_ul:
                entry.width = width_ul
                width_updates += 1
            if width_ul is not None:
                ref_entry = pdg2025_ref_with_charm_by_id.get(entry.pdgid)
                if ref_entry is not None and 0.0 < ref_entry.width <= width_ul and entry.width > ref_entry.width:
                    entry.width = ref_entry.width
                    width_updates += 1

    # Update existing decays where explicit BFX channels map to known daughters.
    known_ids_existing = [entry.pdgid for entry in entries]
    by_absid_existing = {abs(entry.pdgid): entry for entry in entries}
    entries_by_id_existing = {entry.pdgid: entry for entry in entries}
    for parent_id, channels in list(decays.items()):
        parent = parent_map.get(abs(parent_id))
        if not parent:
            continue
        if parent_id not in entries_by_id_existing:
            continue
        parent_entry = entries_by_id_existing[parent_id]
        pdg_channels_raw = fetch_explicit_bfx_modes(cur, parent)
        pdg_channels = channels_with_known_daughters(pdg_channels_raw, known_ids_existing)
        pdg_channels = filter_channels_for_parent(parent_entry, pdg_channels, by_absid_existing)
        if not pdg_channels:
            continue
        new_channels, changed = update_channels_from_pdg(channels, pdg_channels)
        if changed:
            decays[parent_id] = new_channels
            decay_updates += 1

    # Expand by additional light/strange/charmed states from PDG2025.
    expansion_rows = fetch_expansion_rows(cur, ids_order)
    new_entries: List[ParticleEntry] = []
    new_parent_by_id: Dict[int, str] = {}
    existing_abs_ids_all = {abs(pid) for pid in ids_order}
    for row in expansion_rows:
        entry = build_entry_from_row(cur, row, mcd_mass_width_by_absid)
        if entry is None:
            continue
        if abs(entry.pdgid) in EXCLUDED_EXPANSION_MCIDS:
            continue
        new_entries.append(entry)
        new_parent_by_id[entry.pdgid] = str(row["pdgid"])
        existing_abs_ids_all.add(abs(entry.pdgid))

    # Add selected established pdgViewer states that are missing MCIDs in DB rows.
    for spec in MANUAL_ESTABLISHED_STATE_OVERRIDES:
        parent_pdgid = str(spec["parent_pdgid"])
        mcid = int(spec["mcid"])
        if parent_pdgid not in pdglive_established_nodes:
            continue
        if abs(mcid) in existing_abs_ids_all:
            continue
        entry = build_entry_from_manual_override(
            cur=cur,
            mcid=mcid,
            parent_pdgid=parent_pdgid,
            name=str(spec["name"]),
            flags=str(spec["flags"]),
            quantum_j=str(spec["quantum_j"]) if spec["quantum_j"] is not None else None,
            charge=float(spec["charge"]) if spec["charge"] is not None else None,
        )
        if entry is None:
            continue
        ref_entry = pdg2025_ref_with_charm_by_id.get(entry.pdgid)
        if ref_entry is not None and is_charm_sector_entry(entry):
            if ref_entry.mass > 0.0:
                entry.mass = ref_entry.mass
        new_entries.append(entry)
        new_parent_by_id[entry.pdgid] = parent_pdgid
        existing_abs_ids_all.add(abs(entry.pdgid))

    # Deterministic order for new states.
    new_entries.sort(key=lambda item: (item.mass, item.pdgid))
    new_entry_ids = {entry.pdgid for entry in new_entries}
    new_charm_ref_threshold_ids = {
        entry.pdgid
        for entry in new_entries
        if is_charm_sector_entry(entry) and entry.pdgid in pdg2025_ref_with_charm_by_id
    }
    recompute_threshold_ids = new_entry_ids - new_charm_ref_threshold_ids
    entries.extend(new_entries)
    ids_order.extend([entry.pdgid for entry in new_entries])
    added_states = len(new_entries)

    entries_by_id = {entry.pdgid: entry for entry in entries}
    mass_by_absid = {abs(entry.pdgid): entry.mass for entry in entries}
    by_absid_full = {abs(entry.pdgid): entry for entry in entries}
    known_ids_full = [entry.pdgid for entry in entries]

    # Build decays for newly added states.
    for entry in new_entries:
        parent = new_parent_by_id.get(entry.pdgid)
        pdg_channels_raw = fetch_explicit_bfx_modes(cur, parent) if parent else []
        pdg_channels = channels_with_known_daughters(pdg_channels_raw, known_ids_full)
        pdg_channels = filter_channels_for_parent(entry, pdg_channels, by_absid_full)
        template = select_template_decay(entry, entries_by_id, decays)

        merged_channels = merge_pdg_and_template(pdg_channels, template or [])
        if pdg_channels:
            added_states_with_pdg_modes += 1
        if template:
            added_states_with_template += 1

        if not merged_channels:
            merged_channels = fallback_generic_decay(entry, entries_by_id)
            if merged_channels:
                added_states_with_generic_fallback += 1

        if merged_channels:
            normalize_channels(merged_channels)
            decays[entry.pdgid] = merged_channels
            entry.threshold = compute_hadronic_strong_threshold(entry, merged_channels, by_absid_full, mass_by_absid)
        else:
            # If no defensible decay model could be built, keep the state but mark stable.
            entry.stable = 1
            entry.width = 0.0
            entry.threshold = 0.0

        # For manual charm additions with synthetic IDs, keep PDG2025 list
        # conventions where available.
        ref_entry = pdg2025_ref_with_charm_by_id.get(entry.pdgid)
        if ref_entry is not None and is_charm_sector_entry(entry):
            if entry.width <= 0.0 and ref_entry.width > 0.0:
                entry.width = ref_entry.width
            entry.threshold = ref_entry.threshold

    # For selected light-sector states, keep PDG2025 reference decay blocks.
    # This avoids template-induced mode drift for sparse/qualitative API channels.
    for pid in sorted(LIGHT_DECAY_REFERENCE_OVERRIDE_IDS):
        parent_entry = entries_by_id.get(pid)
        ref_channels = pdg2025_ref_decays.get(pid)
        if parent_entry is None or not ref_channels:
            continue
        override_channels = channels_with_known_daughters(ref_channels, known_ids_full)
        override_channels = filter_channels_for_parent(parent_entry, override_channels, by_absid_full)
        if not override_channels or not normalize_channels(override_channels):
            continue
        decays[pid] = override_channels
        if parent_entry.stable == 0:
            parent_entry.threshold = compute_hadronic_strong_threshold(
                parent_entry,
                override_channels,
                by_absid_full,
                mass_by_absid,
            )

    # Explicit manual decay overrides (applied after PDG reference overrides).
    for pid, override_spec in sorted(MANUAL_DECAY_OVERRIDE_CHANNELS.items()):
        parent_entry = entries_by_id.get(pid)
        if parent_entry is None:
            continue
        override_channels = [
            DecayChannel(br=channel.br, daughters=list(channel.daughters))
            for channel in override_spec
        ]
        override_channels = channels_with_known_daughters(override_channels, known_ids_full)
        override_channels = filter_channels_for_parent(parent_entry, override_channels, by_absid_full)
        if not override_channels or not normalize_channels(override_channels):
            continue
        decays[pid] = override_channels
        if parent_entry.stable == 0:
            parent_entry.threshold = compute_hadronic_strong_threshold(
                parent_entry,
                override_channels,
                by_absid_full,
                mass_by_absid,
            )

    conn.close()

    # Build extra (non-default) output set.
    entries_full = clone_particle_entries(entries)
    entries_full, ids_order_full = sorted_entries_and_ids(entries_full, pdg2020_rank_map)
    decays_full = sanitize_decays_for_entries(
        entries_full,
        ids_order_full,
        decays,
        recompute_threshold_ids=recompute_threshold_ids,
    )
    names_full = {entry.pdgid: entry.name for entry in entries_full}

    # Write extra output.
    write_particle_list(
        out_list_extra,
        entries_full,
        source_text,
        description_line="PDG2025 extra hadron list containing established and non-established hadrons (incl. charm) from DB",
    )

    # Build established output set with charm from the extra list.
    entries_established_with_charm: List[ParticleEntry] = []
    for entry in entries_full:
        if abs(entry.pdgid) in established_abs_ids_augmented:
            entries_established_with_charm.append(entry)
            continue
        entry_keys = name_match_keys(entry.name)
        if any(key in established_node_part_keys for key in entry_keys):
            entries_established_with_charm.append(entry)
    entries_established_with_charm = [
        entry
        for entry in entries_established_with_charm
        if abs(entry.pdgid) not in EXCLUDED_ESTABLISHED_WITH_CHARM_PDGS
    ]
    entries_established_with_charm = clone_particle_entries(entries_established_with_charm)
    established_with_charm_id_set = {entry.pdgid for entry in entries_established_with_charm}
    ids_order_established_with_charm = [pid for pid in ids_order_full if pid in established_with_charm_id_set]

    write_particle_list(
        out_list_with_charm,
        entries_established_with_charm,
        source_text,
        description_line="PDG2025 list containing established hadrons including charm",
    )

    # Build established default output without any charm content (open or hidden).
    entries_established = clone_particle_entries(
        [entry for entry in entries_established_with_charm if not is_charm_sector_entry(entry)]
    )
    established_id_set = {entry.pdgid for entry in entries_established}
    ids_order_established = [pid for pid in ids_order_established_with_charm if pid in established_id_set]

    write_particle_list(
        out_list,
        entries_established,
        source_text,
        description_line="PDG2025 list containing established strange and non-strange hadrons only",
    )

    entries_with_nuclei = clone_particle_entries(entries_established)
    ids_order_with_nuclei = list(ids_order_established)
    ids_with_nuclei_set = set(ids_order_with_nuclei)
    for entry in nuclei_entries:
        if entry.pdgid in ids_with_nuclei_set:
            continue
        entries_with_nuclei.append(
            ParticleEntry(
                pdgid=entry.pdgid,
                name=entry.name,
                stable=entry.stable,
                mass=entry.mass,
                degeneracy=entry.degeneracy,
                statistics=entry.statistics,
                baryon=entry.baryon,
                charge=entry.charge,
                strange=entry.strange,
                charm=entry.charm,
                abs_s=entry.abs_s,
                abs_c=entry.abs_c,
                width=entry.width,
                threshold=entry.threshold,
            )
        )
        ids_order_with_nuclei.append(entry.pdgid)
        ids_with_nuclei_set.add(entry.pdgid)

    entries_with_nuclei, ids_order_with_nuclei = sorted_entries_and_ids(entries_with_nuclei, pdg2020_rank_map)

    write_particle_list(
        out_list_with_nuclei,
        entries_with_nuclei,
        source_text,
        description_line="PDG2025 list containing established hadrons and stable light nuclei",
    )

    entries_with_excited_nuclei = clone_particle_entries(entries_with_nuclei)
    ids_order_with_excited_nuclei = list(ids_order_with_nuclei)
    ids_with_excited_nuclei_set = set(ids_order_with_excited_nuclei)
    for entry in excited_nuclei_entries:
        if entry.pdgid in ids_with_excited_nuclei_set:
            continue
        entries_with_excited_nuclei.append(
            ParticleEntry(
                pdgid=entry.pdgid,
                name=entry.name,
                stable=entry.stable,
                mass=entry.mass,
                degeneracy=entry.degeneracy,
                statistics=entry.statistics,
                baryon=entry.baryon,
                charge=entry.charge,
                strange=entry.strange,
                charm=entry.charm,
                abs_s=entry.abs_s,
                abs_c=entry.abs_c,
                width=entry.width,
                threshold=entry.threshold,
            )
        )
        ids_order_with_excited_nuclei.append(entry.pdgid)
        ids_with_excited_nuclei_set.add(entry.pdgid)

    for entry in entries_with_excited_nuclei:
        threshold = pdg2020_nuclei_thresholds.get(entry.pdgid)
        if threshold is not None:
            entry.threshold = threshold

    entries_with_excited_nuclei, ids_order_with_excited_nuclei = sorted_entries_and_ids(
        entries_with_excited_nuclei,
        pdg2020_rank_map,
    )

    write_particle_list(
        out_list_with_excited_nuclei,
        entries_with_excited_nuclei,
        source_text,
        description_line="PDG2025 list containing established hadrons, stable and excited light nuclei",
    )

    # Build a single canonical decays.dat with all generated states:
    # full hadron set + stable nuclei + excited nuclei.
    entries_all = clone_particle_entries(entries_full)
    ids_order_all = list(ids_order_full)
    ids_all_set = set(ids_order_all)
    for entry in nuclei_entries + excited_nuclei_entries:
        if entry.pdgid in ids_all_set:
            continue
        entries_all.append(
            ParticleEntry(
                pdgid=entry.pdgid,
                name=entry.name,
                stable=entry.stable,
                mass=entry.mass,
                degeneracy=entry.degeneracy,
                statistics=entry.statistics,
                baryon=entry.baryon,
                charge=entry.charge,
                strange=entry.strange,
                charm=entry.charm,
                abs_s=entry.abs_s,
                abs_c=entry.abs_c,
                width=entry.width,
                threshold=entry.threshold,
            )
        )
        ids_order_all.append(entry.pdgid)
        ids_all_set.add(entry.pdgid)

    entries_all, ids_order_all = sorted_entries_and_ids(entries_all, pdg2020_rank_map)

    decays_all_seed = clone_decay_map(decays_full)
    for pid in ids_order_all:
        if pid not in decays_all_seed and pid in decays_pdg2020:
            decays_all_seed[pid] = [
                DecayChannel(br=channel.br, daughters=list(channel.daughters))
                for channel in decays_pdg2020[pid]
            ]
    decays_all = sanitize_decays_for_entries(
        entries_all,
        ids_order_all,
        decays_all_seed,
        recompute_threshold_ids=recompute_threshold_ids,
    )
    names_all = {entry.pdgid: entry.name for entry in entries_all}
    write_decay_file(out_decays, decays_all, ids_order_all, names_all, source_text)

    print(f"Wrote {out_list}")
    print(f"Wrote {out_decays}")
    print(f"Wrote {out_list_with_charm}")
    print(f"Wrote {out_list_extra}")
    print(f"Wrote {out_list_with_nuclei}")
    print(f"Wrote {out_list_with_excited_nuclei}")
    print(f"Removed baseline states absent in both DB and mcd: {baseline_removed_unestablished}")
    print(f"Kept baseline states by DB PART name match: {kept_by_name}")
    print(f"Updated masses: {mass_updates}")
    print(f"Updated widths: {width_updates}")
    print(f"Updated existing decay blocks: {decay_updates}")
    print(f"Added states: {added_states}")
    print(f"Extra-list states: {len(entries_full)}")
    print(f"Established-with-charm states: {len(entries_established_with_charm)}")
    print(f"Established-only states: {len(entries_established)}")
    print(f"With stable nuclei states: {len(entries_with_nuclei)}")
    print(f"With excited nuclei states: {len(entries_with_excited_nuclei)}")
    print(f"Established nodes from pdgLive tables: {len(pdglive_established_nodes)}")
    print(f"Added states with explicit PDG modes: {added_states_with_pdg_modes}")
    print(f"Added states with template decays: {added_states_with_template}")
    print(f"Added states with generic fallback decays: {added_states_with_generic_fallback}")


if __name__ == "__main__":
    main()
