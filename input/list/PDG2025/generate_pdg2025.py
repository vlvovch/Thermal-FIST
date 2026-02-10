#!/usr/bin/env python3
"""
Generate PDG2025 particle list and decay list for Thermal-FIST
Based on PDG2020 as template, updated with PDG 2025 mass-width data
"""

import sys
import os
import re
from collections import OrderedDict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)                    # input/list
PDG2020_LIST = os.path.join(BASE_DIR, "PDG2020", "list.dat")
PDG2020_LIST_ALL = os.path.join(BASE_DIR, "PDG2020", "list-all.dat")
PDG2020_DECAYS = os.path.join(BASE_DIR, "PDG2020", "decays.dat")
PDG2025_MASSWIDTH = os.path.join(SCRIPT_DIR, "pdglisting", "mass_width_2025.txt")
OUTPUT_DIR = SCRIPT_DIR
OUTPUT_LIST = os.path.join(OUTPUT_DIR, "list.dat")
OUTPUT_DECAYS = os.path.join(OUTPUT_DIR, "decays.dat")

# =====================================================================
# PDG 2025 mass-width data (parsed from mass_width_2025.txt)
# Only include light-flavor hadrons relevant for Thermal-FIST
# =====================================================================

def parse_pdg2025_masswidth(filename):
    """Parse the PDG 2025 mass-width file"""
    entries = {}
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('*'):
                continue
            line = line.rstrip('\n')
            if len(line) < 107:
                continue
            # Parse PDG IDs (up to 4, 8 chars each)
            ids = []
            for i in range(4):
                s = line[i*8:(i+1)*8].strip()
                if s:
                    ids.append(int(s))
            # Parse mass
            mass_str = line[33:51].strip()
            mass = float(mass_str) if mass_str else None
            # Parse width
            width_str = line[70:88].strip()
            width = float(width_str) if width_str else 0.0
            # Parse name
            name = line[107:].strip() if len(line) > 107 else ""

            for pid in ids:
                entries[pid] = {
                    'mass': mass,
                    'width': width,
                    'name': name,
                    'all_ids': ids
                }
    return entries

def parse_pdg2020_list(filename):
    """Parse the PDG 2020 particle list"""
    particles = []
    with open(filename, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 14:
                continue
            p = {
                'pdgid': int(parts[0]),
                'name': parts[1],
                'stable': int(parts[2]),
                'mass': float(parts[3]),
                'degeneracy': int(parts[4]),
                'statistics': int(parts[5]),
                'B': int(parts[6]),
                'Q': int(parts[7]),
                'S': int(parts[8]),
                'C': int(parts[9]),
                'absS': float(parts[10]),
                'absC': float(parts[11]),
                'width': float(parts[12]),
                'threshold': float(parts[13]),
            }
            particles.append(p)
    return particles

def parse_pdg2020_decays(filename):
    """Parse the PDG 2020 decays file"""
    decays = OrderedDict()
    with open(filename, 'r', encoding='utf-8', errors='replace') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        # Skip empty lines and comments
        if not line or line.startswith('#'):
            i += 1
            continue

        # Parse PDG ID
        parts = line.split('#')[0].strip().split()
        if not parts:
            i += 1
            continue

        try:
            pdgid = int(parts[0])
        except ValueError:
            i += 1
            continue

        i += 1
        # Parse number of decay channels
        if i >= len(lines):
            break
        nch_line = lines[i].strip().split('#')[0].strip().split()
        if not nch_line:
            continue
        nch = int(nch_line[0])
        i += 1

        channels = []
        for j in range(nch):
            if i >= len(lines):
                break
            ch_line = lines[i].strip()
            comment = ""
            if '#' in ch_line:
                comment = ch_line[ch_line.index('#'):]
                ch_line = ch_line[:ch_line.index('#')].strip()
            ch_parts = ch_line.split()
            if len(ch_parts) >= 2:
                br = float(ch_parts[0])
                daughters = [int(x) for x in ch_parts[1:]]
                channels.append({
                    'br': br,
                    'daughters': daughters,
                    'comment': comment
                })
            i += 1

        decays[pdgid] = channels

    return decays

# =====================================================================
# PDG IDs that should be REMOVED (not in PDG2025 established list)
# These are states from PDG2020 that are not in the PDG2025 mass-width file
# =====================================================================

# States to remove: those not established (*** or ****) in PDG2025
# Based on cross-checking PDG 2025 Summary Tables and Review articles
# Note: The mass_width_2025.txt file is INCOMPLETE -- it only has particles
# with well-determined numerical mass/width values. Many established states
# (especially baryons with mass/width given as ranges) are missing from it.
# We cross-checked all states against the actual PDG 2025 Summary Tables.
REMOVE_IDS = set([
    # pi(1)(1400) - OMITTED from PDG2025 Summary Table
    # PDG says coupled-channel analyses favor only one broad 1-+ state: pi(1)(1600)
    9000113, 9000213,
    # Sigma(2250) - only ** (two-star) status in PDG2025, not established
    33114, 33214, 33224,
])

# =====================================================================
# Additional mass/width updates for particles NOT in mass_width_2025.txt
# but still established in PDG 2025 (verified from Summary Tables/listings)
# =====================================================================
ADDITIONAL_UPDATES = {
    # phi(2170): mass 2175->2164 MeV, width 61->88 MeV (significant update, many new BESIII measurements)
    200333: {'mass': 2.164, 'width': 0.088},
    # pi2(1880): width 235->237 MeV (minor)
    20115: {'width': 0.237},
    20215: {'width': 0.237},
    # Omega(2012): width 6.11->6.4 MeV (minor)
    9903334: {'width': 0.0064},
}

# =====================================================================
# NEW states to add from PDG2025
# =====================================================================

# New light unflavored mesons
NEW_MESONS = [
    # f(2)(1565): I^G(J^PC) = 0+(2++), PDG ID 9010225
    {
        'pdgid': 9010225,
        'name': 'f(2)(1565)',
        'stable': 0,
        'mass': 1.571,
        'degeneracy': 5,  # 2J+1 = 5 for J=2
        'statistics': -1,  # boson
        'B': 0, 'Q': 0, 'S': 0, 'C': 0,
        'absS': 0, 'absC': 0,
        'width': 0.133,
        'threshold': 0  # will be computed from decays
    },
    # f(0)(2020): I^G(J^PC) = 0+(0++), PDG ID 9050221
    {
        'pdgid': 9050221,
        'name': 'f(0)(2020)',
        'stable': 0,
        'mass': 1.982,
        'degeneracy': 1,  # 2J+1 = 1 for J=0
        'statistics': -1,
        'B': 0, 'Q': 0, 'S': 0, 'C': 0,
        'absS': 0, 'absC': 0,
        'width': 0.44,
        'threshold': 0
    },
    # f(2)(2150): I^G(J^PC) = 0+(2++), PDG ID 9070225
    {
        'pdgid': 9070225,
        'name': 'f(2)(2150)',
        'stable': 0,
        'mass': 2.157,
        'degeneracy': 5,
        'statistics': -1,
        'B': 0, 'Q': 0, 'S': 0, 'C': 0,
        'absS': 0, 'absC': 0,
        'width': 0.152,
        'threshold': 0
    },
]

# New strange mesons
NEW_STRANGE_MESONS = [
    # K(0)*(700) / kappa: I(J^P) = 1/2(0+), PDG IDs 9000311, 9000321
    # degeneracy=0: cancellation in I=3/2 Kpi repulsive channel, analogous to f(0)(500)/sigma in I=2 pipi
    {
        'pdgid': 9000311,
        'name': 'K(0)*(700)0',
        'stable': 0,
        'mass': 0.838,
        'degeneracy': 0,
        'statistics': -1,
        'B': 0, 'Q': 0, 'S': 1, 'C': 0,
        'absS': 1, 'absC': 0,
        'width': 0.463,
        'threshold': 0
    },
    {
        'pdgid': 9000321,
        'name': 'K(0)*(700)+',
        'stable': 0,
        'mass': 0.838,
        'degeneracy': 0,
        'statistics': -1,
        'B': 0, 'Q': 1, 'S': 1, 'C': 0,
        'absS': 1, 'absC': 0,
        'width': 0.463,
        'threshold': 0
    },
    # K(1)(1650): I(J^P) = 1/2(1+), PDG IDs 9000313, 9000323
    {
        'pdgid': 9000313,
        'name': 'K(1)(1650)0',
        'stable': 0,
        'mass': 1.65,
        'degeneracy': 3,
        'statistics': -1,
        'B': 0, 'Q': 0, 'S': 1, 'C': 0,
        'absS': 1, 'absC': 0,
        'width': 0.15,
        'threshold': 0
    },
    {
        'pdgid': 9000323,
        'name': 'K(1)(1650)+',
        'stable': 0,
        'mass': 1.65,
        'degeneracy': 3,
        'statistics': -1,
        'B': 0, 'Q': 1, 'S': 1, 'C': 0,
        'absS': 1, 'absC': 0,
        'width': 0.15,
        'threshold': 0
    },
    # K(0)*(1950): I(J^P) = 1/2(0+), PDG IDs 9020311, 9020321
    {
        'pdgid': 9020311,
        'name': 'K(0)*(1950)0',
        'stable': 0,
        'mass': 1.957,
        'degeneracy': 1,
        'statistics': -1,
        'B': 0, 'Q': 0, 'S': 1, 'C': 0,
        'absS': 1, 'absC': 0,
        'width': 0.17,
        'threshold': 0
    },
    {
        'pdgid': 9020321,
        'name': 'K(0)*(1950)+',
        'stable': 0,
        'mass': 1.957,
        'degeneracy': 1,
        'statistics': -1,
        'B': 0, 'Q': 1, 'S': 1, 'C': 0,
        'absS': 1, 'absC': 0,
        'width': 0.17,
        'threshold': 0
    },
    # K(2)*(1980): I(J^P) = 1/2(2+), PDG IDs 9010315, 9010325
    {
        'pdgid': 9010315,
        'name': 'K(2)*(1980)0',
        'stable': 0,
        'mass': 1.99,
        'degeneracy': 5,
        'statistics': -1,
        'B': 0, 'Q': 0, 'S': 1, 'C': 0,
        'absS': 1, 'absC': 0,
        'width': 0.348,
        'threshold': 0
    },
    {
        'pdgid': 9010325,
        'name': 'K(2)*(1980)+',
        'stable': 0,
        'mass': 1.99,
        'degeneracy': 5,
        'statistics': -1,
        'B': 0, 'Q': 1, 'S': 1, 'C': 0,
        'absS': 1, 'absC': 0,
        'width': 0.348,
        'threshold': 0
    },
]

# =====================================================================
# NEW decays for new states
# =====================================================================

NEW_DECAYS = {
    # f(2)(1565): I^G(J^PC) = 0+(2++)
    # Similar to f(2)(1270) but heavier. Decays to pi pi, KK, eta eta
    9010225: [
        {'br': 0.50, 'daughters': [211, -211], 'comment': '# f(2)(1565) -> pi+ + pi-'},
        {'br': 0.25, 'daughters': [111, 111], 'comment': '# f(2)(1565) -> pi0 + pi0'},
        {'br': 0.08, 'daughters': [321, -321], 'comment': '# f(2)(1565) -> K+ + K-'},
        {'br': 0.08, 'daughters': [311, -311], 'comment': '# f(2)(1565) -> K0 + anti-K0'},
        {'br': 0.09, 'daughters': [221, 221], 'comment': '# f(2)(1565) -> eta + eta'},
    ],

    # f(0)(2020): I^G(J^PC) = 0+(0++)
    # PDG 2025 seen modes: rho pi pi, pi0 pi0, rho rho, omega omega, eta eta, eta' eta'
    # Note: KKbar is NOT listed in PDG 2025. Replaced with eta'eta'.
    # rho pi pi / rho rho -> 4pi final states
    9050221: [
        {'br': 0.20, 'daughters': [211, -211], 'comment': '# f(0)(2020) -> pi+ + pi-'},
        {'br': 0.10, 'daughters': [111, 111], 'comment': '# f(0)(2020) -> pi0 + pi0'},
        {'br': 0.10, 'daughters': [221, 221], 'comment': '# f(0)(2020) -> eta + eta'},
        {'br': 0.10, 'daughters': [331, 331], 'comment': "# f(0)(2020) -> eta'(958) + eta'(958)"},
        {'br': 0.25, 'daughters': [211, -211, 211, -211], 'comment': '# f(0)(2020) -> pi+ + pi- + pi+ + pi-'},
        {'br': 0.25, 'daughters': [211, -211, 111, 111], 'comment': '# f(0)(2020) -> pi+ + pi- + pi0 + pi0'},
    ],

    # f(2)(2150): I^G(J^PC) = 0+(2++)
    # PDG 2025: seen modes = pipi, phiphi, etaeta, KKbar, f2(1270)eta, a2(1320)pi, ppbar
    # phiphi seen in PDG 2025, added; a2(1320)pi also seen
    9070225: [
        {'br': 0.15, 'daughters': [211, -211], 'comment': '# f(2)(2150) -> pi+ + pi-'},
        {'br': 0.075, 'daughters': [111, 111], 'comment': '# f(2)(2150) -> pi0 + pi0'},
        {'br': 0.10, 'daughters': [321, -321], 'comment': '# f(2)(2150) -> K+ + K-'},
        {'br': 0.10, 'daughters': [311, -311], 'comment': '# f(2)(2150) -> K0 + anti-K0'},
        {'br': 0.10, 'daughters': [221, 221], 'comment': '# f(2)(2150) -> eta + eta'},
        {'br': 0.10, 'daughters': [333, 333], 'comment': '# f(2)(2150) -> phi(1020) + phi(1020)'},
        {'br': 0.10, 'daughters': [225, 221], 'comment': '# f(2)(2150) -> f(2)(1270) + eta'},
        {'br': 0.10, 'daughters': [323, -323], 'comment': '# f(2)(2150) -> K*(892)+ + K*(892)-'},
        {'br': 0.075, 'daughters': [115, 111], 'comment': '# f(2)(2150) -> a(2)(1320)0 + pi0'},
        {'br': 0.05, 'daughters': [215, -211], 'comment': '# f(2)(2150) -> a(2)(1320)+ + pi-'},
        {'br': 0.05, 'daughters': [-215, 211], 'comment': '# f(2)(2150) -> a(2)(1320)- + pi+'},
    ],

    # K(0)*(700)0: I(J^P) = 1/2(0+), kappa
    # Dominant decay to K pi
    9000311: [
        {'br': 0.6667, 'daughters': [321, -211], 'comment': '# K(0)*(700)0 -> K+ + pi-'},
        {'br': 0.3333, 'daughters': [311, 111], 'comment': '# K(0)*(700)0 -> K0 + pi0'},
    ],

    # K(0)*(700)+
    9000321: [
        {'br': 0.6667, 'daughters': [311, 211], 'comment': '# K(0)*(700)+ -> K0 + pi+'},
        {'br': 0.3333, 'daughters': [321, 111], 'comment': '# K(0)*(700)+ -> K+ + pi0'},
    ],

    # K(1)(1650)0: I(J^P) = 1/2(1+)
    # Axial kaon, similar to K(1)(1270) but heavier
    # Decays to K*(892)pi, K*rho, K omega, K phi
    9000313: [
        {'br': 0.20, 'daughters': [313, 111], 'comment': '# K(1)(1650)0 -> K*(892)0 + pi0'},
        {'br': 0.40, 'daughters': [323, -211], 'comment': '# K(1)(1650)0 -> K*(892)+ + pi-'},
        {'br': 0.10, 'daughters': [311, 113], 'comment': '# K(1)(1650)0 -> K0 + rho(770)0'},
        {'br': 0.20, 'daughters': [321, -213], 'comment': '# K(1)(1650)0 -> K+ + rho(770)-'},
        {'br': 0.10, 'daughters': [311, 223], 'comment': '# K(1)(1650)0 -> K0 + omega(782)'},
    ],

    # K(1)(1650)+
    9000323: [
        {'br': 0.40, 'daughters': [313, 211], 'comment': '# K(1)(1650)+ -> K*(892)0 + pi+'},
        {'br': 0.20, 'daughters': [323, 111], 'comment': '# K(1)(1650)+ -> K*(892)+ + pi0'},
        {'br': 0.20, 'daughters': [311, 213], 'comment': '# K(1)(1650)+ -> K0 + rho(770)+'},
        {'br': 0.10, 'daughters': [321, 113], 'comment': '# K(1)(1650)+ -> K+ + rho(770)0'},
        {'br': 0.10, 'daughters': [321, 223], 'comment': '# K(1)(1650)+ -> K+ + omega(782)'},
    ],

    # K(0)*(1950)0: I(J^P) = 1/2(0+)
    # Scalar kaon, decays to K pi, K eta
    9020311: [
        {'br': 0.40, 'daughters': [321, -211], 'comment': '# K(0)*(1950)0 -> K+ + pi-'},
        {'br': 0.20, 'daughters': [311, 111], 'comment': '# K(0)*(1950)0 -> K0 + pi0'},
        {'br': 0.20, 'daughters': [311, 221], 'comment': '# K(0)*(1950)0 -> K0 + eta'},
        {'br': 0.10, 'daughters': [321, -321, 311], 'comment': '# K(0)*(1950)0 -> K+ + K- + K0'},
        {'br': 0.10, 'daughters': [311, 331], 'comment': "# K(0)*(1950)0 -> K0 + eta'(958)"},
    ],

    # K(0)*(1950)+
    9020321: [
        {'br': 0.40, 'daughters': [311, 211], 'comment': '# K(0)*(1950)+ -> K0 + pi+'},
        {'br': 0.20, 'daughters': [321, 111], 'comment': '# K(0)*(1950)+ -> K+ + pi0'},
        {'br': 0.20, 'daughters': [321, 221], 'comment': '# K(0)*(1950)+ -> K+ + eta'},
        {'br': 0.10, 'daughters': [321, -311, 311], 'comment': '# K(0)*(1950)+ -> K+ + anti-K0 + K0'},
        {'br': 0.10, 'daughters': [321, 331], 'comment': "# K(0)*(1950)+ -> K+ + eta'(958)"},
    ],

    # K(2)*(1980)0: I(J^P) = 1/2(2+)
    # PDG 2025 modes: K*(892)pi (possibly seen), K rho (possibly seen), K f2(1270) (possibly seen),
    #                 K phi (seen), K eta (seen)
    # Replaced K omega (not in PDG) with K phi and K eta (both "seen" in PDG 2025)
    9010315: [
        {'br': 0.10, 'daughters': [321, -211], 'comment': '# K(2)*(1980)0 -> K+ + pi-'},
        {'br': 0.05, 'daughters': [311, 111], 'comment': '# K(2)*(1980)0 -> K0 + pi0'},
        {'br': 0.10, 'daughters': [313, 111], 'comment': '# K(2)*(1980)0 -> K*(892)0 + pi0'},
        {'br': 0.20, 'daughters': [323, -211], 'comment': '# K(2)*(1980)0 -> K*(892)+ + pi-'},
        {'br': 0.075, 'daughters': [311, 113], 'comment': '# K(2)*(1980)0 -> K0 + rho(770)0'},
        {'br': 0.15, 'daughters': [321, -213], 'comment': '# K(2)*(1980)0 -> K+ + rho(770)-'},
        {'br': 0.10, 'daughters': [311, 225], 'comment': '# K(2)*(1980)0 -> K0 + f(2)(1270)'},
        {'br': 0.075, 'daughters': [311, 333], 'comment': '# K(2)*(1980)0 -> K0 + phi(1020)'},
        {'br': 0.15, 'daughters': [311, 221], 'comment': '# K(2)*(1980)0 -> K0 + eta'},
    ],

    # K(2)*(1980)+
    9010325: [
        {'br': 0.10, 'daughters': [311, 211], 'comment': '# K(2)*(1980)+ -> K0 + pi+'},
        {'br': 0.05, 'daughters': [321, 111], 'comment': '# K(2)*(1980)+ -> K+ + pi0'},
        {'br': 0.20, 'daughters': [313, 211], 'comment': '# K(2)*(1980)+ -> K*(892)0 + pi+'},
        {'br': 0.10, 'daughters': [323, 111], 'comment': '# K(2)*(1980)+ -> K*(892)+ + pi0'},
        {'br': 0.15, 'daughters': [311, 213], 'comment': '# K(2)*(1980)+ -> K0 + rho(770)+'},
        {'br': 0.075, 'daughters': [321, 113], 'comment': '# K(2)*(1980)+ -> K+ + rho(770)0'},
        {'br': 0.10, 'daughters': [321, 225], 'comment': '# K(2)*(1980)+ -> K+ + f(2)(1270)'},
        {'br': 0.075, 'daughters': [321, 333], 'comment': '# K(2)*(1980)+ -> K+ + phi(1020)'},
        {'br': 0.15, 'daughters': [321, 221], 'comment': '# K(2)*(1980)+ -> K+ + eta'},
    ],
}

# =====================================================================
# Decay updates for EXISTING particles (from PDG2020) based on PDG 2025
# These replace the PDG2020 decay channels entirely
# =====================================================================

DECAY_UPDATES = {
    # eta: Updated branching fractions from PDG 2025
    # PDG 2025: 2gamma (39.36%), 3pi0 (32.56%), pi+pi-pi0 (23.02%), pi+pi-gamma (4.28%),
    #           e+e-gamma (0.70%), mu+mu-gamma (0.031%)
    # Changes vs PDG2020: gamma gamma 39.41->39.36, 3pi0 32.68->32.56, pi+pi-pi0 22.92->23.02,
    #                     pi+pi-gamma 4.22->4.28, e+e-gamma 0.73->0.70, mu+mu-gamma 0.04->0.031
    221: [
        {'br': 0.3256, 'daughters': [111, 111, 111], 'comment': '# eta -> pi0 + pi0 + pi0'},
        {'br': 0.2302, 'daughters': [211, -211, 111], 'comment': '# eta -> pi+ + pi- + pi0'},
        {'br': 0.0428, 'daughters': [211, -211, 22], 'comment': '# eta -> pi+ + pi- + gamma'},
        {'br': 0.3936, 'daughters': [22, 22], 'comment': '# eta -> gamma + gamma'},
        {'br': 0.0070, 'daughters': [-11, 11, 22], 'comment': '# eta -> e+ + e- + gamma'},
        {'br': 0.0008, 'daughters': [-13, 13, 22], 'comment': '# eta -> mu+ + mu- + gamma'},
    ],

    # omega(782): Updated branching fractions from PDG 2025
    # PDG 2025 BFX: pi+pi-pi0 (89.2%), pi0 gamma (8.33%), pi+pi- (1.53%),
    #               pi0 e+e- (0.077%), eta gamma (0.045%)
    # Note: PDG2020 hand-input had eta gamma=0.45% and pi0 e+e-=0.49% (10x errors)
    # Remaining 0.82% absorbed into pi+pi-pi0
    223: [
        {'br': 0.9002, 'daughters': [211, -211, 111], 'comment': '# omega(782) -> pi+ + pi- + pi0'},
        {'br': 0.0833, 'daughters': [111, 22], 'comment': '# omega(782) -> pi0 + gamma'},
        {'br': 0.0153, 'daughters': [211, -211], 'comment': '# omega(782) -> pi+ + pi-'},
        {'br': 0.00045, 'daughters': [221, 22], 'comment': '# omega(782) -> eta + gamma'},
        {'br': 0.00077, 'daughters': [111, -11, 11], 'comment': '# omega(782) -> pi0 + e+ + e-'},
    ],

    # eta'(958): Updated branching fractions from PDG 2025
    # PDG 2025 BFX: pi+pi-eta (42.5%), rho0 gamma (29.48%), pi0pi0eta (22.4%),
    #               omega gamma (2.52%), gamma gamma (2.307%), pi+pi-pi0 (0.361%)
    # Fix: pi+pi-pi0 was 0.8% in PDG2020 hand-input, PDG2025 says 0.361%
    # Remaining 0.73% absorbed into dominant pi+pi-eta channel
    331: [
        {'br': 0.4294, 'daughters': [211, -211, 221], 'comment': "# eta'(958) -> pi+ + pi- + eta"},
        {'br': 0.2948, 'daughters': [113, 22], 'comment': "# eta'(958) -> rho(770)0 + gamma"},
        {'br': 0.224, 'daughters': [111, 111, 221], 'comment': "# eta'(958) -> pi0 + pi0 + eta"},
        {'br': 0.0252, 'daughters': [223, 22], 'comment': "# eta'(958) -> omega(782) + gamma"},
        {'br': 0.023, 'daughters': [22, 22], 'comment': "# eta'(958) -> gamma + gamma"},
        {'br': 0.0036, 'daughters': [211, -211, 111], 'comment': "# eta'(958) -> pi+ + pi- + pi0"},
    ],

    # phi(1020): Updated branching fractions from PDG 2025
    # PDG 2025 BFX: K+K- (49.9%), KL KS (33.9%), rho pi + pi+pi-pi0 (15.2%), eta gamma (1.306%)
    # Fix: K+K- was 49.2% in PDG2020, PDG2025 says 49.9%
    # Adjusted KK and 3pi to match PDG2025 values
    333: [
        {'br': 0.499, 'daughters': [321, -321], 'comment': '# phi(1020) -> K+ + K-'},
        {'br': 0.339, 'daughters': [311, -311], 'comment': '# phi(1020) -> K0 + anti-K0'},
        {'br': 0.147, 'daughters': [211, -211, 111], 'comment': '# phi(1020) -> pi+ + pi- + pi0'},
        {'br': 0.013, 'daughters': [221, 22], 'comment': '# phi(1020) -> eta + gamma'},
        {'br': 0.002, 'daughters': [111, 22], 'comment': '# phi(1020) -> pi0 + gamma'},
    ],

    # phi(2170): Broadened decay channels based on PDG 2025 observed modes
    # PDG 2025: phi f0(980) seen, phi eta seen, phi eta' seen, phi pi pi seen,
    #           K+K- pi+pi- seen, KS KL seen, K*(892)+K- + c.c. seen, K1(1270)+K- + c.c. seen
    200333: [
        {'br': 0.20, 'daughters': [333, 9010221], 'comment': '# phi(2170) -> phi(1020) + f(0)(980)'},
        {'br': 0.15, 'daughters': [333, 221], 'comment': '# phi(2170) -> phi(1020) + eta'},
        {'br': 0.10, 'daughters': [333, 331], 'comment': "# phi(2170) -> phi(1020) + eta'(958)"},
        {'br': 0.10, 'daughters': [321, -321, 211, -211], 'comment': '# phi(2170) -> K+ + K- + pi+ + pi-'},
        {'br': 0.10, 'daughters': [321, -321], 'comment': '# phi(2170) -> K+ + K-'},
        {'br': 0.10, 'daughters': [311, -311], 'comment': '# phi(2170) -> K0 + anti-K0'},
        {'br': 0.075, 'daughters': [323, -321], 'comment': '# phi(2170) -> K*(892)+ + K-'},
        {'br': 0.075, 'daughters': [-323, 321], 'comment': '# phi(2170) -> K*(892)- + K+'},
        {'br': 0.05, 'daughters': [313, -311], 'comment': '# phi(2170) -> K*(892)0 + anti-K0'},
        {'br': 0.05, 'daughters': [-313, 311], 'comment': '# phi(2170) -> anti-K*(892)0 + K0'},
    ],

    # =====================================================================
    # Baryon decay updates based on PDG 2025 cross-check
    # =====================================================================

    # N(1675)0: PDG 2025 ranges: Npi (38-42%), Delta(1232)pi (23-37%), Nsigma (3-7%)
    # Current: Npi 42%, Deltapi 54%, Nsigma 4% -> Deltapi too high
    # Updated: Npi 40%, Deltapi 30%, Nsigma 5%, N(1440)pi 25%
    # N(1440)pi is kinematically open (1440+140=1580 < 1675)
    # C-G for I=1/2 N* -> I=1/2 N(1440) + I=1 pi: 1/3, 2/3
    2116: [  # N(1675)0
        {'br': 0.133, 'daughters': [2112, 111], 'comment': '# N(1675)0 -> n + pi0'},
        {'br': 0.267, 'daughters': [2212, -211], 'comment': '# N(1675)0 -> p + pi-'},
        {'br': 0.05, 'daughters': [2214, -211], 'comment': '# N(1675)0 -> Delta(1232)+ + pi-'},
        {'br': 0.10, 'daughters': [2114, 111], 'comment': '# N(1675)0 -> Delta(1232)0 + pi0'},
        {'br': 0.15, 'daughters': [1114, 211], 'comment': '# N(1675)0 -> Delta(1232)- + pi+'},
        {'br': 0.05, 'daughters': [2112, 9000221], 'comment': '# N(1675)0 -> n + f(0)(500)'},
        {'br': 0.083, 'daughters': [12112, 111], 'comment': '# N(1675)0 -> N(1440)0 + pi0'},
        {'br': 0.167, 'daughters': [12212, -211], 'comment': '# N(1675)0 -> N(1440)+ + pi-'},
    ],
    2216: [  # N(1675)+
        {'br': 0.267, 'daughters': [2112, 211], 'comment': '# N(1675)+ -> n + pi+'},
        {'br': 0.133, 'daughters': [2212, 111], 'comment': '# N(1675)+ -> p + pi0'},
        {'br': 0.15, 'daughters': [2224, -211], 'comment': '# N(1675)+ -> Delta(1232)++ + pi-'},
        {'br': 0.10, 'daughters': [2214, 111], 'comment': '# N(1675)+ -> Delta(1232)+ + pi0'},
        {'br': 0.05, 'daughters': [2114, 211], 'comment': '# N(1675)+ -> Delta(1232)0 + pi+'},
        {'br': 0.05, 'daughters': [2212, 9000221], 'comment': '# N(1675)+ -> p + f(0)(500)'},
        {'br': 0.167, 'daughters': [12112, 211], 'comment': '# N(1675)+ -> N(1440)0 + pi+'},
        {'br': 0.083, 'daughters': [12212, 111], 'comment': '# N(1675)+ -> N(1440)+ + pi0'},
    ],

    # N(1710)0: PDG 2025 ranges: Npi (5-20%), Neta (10-50%), LambdaK (5-25%),
    #           Deltapi (3-9%), N(1535)pi (9-21%)
    # Current: Npi 24%, Neta 30%, LambdaK 16%, Deltapi 12%, N(1535)pi 18% -> Npi, Deltapi too high
    # Updated: Npi 15%, Neta 30%, LambdaK 16%, Deltapi 6%, N(1535)pi 18%, Nsigma 15%
    42112: [  # N(1710)0
        {'br': 0.05, 'daughters': [2112, 111], 'comment': '# N(1710)0 -> n + pi0'},
        {'br': 0.10, 'daughters': [2212, -211], 'comment': '# N(1710)0 -> p + pi-'},
        {'br': 0.30, 'daughters': [2112, 221], 'comment': '# N(1710)0 -> n + eta'},
        {'br': 0.16, 'daughters': [3122, 311], 'comment': '# N(1710)0 -> Lambda + K0'},
        {'br': 0.01, 'daughters': [2214, -211], 'comment': '# N(1710)0 -> Delta(1232)+ + pi-'},
        {'br': 0.02, 'daughters': [2114, 111], 'comment': '# N(1710)0 -> Delta(1232)0 + pi0'},
        {'br': 0.03, 'daughters': [1114, 211], 'comment': '# N(1710)0 -> Delta(1232)- + pi+'},
        {'br': 0.06, 'daughters': [22112, 111], 'comment': '# N(1710)0 -> N(1535)0 + pi0'},
        {'br': 0.12, 'daughters': [22212, -211], 'comment': '# N(1710)0 -> N(1535)+ + pi-'},
        {'br': 0.15, 'daughters': [2112, 9000221], 'comment': '# N(1710)0 -> n + f(0)(500)'},
    ],
    42212: [  # N(1710)+
        {'br': 0.10, 'daughters': [2112, 211], 'comment': '# N(1710)+ -> n + pi+'},
        {'br': 0.05, 'daughters': [2212, 111], 'comment': '# N(1710)+ -> p + pi0'},
        {'br': 0.30, 'daughters': [2212, 221], 'comment': '# N(1710)+ -> p + eta'},
        {'br': 0.16, 'daughters': [3122, 321], 'comment': '# N(1710)+ -> Lambda + K+'},
        {'br': 0.03, 'daughters': [2224, -211], 'comment': '# N(1710)+ -> Delta(1232)++ + pi-'},
        {'br': 0.02, 'daughters': [2214, 111], 'comment': '# N(1710)+ -> Delta(1232)+ + pi0'},
        {'br': 0.01, 'daughters': [2114, 211], 'comment': '# N(1710)+ -> Delta(1232)0 + pi+'},
        {'br': 0.12, 'daughters': [22112, 211], 'comment': '# N(1710)+ -> N(1535)0 + pi+'},
        {'br': 0.06, 'daughters': [22212, 111], 'comment': '# N(1710)+ -> N(1535)+ + pi0'},
        {'br': 0.15, 'daughters': [2212, 9000221], 'comment': '# N(1710)+ -> p + f(0)(500)'},
    ],

    # N(1875)0: PDG 2025 ranges: Npi (3-11%), Neta (3-16%), Nomega (15-25%),
    #           Nsigma (16-60%), Deltapi (not listed separately), LambdaK (1-2%)
    # Current: Npi 9%, Neta 1%, Nomega 20%, Deltapi 30%, Nsigma 40% -> Neta too low
    # Updated: Npi 9%, Neta 8%, Nomega 20%, Deltapi 25%, Nsigma 37%, LambdaK 1%
    9902114: [  # N(1875)0
        {'br': 0.03, 'daughters': [2112, 111], 'comment': '# N(1875)0 -> n + pi0'},
        {'br': 0.06, 'daughters': [2212, -211], 'comment': '# N(1875)0 -> p + pi-'},
        {'br': 0.08, 'daughters': [2112, 221], 'comment': '# N(1875)0 -> n + eta'},
        {'br': 0.20, 'daughters': [2112, 223], 'comment': '# N(1875)0 -> n + omega(782)'},
        {'br': 0.125, 'daughters': [1114, 211], 'comment': '# N(1875)0 -> Delta(1232)- + pi+'},
        {'br': 0.083, 'daughters': [2114, 111], 'comment': '# N(1875)0 -> Delta(1232)0 + pi0'},
        {'br': 0.042, 'daughters': [2214, -211], 'comment': '# N(1875)0 -> Delta(1232)+ + pi-'},
        {'br': 0.37, 'daughters': [2112, 9000221], 'comment': '# N(1875)0 -> n + f(0)(500)'},
        {'br': 0.01, 'daughters': [3122, 311], 'comment': '# N(1875)0 -> Lambda + K0'},
    ],
    9902214: [  # N(1875)+
        {'br': 0.03, 'daughters': [2212, 111], 'comment': '# N(1875)+ -> p + pi0'},
        {'br': 0.06, 'daughters': [2112, 211], 'comment': '# N(1875)+ -> n + pi+'},
        {'br': 0.08, 'daughters': [2212, 221], 'comment': '# N(1875)+ -> p + eta'},
        {'br': 0.20, 'daughters': [2212, 223], 'comment': '# N(1875)+ -> p + omega(782)'},
        {'br': 0.083, 'daughters': [2214, 111], 'comment': '# N(1875)+ -> Delta(1232)+ + pi0'},
        {'br': 0.042, 'daughters': [2114, 211], 'comment': '# N(1875)+ -> Delta(1232)0 + pi+'},
        {'br': 0.125, 'daughters': [2224, -211], 'comment': '# N(1875)+ -> Delta(1232)++ + pi-'},
        {'br': 0.37, 'daughters': [2212, 9000221], 'comment': '# N(1875)+ -> p + f(0)(500)'},
        {'br': 0.01, 'daughters': [3122, 321], 'comment': '# N(1875)+ -> Lambda + K+'},
    ],

    # Delta(1600)-: PDG 2025 ranges: Npi (8-24%), Delta(1232)pi (72-82%), N(1440)pi (17-27%)
    # Updated: Npi 8%, Deltapi 72% (low end), N(1440)pi 20% (midpoint)
    # Total = 100%. Deltapi is constrained by sum, using lower end to fit N(1440)pi
    # Delta(1600) is P33 (I=3/2)
    # Isospin C-G for I=3/2(Q=-1) -> I=3/2 Delta + I=1 pi: Delta0 pi- (2/5), Delta- pi0 (3/5)
    31114: [  # Delta(1600)-
        {'br': 0.08, 'daughters': [2112, -211], 'comment': '# Delta(1600)- -> n + pi-'},
        {'br': 0.288, 'daughters': [2114, -211], 'comment': '# Delta(1600)- -> Delta(1232)0 + pi-'},
        {'br': 0.432, 'daughters': [1114, 111], 'comment': '# Delta(1600)- -> Delta(1232)- + pi0'},
        {'br': 0.20, 'daughters': [12112, -211], 'comment': '# Delta(1600)- -> N(1440)0 + pi-'},
    ],
    32114: [  # Delta(1600)0
        {'br': 0.053, 'daughters': [2112, 111], 'comment': '# Delta(1600)0 -> n + pi0'},
        {'br': 0.027, 'daughters': [2212, -211], 'comment': '# Delta(1600)0 -> p + pi-'},
        {'br': 0.384, 'daughters': [2214, -211], 'comment': '# Delta(1600)0 -> Delta(1232)+ + pi-'},
        {'br': 0.048, 'daughters': [2114, 111], 'comment': '# Delta(1600)0 -> Delta(1232)0 + pi0'},
        {'br': 0.288, 'daughters': [1114, 211], 'comment': '# Delta(1600)0 -> Delta(1232)- + pi+'},
        {'br': 0.133, 'daughters': [12112, 111], 'comment': '# Delta(1600)0 -> N(1440)0 + pi0'},
        {'br': 0.067, 'daughters': [12212, -211], 'comment': '# Delta(1600)0 -> N(1440)+ + pi-'},
    ],
    32214: [  # Delta(1600)+
        {'br': 0.027, 'daughters': [2112, 211], 'comment': '# Delta(1600)+ -> n + pi+'},
        {'br': 0.053, 'daughters': [2212, 111], 'comment': '# Delta(1600)+ -> p + pi0'},
        {'br': 0.288, 'daughters': [2224, -211], 'comment': '# Delta(1600)+ -> Delta(1232)++ + pi-'},
        {'br': 0.048, 'daughters': [2214, 111], 'comment': '# Delta(1600)+ -> Delta(1232)+ + pi0'},
        {'br': 0.384, 'daughters': [2114, 211], 'comment': '# Delta(1600)+ -> Delta(1232)0 + pi+'},
        {'br': 0.067, 'daughters': [12112, 211], 'comment': '# Delta(1600)+ -> N(1440)0 + pi+'},
        {'br': 0.133, 'daughters': [12212, 111], 'comment': '# Delta(1600)+ -> N(1440)+ + pi0'},
    ],
    32224: [  # Delta(1600)++
        {'br': 0.08, 'daughters': [2212, 211], 'comment': '# Delta(1600)++ -> p + pi+'},
        {'br': 0.432, 'daughters': [2224, 111], 'comment': '# Delta(1600)++ -> Delta(1232)++ + pi0'},
        {'br': 0.288, 'daughters': [2214, 211], 'comment': '# Delta(1600)++ -> Delta(1232)+ + pi+'},
        {'br': 0.20, 'daughters': [12212, 211], 'comment': '# Delta(1600)++ -> N(1440)+ + pi+'},
    ],

    # Delta(1900)-: PDG 2025 ranges: Npi (4-12%), Deltapi (30-70%), Nrho (22-60%),
    #               N(1440)pi (3-32%), Delta eta (not listed)
    # Current: Npi 12%, Deltapi 60%, N(1440)pi 27%, Delta eta 1%, Nrho MISSING
    # Updated: Npi 10%, Deltapi 40%, Nrho 30%, N(1440)pi 18%, Delta eta 2%
    11112: [  # Delta(1900)-
        {'br': 0.10, 'daughters': [2112, -211], 'comment': '# Delta(1900)- -> n + pi-'},
        {'br': 0.16, 'daughters': [2114, -211], 'comment': '# Delta(1900)- -> Delta(1232)0 + pi-'},
        {'br': 0.24, 'daughters': [1114, 111], 'comment': '# Delta(1900)- -> Delta(1232)- + pi0'},
        {'br': 0.30, 'daughters': [2112, -213], 'comment': '# Delta(1900)- -> n + rho(770)-'},
        {'br': 0.18, 'daughters': [12112, -211], 'comment': '# Delta(1900)- -> N(1440)0 + pi-'},
        {'br': 0.02, 'daughters': [1114, 221], 'comment': '# Delta(1900)- -> Delta(1232)- + eta'},
    ],
    11212: [  # Delta(1900)0
        {'br': 0.067, 'daughters': [2112, 111], 'comment': '# Delta(1900)0 -> n + pi0'},
        {'br': 0.033, 'daughters': [2212, -211], 'comment': '# Delta(1900)0 -> p + pi-'},
        {'br': 0.213, 'daughters': [2214, -211], 'comment': '# Delta(1900)0 -> Delta(1232)+ + pi-'},
        {'br': 0.027, 'daughters': [2114, 111], 'comment': '# Delta(1900)0 -> Delta(1232)0 + pi0'},
        {'br': 0.16, 'daughters': [1114, 211], 'comment': '# Delta(1900)0 -> Delta(1232)- + pi+'},
        {'br': 0.20, 'daughters': [2112, 113], 'comment': '# Delta(1900)0 -> n + rho(770)0'},
        {'br': 0.10, 'daughters': [2212, -213], 'comment': '# Delta(1900)0 -> p + rho(770)-'},
        {'br': 0.06, 'daughters': [12212, -211], 'comment': '# Delta(1900)0 -> N(1440)+ + pi-'},
        {'br': 0.12, 'daughters': [12112, 111], 'comment': '# Delta(1900)0 -> N(1440)0 + pi0'},
        {'br': 0.02, 'daughters': [2114, 221], 'comment': '# Delta(1900)0 -> Delta(1232)0 + eta'},
    ],
    12122: [  # Delta(1900)+
        {'br': 0.033, 'daughters': [2112, 211], 'comment': '# Delta(1900)+ -> n + pi+'},
        {'br': 0.067, 'daughters': [2212, 111], 'comment': '# Delta(1900)+ -> p + pi0'},
        {'br': 0.16, 'daughters': [2224, -211], 'comment': '# Delta(1900)+ -> Delta(1232)++ + pi-'},
        {'br': 0.027, 'daughters': [2214, 111], 'comment': '# Delta(1900)+ -> Delta(1232)+ + pi0'},
        {'br': 0.213, 'daughters': [2114, 211], 'comment': '# Delta(1900)+ -> Delta(1232)0 + pi+'},
        {'br': 0.10, 'daughters': [2112, 213], 'comment': '# Delta(1900)+ -> n + rho(770)+'},
        {'br': 0.20, 'daughters': [2212, 113], 'comment': '# Delta(1900)+ -> p + rho(770)0'},
        {'br': 0.12, 'daughters': [12212, 111], 'comment': '# Delta(1900)+ -> N(1440)+ + pi0'},
        {'br': 0.06, 'daughters': [12112, 211], 'comment': '# Delta(1900)+ -> N(1440)0 + pi+'},
        {'br': 0.02, 'daughters': [2214, 221], 'comment': '# Delta(1900)+ -> Delta(1232)+ + eta'},
    ],
    12222: [  # Delta(1900)++
        {'br': 0.10, 'daughters': [2212, 211], 'comment': '# Delta(1900)++ -> p + pi+'},
        {'br': 0.24, 'daughters': [2224, 111], 'comment': '# Delta(1900)++ -> Delta(1232)++ + pi0'},
        {'br': 0.16, 'daughters': [2214, 211], 'comment': '# Delta(1900)++ -> Delta(1232)+ + pi+'},
        {'br': 0.30, 'daughters': [2212, 213], 'comment': '# Delta(1900)++ -> p + rho(770)+'},
        {'br': 0.18, 'daughters': [12212, 211], 'comment': '# Delta(1900)++ -> N(1440)+ + pi+'},
        {'br': 0.02, 'daughters': [2224, 221], 'comment': '# Delta(1900)++ -> Delta(1232)++ + eta'},
    ],

    # Delta(1910)-: PDG 2025 ranges: Npi (10-30%), Deltapi (34-66%), SigmaK (4-14%),
    #               N(1440)pi (3-45%), Delta eta (5-13%)
    # Current: Npi 24%, Deltapi 60%, N(1440)pi 9%, Delta eta 7%, SigmaK MISSING
    # Updated: Npi 20%, Deltapi 50%, SigmaK 8%, N(1440)pi 15%, Delta eta 7%
    # SigmaK: I=3/2 -> I=1 Sigma + I=1/2 K
    # For Q=-1: Sigma- K0 (only option)
    21112: [  # Delta(1910)-
        {'br': 0.20, 'daughters': [2112, -211], 'comment': '# Delta(1910)- -> n + pi-'},
        {'br': 0.20, 'daughters': [2114, -211], 'comment': '# Delta(1910)- -> Delta(1232)0 + pi-'},
        {'br': 0.30, 'daughters': [1114, 111], 'comment': '# Delta(1910)- -> Delta(1232)- + pi0'},
        {'br': 0.08, 'daughters': [3112, 311], 'comment': '# Delta(1910)- -> Sigma- + K0'},
        {'br': 0.15, 'daughters': [12112, -211], 'comment': '# Delta(1910)- -> N(1440)0 + pi-'},
        {'br': 0.07, 'daughters': [1114, 221], 'comment': '# Delta(1910)- -> Delta(1232)- + eta'},
    ],
    21212: [  # Delta(1910)0
        {'br': 0.133, 'daughters': [2112, 111], 'comment': '# Delta(1910)0 -> n + pi0'},
        {'br': 0.067, 'daughters': [2212, -211], 'comment': '# Delta(1910)0 -> p + pi-'},
        {'br': 0.267, 'daughters': [2214, -211], 'comment': '# Delta(1910)0 -> Delta(1232)+ + pi-'},
        {'br': 0.033, 'daughters': [2114, 111], 'comment': '# Delta(1910)0 -> Delta(1232)0 + pi0'},
        {'br': 0.20, 'daughters': [1114, 211], 'comment': '# Delta(1910)0 -> Delta(1232)- + pi+'},
        {'br': 0.04, 'daughters': [3112, 321], 'comment': '# Delta(1910)0 -> Sigma- + K+'},
        {'br': 0.04, 'daughters': [3212, 311], 'comment': '# Delta(1910)0 -> Sigma0 + K0'},
        {'br': 0.05, 'daughters': [12212, -211], 'comment': '# Delta(1910)0 -> N(1440)+ + pi-'},
        {'br': 0.10, 'daughters': [12112, 111], 'comment': '# Delta(1910)0 -> N(1440)0 + pi0'},
        {'br': 0.07, 'daughters': [2114, 221], 'comment': '# Delta(1910)0 -> Delta(1232)0 + eta'},
    ],
    22122: [  # Delta(1910)+
        {'br': 0.067, 'daughters': [2112, 211], 'comment': '# Delta(1910)+ -> n + pi+'},
        {'br': 0.133, 'daughters': [2212, 111], 'comment': '# Delta(1910)+ -> p + pi0'},
        {'br': 0.20, 'daughters': [2224, -211], 'comment': '# Delta(1910)+ -> Delta(1232)++ + pi-'},
        {'br': 0.033, 'daughters': [2214, 111], 'comment': '# Delta(1910)+ -> Delta(1232)+ + pi0'},
        {'br': 0.267, 'daughters': [2114, 211], 'comment': '# Delta(1910)+ -> Delta(1232)0 + pi+'},
        {'br': 0.04, 'daughters': [3222, 311], 'comment': '# Delta(1910)+ -> Sigma+ + K0'},
        {'br': 0.04, 'daughters': [3212, 321], 'comment': '# Delta(1910)+ -> Sigma0 + K+'},
        {'br': 0.10, 'daughters': [12212, 111], 'comment': '# Delta(1910)+ -> N(1440)+ + pi0'},
        {'br': 0.05, 'daughters': [12112, 211], 'comment': '# Delta(1910)+ -> N(1440)0 + pi+'},
        {'br': 0.07, 'daughters': [2214, 221], 'comment': '# Delta(1910)+ -> Delta(1232)+ + eta'},
    ],
    22222: [  # Delta(1910)++
        {'br': 0.20, 'daughters': [2212, 211], 'comment': '# Delta(1910)++ -> p + pi+'},
        {'br': 0.30, 'daughters': [2224, 111], 'comment': '# Delta(1910)++ -> Delta(1232)++ + pi0'},
        {'br': 0.20, 'daughters': [2214, 211], 'comment': '# Delta(1910)++ -> Delta(1232)+ + pi+'},
        {'br': 0.08, 'daughters': [3222, 321], 'comment': '# Delta(1910)++ -> Sigma+ + K+'},
        {'br': 0.15, 'daughters': [12212, 211], 'comment': '# Delta(1910)++ -> N(1440)+ + pi+'},
        {'br': 0.07, 'daughters': [2224, 221], 'comment': '# Delta(1910)++ -> Delta(1232)++ + eta'},
    ],

    # Delta(1950)-: PDG 2025 ranges: Npi (35-45%), Deltapi (1-9%), N(1680)pi (3-9%)
    # Current: Npi 45%, Deltapi 45%, N(1680)pi 10% -> Deltapi MASSIVELY too high
    # Updated: Npi 40%, Deltapi 5%, N(1680)pi 6%, Nrho 42%, Delta eta 7%
    # PDG 2025 also indicates significant N rho contribution
    1118: [  # Delta(1950)-
        {'br': 0.40, 'daughters': [2112, -211], 'comment': '# Delta(1950)- -> n + pi-'},
        {'br': 0.02, 'daughters': [2114, -211], 'comment': '# Delta(1950)- -> Delta(1232)0 + pi-'},
        {'br': 0.03, 'daughters': [1114, 111], 'comment': '# Delta(1950)- -> Delta(1232)- + pi0'},
        {'br': 0.42, 'daughters': [2112, -213], 'comment': '# Delta(1950)- -> n + rho(770)-'},
        {'br': 0.06, 'daughters': [12116, -211], 'comment': '# Delta(1950)- -> N(1680)0 + pi-'},
        {'br': 0.07, 'daughters': [1114, 221], 'comment': '# Delta(1950)- -> Delta(1232)- + eta'},
    ],
    2118: [  # Delta(1950)0
        {'br': 0.267, 'daughters': [2112, 111], 'comment': '# Delta(1950)0 -> n + pi0'},
        {'br': 0.133, 'daughters': [2212, -211], 'comment': '# Delta(1950)0 -> p + pi-'},
        {'br': 0.027, 'daughters': [2214, -211], 'comment': '# Delta(1950)0 -> Delta(1232)+ + pi-'},
        {'br': 0.003, 'daughters': [2114, 111], 'comment': '# Delta(1950)0 -> Delta(1232)0 + pi0'},
        {'br': 0.02, 'daughters': [1114, 211], 'comment': '# Delta(1950)0 -> Delta(1232)- + pi+'},
        {'br': 0.28, 'daughters': [2112, 113], 'comment': '# Delta(1950)0 -> n + rho(770)0'},
        {'br': 0.14, 'daughters': [2212, -213], 'comment': '# Delta(1950)0 -> p + rho(770)-'},
        {'br': 0.04, 'daughters': [12116, 111], 'comment': '# Delta(1950)0 -> N(1680)0 + pi0'},
        {'br': 0.02, 'daughters': [12216, -211], 'comment': '# Delta(1950)0 -> N(1680)+ + pi-'},
        {'br': 0.07, 'daughters': [2114, 221], 'comment': '# Delta(1950)0 -> Delta(1232)0 + eta'},
    ],
    2218: [  # Delta(1950)+
        {'br': 0.133, 'daughters': [2112, 211], 'comment': '# Delta(1950)+ -> n + pi+'},
        {'br': 0.267, 'daughters': [2212, 111], 'comment': '# Delta(1950)+ -> p + pi0'},
        {'br': 0.02, 'daughters': [2224, -211], 'comment': '# Delta(1950)+ -> Delta(1232)++ + pi-'},
        {'br': 0.003, 'daughters': [2214, 111], 'comment': '# Delta(1950)+ -> Delta(1232)+ + pi0'},
        {'br': 0.027, 'daughters': [2114, 211], 'comment': '# Delta(1950)+ -> Delta(1232)0 + pi+'},
        {'br': 0.14, 'daughters': [2112, 213], 'comment': '# Delta(1950)+ -> n + rho(770)+'},
        {'br': 0.28, 'daughters': [2212, 113], 'comment': '# Delta(1950)+ -> p + rho(770)0'},
        {'br': 0.02, 'daughters': [12116, 211], 'comment': '# Delta(1950)+ -> N(1680)0 + pi+'},
        {'br': 0.04, 'daughters': [12216, 111], 'comment': '# Delta(1950)+ -> N(1680)+ + pi0'},
        {'br': 0.07, 'daughters': [2214, 221], 'comment': '# Delta(1950)+ -> Delta(1232)+ + eta'},
    ],
    2228: [  # Delta(1950)++
        {'br': 0.40, 'daughters': [2212, 211], 'comment': '# Delta(1950)++ -> p + pi+'},
        {'br': 0.03, 'daughters': [2224, 111], 'comment': '# Delta(1950)++ -> Delta(1232)++ + pi0'},
        {'br': 0.02, 'daughters': [2214, 211], 'comment': '# Delta(1950)++ -> Delta(1232)+ + pi+'},
        {'br': 0.42, 'daughters': [2212, 213], 'comment': '# Delta(1950)++ -> p + rho(770)+'},
        {'br': 0.06, 'daughters': [12216, 211], 'comment': '# Delta(1950)++ -> N(1680)+ + pi+'},
        {'br': 0.07, 'daughters': [2224, 221], 'comment': '# Delta(1950)++ -> Delta(1232)++ + eta'},
    ],

    # =====================================================================
    # New meson decay refinements based on PDG 2025 cross-check
    # =====================================================================

    # rho(770)0: PDG 2025 BFX: pi+pi- (98.918%), pi0 gamma (0.047%), pi+pi-gamma (0.99%)
    # PDG2020 had 100% pi+pi-, add radiative channels from PDG 2025
    113: [
        {'br': 0.98918, 'daughters': [211, -211], 'comment': '# rho(770)0 -> pi+ + pi-'},
        {'br': 0.00047, 'daughters': [111, 22], 'comment': '# rho(770)0 -> pi0 + gamma'},
        {'br': 0.0099, 'daughters': [211, -211, 22], 'comment': '# rho(770)0 -> pi+ + pi- + gamma'},
        {'br': 0.00045, 'daughters': [221, 22], 'comment': '# rho(770)0 -> eta + gamma'},
    ],

    # f(2)'(1525): Updated branching fractions from PDG 2025
    # PDG 2025 BFX: KKbar (87.6%), eta eta (10.3%), pi pi (0.82%)
    # Fix: eta eta was 11.5% in PDG2020, PDG2025 says 10.3%
    335: [
        {'br': 0.4396, 'daughters': [321, -321], 'comment': "# f(2)'(1525) -> K+ + K-"},
        {'br': 0.4396, 'daughters': [311, -311], 'comment': "# f(2)'(1525) -> K0 + anti-K0"},
        {'br': 0.103, 'daughters': [221, 221], 'comment': "# f(2)'(1525) -> eta + eta"},
        {'br': 0.003, 'daughters': [111, 111], 'comment': "# f(2)'(1525) -> pi0 + pi0"},
        {'br': 0.006, 'daughters': [211, -211], 'comment': "# f(2)'(1525) -> pi+ + pi-"},
        {'br': 0.009, 'daughters': [22, 22], 'comment': "# f(2)'(1525) -> gamma + gamma"},
    ],

    # f(1)(1285): Updated branching fractions from PDG 2025
    # PDG 2025 BFX: eta pi pi (34.9%), 4pi (32.8% total), KKbar pi (9%)
    # Note: a0(980)pi is a subset of eta pi pi (via a0->eta pi)
    # In FIST the eta pi+pi- channel and a0(980) pi channel overlap in final state
    # PDG 2025 also gives pi0pi0pi+pi- (21.8% from BR ratio)
    # Redistributed: eta pipi 35%, 4pi channels 33%, KKbarpi 9%, rho gamma 6%,
    # remaining ~17% in a0(980)pi which feeds into eta pipi final state
    20223: [
        {'br': 0.349, 'daughters': [221, 211, -211], 'comment': '# f(1)(1285) -> eta + pi+ + pi-'},
        {'br': 0.218, 'daughters': [111, 111, 211, -211], 'comment': '# f(1)(1285) -> pi0 + pi0 + pi+ + pi-'},
        {'br': 0.110, 'daughters': [211, -211, 113], 'comment': '# f(1)(1285) -> pi+ + pi- + rho(770)0'},
        {'br': 0.024, 'daughters': [321, -321, 111], 'comment': '# f(1)(1285) -> K+ + K- + pi0'},
        {'br': 0.024, 'daughters': [321, -311, -211], 'comment': '# f(1)(1285) -> K+ + anti-K0 + pi-'},
        {'br': 0.024, 'daughters': [311, -321, 211], 'comment': '# f(1)(1285) -> K0 + K- + pi+'},
        {'br': 0.024, 'daughters': [311, -311, 111], 'comment': '# f(1)(1285) -> K0 + anti-K0 + pi0'},
        {'br': 0.167, 'daughters': [221, 111, 111], 'comment': '# f(1)(1285) -> eta + pi0 + pi0'},
        {'br': 0.060, 'daughters': [113, 22], 'comment': '# f(1)(1285) -> rho(770)0 + gamma'},
    ],

    # =====================================================================
    # Lambda, Sigma, Xi, Omega decay updates based on PDG 2025 cross-check
    # =====================================================================

    # Lambda(1520): PDG 2025 BFX: N Kbar (45%), Sigma pi (42%), Lambda pipi (10%),
    #               Sigma(1385) pipi (0.9%), Lambda gamma (0.85%)
    # Fix: Lambda gamma was 1.0% in PDG2020, PDG2025 says 0.85%
    3124: [
        {'br': 0.230, 'daughters': [2212, -321], 'comment': '# Lambda(1520) -> p + K-'},
        {'br': 0.230, 'daughters': [2112, -311], 'comment': '# Lambda(1520) -> n + anti-K0'},
        {'br': 0.140, 'daughters': [3222, -211], 'comment': '# Lambda(1520) -> Sigma+ + pi-'},
        {'br': 0.140, 'daughters': [3212, 111], 'comment': '# Lambda(1520) -> Sigma0 + pi0'},
        {'br': 0.140, 'daughters': [3112, 211], 'comment': '# Lambda(1520) -> Sigma- + pi+'},
        {'br': 0.068, 'daughters': [3122, 211, -211], 'comment': '# Lambda(1520) -> Lambda + pi+ + pi-'},
        {'br': 0.034, 'daughters': [3122, 111, 111], 'comment': '# Lambda(1520) -> Lambda + pi0 + pi0'},
        {'br': 0.003, 'daughters': [3222, -211, 111], 'comment': '# Lambda(1520) -> Sigma+ + pi- + pi0'},
        {'br': 0.003, 'daughters': [3212, 211, -211], 'comment': '# Lambda(1520) -> Sigma0 + pi+ + pi-'},
        {'br': 0.001, 'daughters': [3212, 111, 111], 'comment': '# Lambda(1520) -> Sigma0 + pi0 + pi0'},
        {'br': 0.003, 'daughters': [3112, 211, 111], 'comment': '# Lambda(1520) -> Sigma- + pi+ + pi0'},
        {'br': 0.0085, 'daughters': [3122, 22], 'comment': '# Lambda(1520) -> Lambda + gamma'},
    ],

    # Lambda(1600): PDG 2025 BFX: N Kbar (15-30%), Sigma pi (10-60%), Lambda sigma (19%)
    # Fix: Lambda sigma was 16% in PDG2020, PDG2025 says 19%
    23122: [
        {'br': 0.105, 'daughters': [2212, -321], 'comment': '# Lambda(1600) -> p + K-'},
        {'br': 0.105, 'daughters': [2112, -311], 'comment': '# Lambda(1600) -> n + anti-K0'},
        {'br': 0.200, 'daughters': [3222, -211], 'comment': '# Lambda(1600) -> Sigma+ + pi-'},
        {'br': 0.200, 'daughters': [3212, 111], 'comment': '# Lambda(1600) -> Sigma0 + pi0'},
        {'br': 0.200, 'daughters': [3112, 211], 'comment': '# Lambda(1600) -> Sigma- + pi+'},
        {'br': 0.190, 'daughters': [3122, 9000221], 'comment': '# Lambda(1600) -> Lambda + f(0)(500)'},
    ],

    # Lambda(1670): PDG 2025 BFX: N Kbar (20-30%), Sigma pi (25-55%), Lambda eta (20%)
    # Fix: Lambda eta was 19% in PDG2020, PDG2025 says 20%
    33122: [
        {'br': 0.145, 'daughters': [2212, -321], 'comment': '# Lambda(1670) -> p + K-'},
        {'br': 0.145, 'daughters': [2112, -311], 'comment': '# Lambda(1670) -> n + anti-K0'},
        {'br': 0.170, 'daughters': [3222, -211], 'comment': '# Lambda(1670) -> Sigma+ + pi-'},
        {'br': 0.170, 'daughters': [3212, 111], 'comment': '# Lambda(1670) -> Sigma0 + pi0'},
        {'br': 0.170, 'daughters': [3112, 211], 'comment': '# Lambda(1670) -> Sigma- + pi+'},
        {'br': 0.200, 'daughters': [3122, 221], 'comment': '# Lambda(1670) -> Lambda + eta'},
    ],

    # Lambda(2100): PDG 2025 ranges: N Kbar (25-35%), Sigma pi (PDG2020 ~9%),
    #               Lambda omega (seen), N K* (seen), Sigma(1385) pi (seen)
    # Fix: N Kbar was 40% in PDG2020, PDG2025 says 25-35%; use 30%
    3128: [
        {'br': 0.150, 'daughters': [2212, -321], 'comment': '# Lambda(2100) -> p + K-'},
        {'br': 0.150, 'daughters': [2112, -311], 'comment': '# Lambda(2100) -> n + anti-K0'},
        {'br': 0.030, 'daughters': [3222, -211], 'comment': '# Lambda(2100) -> Sigma+ + pi-'},
        {'br': 0.030, 'daughters': [3212, 111], 'comment': '# Lambda(2100) -> Sigma0 + pi0'},
        {'br': 0.030, 'daughters': [3112, 211], 'comment': '# Lambda(2100) -> Sigma- + pi+'},
        {'br': 0.030, 'daughters': [3122, 221], 'comment': '# Lambda(2100) -> Lambda + eta'},
        {'br': 0.015, 'daughters': [3322, 311], 'comment': '# Lambda(2100) -> Xi0 + K0'},
        {'br': 0.015, 'daughters': [3312, 321], 'comment': '# Lambda(2100) -> Xi- + K+'},
        {'br': 0.070, 'daughters': [3122, 223], 'comment': '# Lambda(2100) -> Lambda + omega(782)'},
        {'br': 0.100, 'daughters': [2212, -323], 'comment': '# Lambda(2100) -> p + K*(892)-'},
        {'br': 0.100, 'daughters': [2112, -313], 'comment': '# Lambda(2100) -> n + anti-K*(892)0'},
        {'br': 0.060, 'daughters': [3224, -211], 'comment': '# Lambda(2100) -> Sigma(1385)+ + pi-'},
        {'br': 0.060, 'daughters': [3214, 111], 'comment': '# Lambda(2100) -> Sigma(1385)0 + pi0'},
        {'br': 0.060, 'daughters': [3114, 211], 'comment': '# Lambda(2100) -> Sigma(1385)- + pi+'},
        {'br': 0.100, 'daughters': [3122, 9000221], 'comment': '# Lambda(2100) -> Lambda + f(0)(500)'},
    ],

    # Omega-: PDG 2025 BFX: Lambda K- (67.8%), Xi0 pi- (24.3%), Xi- pi0 (8.6%)
    # Fix: Xi0 pi- was 23.6% in PDG2020, PDG2025 says 24.3%
    3334: [
        {'br': 0.671, 'daughters': [3122, -321], 'comment': '# Omega -> Lambda + K-'},
        {'br': 0.243, 'daughters': [3322, -211], 'comment': '# Omega -> Xi0 + pi-'},
        {'br': 0.086, 'daughters': [3312, 111], 'comment': '# Omega -> Xi- + pi0'},
    ],

    # =====================================================================
    # Delta(1905): PDG 2025 ranges: N pi (9-15%), Delta pi F-wave (40-58%),
    #              Delta pi P-wave (8-43%), Delta eta (2-6%)
    # Fix: Delta pi was 90% total in PDG2020 (not split by partial wave)
    # For FIST: P-wave and F-wave go to same final state (Delta+pi),
    # so we use total Delta pi = 75% (midpoints: F=49% + P=26%)
    # Updated: N pi 12%, Delta pi 75%, Delta eta 4%, N rho 9%
    # =====================================================================

    # Delta(1905)-
    1116: [
        {'br': 0.12, 'daughters': [2112, -211], 'comment': '# Delta(1905)- -> n + pi-'},
        {'br': 0.300, 'daughters': [2114, -211], 'comment': '# Delta(1905)- -> Delta(1232)0 + pi-'},
        {'br': 0.450, 'daughters': [1114, 111], 'comment': '# Delta(1905)- -> Delta(1232)- + pi0'},
        {'br': 0.09, 'daughters': [2112, -213], 'comment': '# Delta(1905)- -> n + rho(770)-'},
        {'br': 0.04, 'daughters': [1114, 221], 'comment': '# Delta(1905)- -> Delta(1232)- + eta'},
    ],
    # Delta(1905)0
    1216: [
        {'br': 0.08, 'daughters': [2112, 111], 'comment': '# Delta(1905)0 -> n + pi0'},
        {'br': 0.04, 'daughters': [2212, -211], 'comment': '# Delta(1905)0 -> p + pi-'},
        {'br': 0.400, 'daughters': [2214, -211], 'comment': '# Delta(1905)0 -> Delta(1232)+ + pi-'},
        {'br': 0.050, 'daughters': [2114, 111], 'comment': '# Delta(1905)0 -> Delta(1232)0 + pi0'},
        {'br': 0.300, 'daughters': [1114, 211], 'comment': '# Delta(1905)0 -> Delta(1232)- + pi+'},
        {'br': 0.06, 'daughters': [2112, 113], 'comment': '# Delta(1905)0 -> n + rho(770)0'},
        {'br': 0.03, 'daughters': [2212, -213], 'comment': '# Delta(1905)0 -> p + rho(770)-'},
        {'br': 0.04, 'daughters': [2114, 221], 'comment': '# Delta(1905)0 -> Delta(1232)0 + eta'},
    ],
    # Delta(1905)+
    2126: [
        {'br': 0.04, 'daughters': [2112, 211], 'comment': '# Delta(1905)+ -> n + pi+'},
        {'br': 0.08, 'daughters': [2212, 111], 'comment': '# Delta(1905)+ -> p + pi0'},
        {'br': 0.300, 'daughters': [2224, -211], 'comment': '# Delta(1905)+ -> Delta(1232)++ + pi-'},
        {'br': 0.050, 'daughters': [2214, 111], 'comment': '# Delta(1905)+ -> Delta(1232)+ + pi0'},
        {'br': 0.400, 'daughters': [2114, 211], 'comment': '# Delta(1905)+ -> Delta(1232)0 + pi+'},
        {'br': 0.03, 'daughters': [2112, 213], 'comment': '# Delta(1905)+ -> n + rho(770)+'},
        {'br': 0.06, 'daughters': [2212, 113], 'comment': '# Delta(1905)+ -> p + rho(770)0'},
        {'br': 0.04, 'daughters': [2214, 221], 'comment': '# Delta(1905)+ -> Delta(1232)+ + eta'},
    ],
    # Delta(1905)++
    2226: [
        {'br': 0.12, 'daughters': [2212, 211], 'comment': '# Delta(1905)++ -> p + pi+'},
        {'br': 0.450, 'daughters': [2224, 111], 'comment': '# Delta(1905)++ -> Delta(1232)++ + pi0'},
        {'br': 0.300, 'daughters': [2214, 211], 'comment': '# Delta(1905)++ -> Delta(1232)+ + pi+'},
        {'br': 0.09, 'daughters': [2212, 213], 'comment': '# Delta(1905)++ -> p + rho(770)+'},
        {'br': 0.04, 'daughters': [2224, 221], 'comment': '# Delta(1905)++ -> Delta(1232)++ + eta'},
    ],

    # =====================================================================
    # Delta(1920): PDG 2025 ranges: N pi (5-20%), Delta pi F-wave (44-72%),
    #              Delta pi P-wave (2-28%), Delta eta (5-17%)
    # Fix: Delta pi was 75% total in PDG2020 (not split)
    # For FIST: total Delta pi = 73% (midpoints: F=58% + P=15%)
    # Updated: N pi 13%, Delta pi 73%, Delta eta 10%, N rho 4%
    # =====================================================================

    # Delta(1920)-
    21114: [
        {'br': 0.13, 'daughters': [2112, -211], 'comment': '# Delta(1920)- -> n + pi-'},
        {'br': 0.292, 'daughters': [2114, -211], 'comment': '# Delta(1920)- -> Delta(1232)0 + pi-'},
        {'br': 0.438, 'daughters': [1114, 111], 'comment': '# Delta(1920)- -> Delta(1232)- + pi0'},
        {'br': 0.04, 'daughters': [2112, -213], 'comment': '# Delta(1920)- -> n + rho(770)-'},
        {'br': 0.10, 'daughters': [1114, 221], 'comment': '# Delta(1920)- -> Delta(1232)- + eta'},
    ],
    # Delta(1920)0
    22114: [
        {'br': 0.087, 'daughters': [2112, 111], 'comment': '# Delta(1920)0 -> n + pi0'},
        {'br': 0.043, 'daughters': [2212, -211], 'comment': '# Delta(1920)0 -> p + pi-'},
        {'br': 0.389, 'daughters': [2214, -211], 'comment': '# Delta(1920)0 -> Delta(1232)+ + pi-'},
        {'br': 0.049, 'daughters': [2114, 111], 'comment': '# Delta(1920)0 -> Delta(1232)0 + pi0'},
        {'br': 0.292, 'daughters': [1114, 211], 'comment': '# Delta(1920)0 -> Delta(1232)- + pi+'},
        {'br': 0.04, 'daughters': [2112, 113], 'comment': '# Delta(1920)0 -> n + rho(770)0'},
        {'br': 0.10, 'daughters': [2114, 221], 'comment': '# Delta(1920)0 -> Delta(1232)0 + eta'},
    ],
    # Delta(1920)+
    22214: [
        {'br': 0.043, 'daughters': [2112, 211], 'comment': '# Delta(1920)+ -> n + pi+'},
        {'br': 0.087, 'daughters': [2212, 111], 'comment': '# Delta(1920)+ -> p + pi0'},
        {'br': 0.292, 'daughters': [2224, -211], 'comment': '# Delta(1920)+ -> Delta(1232)++ + pi-'},
        {'br': 0.049, 'daughters': [2214, 111], 'comment': '# Delta(1920)+ -> Delta(1232)+ + pi0'},
        {'br': 0.389, 'daughters': [2114, 211], 'comment': '# Delta(1920)+ -> Delta(1232)0 + pi+'},
        {'br': 0.04, 'daughters': [2212, 113], 'comment': '# Delta(1920)+ -> p + rho(770)0'},
        {'br': 0.10, 'daughters': [2214, 221], 'comment': '# Delta(1920)+ -> Delta(1232)+ + eta'},
    ],
    # Delta(1920)++
    22224: [
        {'br': 0.13, 'daughters': [2212, 211], 'comment': '# Delta(1920)++ -> p + pi+'},
        {'br': 0.438, 'daughters': [2224, 111], 'comment': '# Delta(1920)++ -> Delta(1232)++ + pi0'},
        {'br': 0.292, 'daughters': [2214, 211], 'comment': '# Delta(1920)++ -> Delta(1232)+ + pi+'},
        {'br': 0.04, 'daughters': [2212, 213], 'comment': '# Delta(1920)++ -> p + rho(770)+'},
        {'br': 0.10, 'daughters': [2224, 221], 'comment': '# Delta(1920)++ -> Delta(1232)++ + eta'},
    ],

    # =====================================================================
    # N(2190): PDG 2025 ranges: N pi (10-20%), Delta pi D-wave (19-31%),
    #          N omega (8-20%), N sigma (3-9%), N eta (1-5%)
    # Fix: Delta pi was 36% (above 31%), reduce to 25% (midpoint)
    # N pi pi was too low at 15%, increase via N sigma
    # =====================================================================

    # N(2190)0
    # PDG 2025: N pi (10-20%), Delta pi D-wave (19-31%), N omega (8-20%),
    #           N sigma (3-9%), N eta (1-5%), Lambda K (~0.5%), N rho significant
    # Updated: N pi 15%, Delta pi 25%, N omega 15%, N sigma 6%, N eta 5%,
    #          N rho 33%, Lambda K 1%
    1218: [
        {'br': 0.050, 'daughters': [2112, 111], 'comment': '# N(2190)0 -> n + pi0'},
        {'br': 0.100, 'daughters': [2212, -211], 'comment': '# N(2190)0 -> p + pi-'},
        {'br': 0.150, 'daughters': [2112, 223], 'comment': '# N(2190)0 -> n + omega(782)'},
        {'br': 0.042, 'daughters': [2214, -211], 'comment': '# N(2190)0 -> Delta(1232)+ + pi-'},
        {'br': 0.083, 'daughters': [2114, 111], 'comment': '# N(2190)0 -> Delta(1232)0 + pi0'},
        {'br': 0.125, 'daughters': [1114, 211], 'comment': '# N(2190)0 -> Delta(1232)- + pi+'},
        {'br': 0.140, 'daughters': [2112, 113], 'comment': '# N(2190)0 -> n + rho(770)0'},
        {'br': 0.190, 'daughters': [2212, -213], 'comment': '# N(2190)0 -> p + rho(770)-'},
        {'br': 0.050, 'daughters': [2112, 221], 'comment': '# N(2190)0 -> n + eta'},
        {'br': 0.060, 'daughters': [2112, 9000221], 'comment': '# N(2190)0 -> n + f(0)(500)'},
        {'br': 0.010, 'daughters': [3122, 311], 'comment': '# N(2190)0 -> Lambda + K0'},
    ],
    # N(2190)+
    2128: [
        {'br': 0.100, 'daughters': [2112, 211], 'comment': '# N(2190)+ -> n + pi+'},
        {'br': 0.050, 'daughters': [2212, 111], 'comment': '# N(2190)+ -> p + pi0'},
        {'br': 0.150, 'daughters': [2212, 223], 'comment': '# N(2190)+ -> p + omega(782)'},
        {'br': 0.042, 'daughters': [2114, 211], 'comment': '# N(2190)+ -> Delta(1232)0 + pi+'},
        {'br': 0.083, 'daughters': [2214, 111], 'comment': '# N(2190)+ -> Delta(1232)+ + pi0'},
        {'br': 0.125, 'daughters': [2224, -211], 'comment': '# N(2190)+ -> Delta(1232)++ + pi-'},
        {'br': 0.190, 'daughters': [2212, 113], 'comment': '# N(2190)+ -> p + rho(770)0'},
        {'br': 0.140, 'daughters': [2112, 213], 'comment': '# N(2190)+ -> n + rho(770)+'},
        {'br': 0.050, 'daughters': [2212, 221], 'comment': '# N(2190)+ -> p + eta'},
        {'br': 0.060, 'daughters': [2212, 9000221], 'comment': '# N(2190)+ -> p + f(0)(500)'},
        {'br': 0.010, 'daughters': [3122, 321], 'comment': '# N(2190)+ -> Lambda + K+'},
    ],

    # pi2(1880)0: Add f(1)(1285)pi channel (seen in PDG 2025)
    # Updated: a0(980)eta 40%, a2(1320)eta 30%, f0(1500)pi 10%, f1(1285)pi 20%
    20115: [  # pi2(1880)0
        {'br': 0.40, 'daughters': [9000111, 221], 'comment': '# pi2(1880)0 -> a(0)(980)0 + eta'},
        {'br': 0.30, 'daughters': [115, 221], 'comment': '# pi2(1880)0 -> a(2)(1320)0 + eta'},
        {'br': 0.10, 'daughters': [9030221, 111], 'comment': '# pi2(1880)0 -> f(0)(1500) + pi0'},
        {'br': 0.20, 'daughters': [20223, 111], 'comment': '# pi2(1880)0 -> f(1)(1285) + pi0'},
    ],
    20215: [  # pi2(1880)+
        {'br': 0.40, 'daughters': [9000211, 221], 'comment': '# pi2(1880)+ -> a(0)(980)+ + eta'},
        {'br': 0.30, 'daughters': [215, 221], 'comment': '# pi2(1880)+ -> a(2)(1320)+ + eta'},
        {'br': 0.10, 'daughters': [9030221, 211], 'comment': '# pi2(1880)+ -> f(0)(1500) + pi+'},
        {'br': 0.20, 'daughters': [20223, 211], 'comment': '# pi2(1880)+ -> f(1)(1285) + pi+'},
    ],
}

# =====================================================================
# Charm decay updates from PDG2025
# These replace charm/charmonium decays passed through from PDG2020
# All values from PDG 2025 Review (BFX exclusive modes)
# For charmonia where PDG exclusive modes sum << 1.0, remaining BR
# goes to a catch-all multi-pion channel
# =====================================================================

CHARM_DECAY_UPDATES = {
    # -----------------------------------------------------------------
    # eta_c(1S):  PDG2025 exclusive sum ~62.8%. Updated from PDG2020.
    # Key changes: restructured with best PDG2025 values, added new modes
    # -----------------------------------------------------------------
    441: [
        {'br': 0.1589, 'daughters': [211, 211, -211, -211, 111, 111], 'comment': '# eta_c -> 3(pi+pi-) pi0pi0 (PDG 15.89%)'},
        {'br': 0.0456, 'daughters': [211, -211, 111, 111], 'comment': '# eta_c -> pi+pi- 2pi0 (PDG 4.55%)'},
        {'br': 0.0430, 'daughters': [221, 211, 211, -211, -211], 'comment': '# eta_c -> eta 2(pi+pi-) (PDG 4.30%)'},
        {'br': 0.0340, 'daughters': [321, -321, 211, -211, 111], 'comment': '# eta_c -> K+K- pi+pi- pi0 (PDG 3.40%)'},
        {'br': 0.0365, 'daughters': [321, -321, 111], 'comment': '# eta_c -> K+K- pi0  ~half of KKbar pi (PDG 7.13%)'},
        {'br': 0.0365, 'daughters': [311, -311, 111], 'comment': '# eta_c -> K0 Kbar0 pi0  ~half of KKbar pi'},
        {'br': 0.0200, 'daughters': [331, 211, -211], 'comment': '# eta_c -> eta\'(958) pi+pi- (PDG 2.02%)'},
        {'br': 0.0189, 'daughters': [211, 211, 211, -211, -211, -211], 'comment': '# eta_c -> 3(pi+pi-) (PDG 1.89%)'},
        {'br': 0.0162, 'daughters': [221, 211, -211], 'comment': '# eta_c -> eta pi+pi- (PDG 1.62%)'},
        {'br': 0.0135, 'daughters': [313, -313, 211, -211], 'comment': '# eta_c -> K*0 Kbar*0 pi+pi- (PDG 1.35%)'},
        {'br': 0.0108, 'daughters': [225, 225], 'comment': '# eta_c -> f2(1270) f2(1270) (PDG 1.08%)'},
        {'br': 0.0096, 'daughters': [211, 211, -211, -211], 'comment': '# eta_c -> 2(pi+pi-) (PDG 0.96%)'},
        {'br': 0.0084, 'daughters': [321, -321, 211, 211, -211, -211], 'comment': '# eta_c -> K+K- 2(pi+pi-) (PDG 0.84%)'},
        {'br': 0.0083, 'daughters': [321, -321, 211, -211], 'comment': '# eta_c -> K+K- pi+pi- (PDG 0.83%)'},
        {'br': 0.0054, 'daughters': [311, -321, 211, -211, 211], 'comment': '# eta_c -> K0 K- pi+pi-pi+ (PDG ~half of 5.39%)'},
        {'br': 0.0054, 'daughters': [-311, 321, -211, 211, -211], 'comment': '# eta_c -> Kbar0 K+ pi-pi+pi- (PDG ~half of 5.39%)'},
        {'br': 0.0037, 'daughters': [2212, -2212, 211, -211], 'comment': '# eta_c -> p pbar pi+pi- (PDG 0.37%)'},
        {'br': 0.0034, 'daughters': [2212, -2212, 111], 'comment': '# eta_c -> p pbar pi0 (PDG 0.34%)'},
        {'br': 0.0027, 'daughters': [223, 223], 'comment': '# eta_c -> omega omega (PDG 0.27%)'},
        {'br': 0.0018, 'daughters': [333, 333], 'comment': '# eta_c -> phi phi (PDG 0.18%)'},
        {'br': 0.0014, 'daughters': [321, 321, -321, -321], 'comment': '# eta_c -> 2(K+K-) (PDG 0.14%)'},
        {'br': 0.0013, 'daughters': [2212, -2212], 'comment': '# eta_c -> p pbar (PDG 0.13%)'},
        # Catch-all: remaining to generic light hadrons
        {'br': 0.5147, 'daughters': [211, -211, 211, -211, 111, 111, 111], 'comment': '# eta_c -> catch-all (unmeasured modes)'},
    ],

    # -----------------------------------------------------------------
    # J/psi(1S):  PDG2025 "hadrons" = 87.7%, exclusive sum overlaps.
    # Strategy: use exclusive modes where measured, rest to catch-all.
    # leptonic e+e- 5.97%, mu+mu- 5.96%, gamma eta_c 1.41%
    # hadronic exclusive modes sum ~20-25% of total
    # -----------------------------------------------------------------
    443: [
        {'br': 0.0597, 'daughters': [-11, 11], 'comment': '# J/psi -> e+e- (PDG 5.97%)'},
        {'br': 0.0596, 'daughters': [-13, 13], 'comment': '# J/psi -> mu+mu- (PDG 5.96%)'},
        {'br': 0.0141, 'daughters': [22, 441], 'comment': '# J/psi -> gamma eta_c (PDG 1.41%)'},
        {'br': 0.0624, 'daughters': [211, 211, -211, -211, 111, 111, 111], 'comment': '# J/psi -> 2(pi+pi-) 3pi0 (PDG 6.24%)'},
        {'br': 0.0422, 'daughters': [211, 211, -211, -211, 111], 'comment': '# J/psi -> 2(pi+pi-) pi0 (PDG 4.22%)'},
        {'br': 0.0286, 'daughters': [211, 211, 211, -211, -211, -211, 111], 'comment': '# J/psi -> 3(pi+pi-) pi0 (PDG 2.86%)'},
        {'br': 0.0200, 'daughters': [211, -211, 111], 'comment': '# J/psi -> pi+pi- pi0 (PDG 2.00%)'},
        {'br': 0.0188, 'daughters': [211, -211, 111, 111, 111], 'comment': '# J/psi -> pi+pi- 3pi0 (PDG 1.88%)'},
        {'br': 0.0161, 'daughters': [211, 211, -211, -211, 111, 111], 'comment': '# J/psi -> 2(pi+pi-) 2pi0 (PDG 1.61%)'},
        {'br': 0.0152, 'daughters': [211, -211, 111, 321, -321], 'comment': '# J/psi -> pi+pi- pi0 K+K- (PDG 1.52%)'},
        {'br': 0.0117, 'daughters': [221, 211, -211, 111], 'comment': '# J/psi -> eta pi+pi- pi0 (PDG 1.17%)'},
        {'br': 0.0112, 'daughters': [323, -323, 111], 'comment': '# J/psi -> K*(892)+ K*(892)- pi0 (PDG 1.12%)'},
        {'br': 0.0090, 'daughters': [211, 211, 211, 211, -211, -211, -211, -211, 111], 'comment': '# J/psi -> 4(pi+pi-) pi0 (PDG 0.90%)'},
        {'br': 0.0083, 'daughters': [22, 211, -211, 111, 111], 'comment': '# J/psi -> gamma pi+pi- 2pi0 (PDG 0.83%)'},
        {'br': 0.0070, 'daughters': [321, -321, 211, -211], 'comment': '# J/psi -> K+K- pi+pi- (PDG 0.70%)'},
        {'br': 0.0060, 'daughters': [2212, -2212, 211, -211], 'comment': '# J/psi -> p pbar pi+pi- (PDG 0.60%)'},
        {'br': 0.0053, 'daughters': [22, 331], 'comment': '# J/psi -> gamma eta\'(958) (PDG 0.53%)'},
        {'br': 0.0047, 'daughters': [321, -321, 211, -211, 221], 'comment': '# J/psi -> K+K- pi+pi- eta (PDG 0.47%)'},
        {'br': 0.0043, 'daughters': [211, 211, 211, -211, -211, -211], 'comment': '# J/psi -> 3(pi+pi-) (PDG 0.43%)'},
        {'br': 0.0038, 'daughters': [2112, -2112, 211, -211], 'comment': '# J/psi -> n nbar pi+pi- (PDG 0.38%)'},
        {'br': 0.0032, 'daughters': [211, 211, -211, -211], 'comment': '# J/psi -> 2(pi+pi-) (PDG 0.32%)'},
        {'br': 0.0031, 'daughters': [321, -321, 211, 211, -211, -211], 'comment': '# J/psi -> K+K- 2(pi+pi-) (PDG 0.31%)'},
        {'br': 0.0021, 'daughters': [2212, -2212], 'comment': '# J/psi -> p pbar (PDG 0.21%)'},
        {'br': 0.0021, 'daughters': [2112, -2112], 'comment': '# J/psi -> n nbar (PDG 0.21%)'},
        {'br': 0.0019, 'daughters': [3122, -3122], 'comment': '# J/psi -> Lambda Lambdabar (PDG 0.19%)'},
        {'br': 0.0016, 'daughters': [22, 225], 'comment': '# J/psi -> gamma f2(1270) (PDG 0.16%)'},
        {'br': 0.0011, 'daughters': [22, 221], 'comment': '# J/psi -> gamma eta (PDG 0.11%)'},
        # Catch-all: remaining to generic hadronic multi-body
        {'br': 0.5769, 'daughters': [211, -211, 211, -211, 111, 111], 'comment': '# J/psi -> catch-all (unmeasured hadronic modes)'},
    ],

    # -----------------------------------------------------------------
    # chi_c0(1P):  PDG2025 exclusive sum ~22%. Dominated by hadronic modes.
    # gamma J/psi = 1.41%
    # -----------------------------------------------------------------
    10441: [
        {'br': 0.0141, 'daughters': [22, 443], 'comment': '# chi_c0 -> gamma J/psi (PDG 1.41%)'},
        {'br': 0.0335, 'daughters': [211, -211, 111, 111], 'comment': '# chi_c0 -> pi+pi- 2pi0 (PDG 3.35%)'},
        {'br': 0.0218, 'daughters': [211, 211, -211, -211], 'comment': '# chi_c0 -> 2(pi+pi-) (PDG 2.18%)'},
        {'br': 0.0196, 'daughters': [211, 211, 211, -211, -211, -211], 'comment': '# chi_c0 -> 3(pi+pi-) (PDG 1.96%)'},
        {'br': 0.0181, 'daughters': [211, -211, 321, -321], 'comment': '# chi_c0 -> pi+pi- K+K- (PDG 1.81%)'},
        {'br': 0.0086, 'daughters': [321, -321, 211, -211, 111], 'comment': '# chi_c0 -> K+K- pi+pi- pi0 (PDG 0.86%)'},
        {'br': 0.0061, 'daughters': [321, -321], 'comment': '# chi_c0 -> K+K- (PDG 0.61%)'},
        {'br': 0.0126, 'daughters': [311, -311], 'comment': '# chi_c0 -> K0 Kbar0 (PDG ~2*KsKs 0.32% + KsKL) '},
        {'br': 0.0056, 'daughters': [321, -321, 111, 111], 'comment': '# chi_c0 -> K+K- 2pi0 (PDG 0.56%)'},
        {'br': 0.0033, 'daughters': [111, 111, 111, 111], 'comment': '# chi_c0 -> 4pi0 (PDG 0.33%)'},
        {'br': 0.0030, 'daughters': [221, 221], 'comment': '# chi_c0 -> eta eta (PDG 0.30%)'},
        {'br': 0.0028, 'daughters': [321, 321, -321, -321], 'comment': '# chi_c0 -> 2(K+K-) (PDG 0.28%)'},
        {'br': 0.0125, 'daughters': [321, -211, -311, 111], 'comment': '# chi_c0 -> K+ pi- Kbar0 pi0 (~half of 2.50%)'},
        {'br': 0.0125, 'daughters': [-321, 211, 311, 111], 'comment': '# chi_c0 -> K- pi+ K0 pi0 (~half of 2.50%)'},
        {'br': 0.0021, 'daughters': [2212, -2212, 211, -211], 'comment': '# chi_c0 -> p pbar pi+pi- (PDG 0.21%)'},
        {'br': 0.0017, 'daughters': [313, -313], 'comment': '# chi_c0 -> K*0 Kbar*0 (PDG 0.17%)'},
        # Catch-all: remaining ~82%
        {'br': 0.8221, 'daughters': [211, -211, 211, -211, 111, 111], 'comment': '# chi_c0 -> catch-all (unmeasured modes)'},
    ],

    # -----------------------------------------------------------------
    # chi_c1(1P):  PDG2025 exclusive sum ~44%. Dominated by gamma J/psi 34.3%.
    # -----------------------------------------------------------------
    20443: [
        {'br': 0.3426, 'daughters': [22, 443], 'comment': '# chi_c1 -> gamma J/psi (PDG 34.26%)'},
        {'br': 0.0119, 'daughters': [211, -211, 111, 111], 'comment': '# chi_c1 -> pi+pi- 2pi0 (PDG 1.19%)'},
        {'br': 0.0115, 'daughters': [321, -321, 211, -211, 111], 'comment': '# chi_c1 -> K+K- pi+pi- pi0 (PDG 1.15%)'},
        {'br': 0.0104, 'daughters': [211, 211, 211, -211, -211, -211], 'comment': '# chi_c1 -> 3(pi+pi-) (PDG 1.04%)'},
        {'br': 0.0076, 'daughters': [211, 211, -211, -211], 'comment': '# chi_c1 -> 2(pi+pi-) (PDG 0.76%)'},
        {'br': 0.0046, 'daughters': [221, 211, -211], 'comment': '# chi_c1 -> eta pi+pi- (PDG 0.46%)'},
        {'br': 0.0045, 'daughters': [211, -211, 321, -321], 'comment': '# chi_c1 -> pi+pi- K+K- (PDG 0.45%)'},
        {'br': 0.0035, 'daughters': [-11, 11, 443], 'comment': '# chi_c1 -> e+e- J/psi (PDG 0.35%)'},
        {'br': 0.0018, 'daughters': [321, -321, 111], 'comment': '# chi_c1 -> K+K- pi0 (PDG 0.18%)'},
        {'br': 0.0014, 'daughters': [313, -313], 'comment': '# chi_c1 -> K*0 Kbar*0 (PDG 0.14%)'},
        {'br': 0.0011, 'daughters': [321, -321, 111, 111], 'comment': '# chi_c1 -> K+K- 2pi0 (PDG 0.11%)'},
        {'br': 0.0011, 'daughters': [321, -321, 221, 111], 'comment': '# chi_c1 -> K+K- eta pi0 (PDG 0.11%)'},
        # Catch-all: remaining ~60%
        {'br': 0.5980, 'daughters': [211, -211, 211, -211, 111], 'comment': '# chi_c1 -> catch-all (unmeasured modes)'},
    ],

    # -----------------------------------------------------------------
    # h_c(1P):  PDG2025 exclusive sum ~66%. Dominated by gamma eta_c 60.0%.
    # -----------------------------------------------------------------
    10443: [
        {'br': 0.6000, 'daughters': [22, 441], 'comment': '# h_c -> gamma eta_c(1S) (PDG 60.0%)'},
        {'br': 0.0094, 'daughters': [211, 211, -211, -211, 111], 'comment': '# h_c -> 2(pi+pi-) pi0 (PDG 0.94%)'},
        {'br': 0.0091, 'daughters': [211, 211, 211, -211, -211, -211, 111], 'comment': '# h_c -> 3(pi+pi-) pi0 (PDG 0.91%)'},
        {'br': 0.0083, 'daughters': [211, -211, 111, 221], 'comment': '# h_c -> pi+pi- pi0 eta (PDG 0.83%)'},
        {'br': 0.0072, 'daughters': [211, 211, -211, -211, 111, 221], 'comment': '# h_c -> 2(pi+pi-) pi0 eta (PDG 0.72%)'},
        {'br': 0.0044, 'daughters': [2212, -2212, 211, -211, 111], 'comment': '# h_c -> p pbar pi+pi- pi0 (PDG 0.44%)'},
        {'br': 0.0038, 'daughters': [321, -321, 211, -211, 111], 'comment': '# h_c -> K+K- pi+pi- pi0 (PDG 0.38%)'},
        {'br': 0.0035, 'daughters': [-11, 11, 441], 'comment': '# h_c -> e+e- eta_c (PDG 0.35%)'},
        {'br': 0.0034, 'daughters': [2212, -2212, 211, -211], 'comment': '# h_c -> p pbar pi+pi- (PDG 0.34%)'},
        {'br': 0.0016, 'daughters': [211, -211, 111], 'comment': '# h_c -> pi+pi- pi0 (PDG 0.16%)'},
        {'br': 0.0014, 'daughters': [22, 331], 'comment': '# h_c -> gamma eta\'(958) (PDG 0.14%)'},
        # Catch-all: remaining ~35%
        {'br': 0.3479, 'daughters': [211, -211, 211, -211, 111, 111], 'comment': '# h_c -> catch-all (unmeasured modes)'},
    ],

    # -----------------------------------------------------------------
    # chi_c2(1P):  PDG2025 exclusive sum ~33%. Dominated by gamma J/psi 19.5%.
    # -----------------------------------------------------------------
    445: [
        {'br': 0.1949, 'daughters': [22, 443], 'comment': '# chi_c2 -> gamma J/psi (PDG 19.49%)'},
        {'br': 0.0186, 'daughters': [211, -211, 111, 111], 'comment': '# chi_c2 -> pi+pi- 2pi0 (PDG 1.86%)'},
        {'br': 0.0153, 'daughters': [211, 211, 211, -211, -211, -211], 'comment': '# chi_c2 -> 3(pi+pi-) (PDG 1.53%)'},
        {'br': 0.0117, 'daughters': [321, -321, 211, -211, 111], 'comment': '# chi_c2 -> K+K- pi+pi- pi0 (PDG 1.17%)'},
        {'br': 0.0112, 'daughters': [211, 211, -211, -211], 'comment': '# chi_c2 -> 2(pi+pi-) (PDG 1.12%)'},
        {'br': 0.0084, 'daughters': [321, -321, 211, -211], 'comment': '# chi_c2 -> K+K- pi+pi- (PDG 0.84%)'},
        {'br': 0.0070, 'daughters': [321, -211, -311, 111], 'comment': '# chi_c2 -> K+ pi- Kbar0 pi0 (~half of 1.40%)'},
        {'br': 0.0070, 'daughters': [-321, 211, 311, 111], 'comment': '# chi_c2 -> K- pi+ K0 pi0 (~half of 1.40%)'},
        {'br': 0.0022, 'daughters': [313, -313], 'comment': '# chi_c2 -> K*0 Kbar*0 (PDG 0.22%)'},
        {'br': 0.0022, 'daughters': [-11, 11, 443], 'comment': '# chi_c2 -> e+e- J/psi (PDG 0.22%)'},
        {'br': 0.0021, 'daughters': [321, -321, 111, 111], 'comment': '# chi_c2 -> K+K- 2pi0 (PDG 0.21%)'},
        {'br': 0.0017, 'daughters': [321, 321, -321, -321], 'comment': '# chi_c2 -> 2(K+K-) (PDG 0.17%)'},
        {'br': 0.0013, 'daughters': [2212, -2212, 211, -211], 'comment': '# chi_c2 -> p pbar pi+pi- (PDG 0.13%)'},
        {'br': 0.0013, 'daughters': [321, -321, 221, 111], 'comment': '# chi_c2 -> K+K- eta pi0 (PDG 0.13%)'},
        {'br': 0.0012, 'daughters': [333, 333], 'comment': '# chi_c2 -> phi phi (PDG 0.12%)'},
        {'br': 0.0011, 'daughters': [321, -321], 'comment': '# chi_c2 -> K+K- (PDG 0.10%)'},
        {'br': 0.0011, 'daughters': [111, 111, 111, 111], 'comment': '# chi_c2 -> 4pi0 (PDG 0.11%)'},
        # Catch-all: remaining ~71%
        {'br': 0.7117, 'daughters': [211, -211, 211, -211, 111, 111], 'comment': '# chi_c2 -> catch-all (unmeasured modes)'},
    ],

    # -----------------------------------------------------------------
    # eta_c(2S):  PDG2025 exclusive sum ~8%. Very poorly known.
    # -----------------------------------------------------------------
    100441: [
        {'br': 0.0190, 'daughters': [321, -321, 111], 'comment': '# eta_c(2S) -> KKbar pi ~half (PDG 1.90% total)'},
        {'br': 0.0190, 'daughters': [311, -311, 111], 'comment': '# eta_c(2S) -> KKbar pi ~half'},
        {'br': 0.0169, 'daughters': [211, 211, 211, -211, -211, -211], 'comment': '# eta_c(2S) -> 3(pi+pi-) (PDG 1.69%)'},
        {'br': 0.0147, 'daughters': [321, -321, 211, -211, 111], 'comment': '# eta_c(2S) -> K+K- pi+pi- pi0 (PDG 1.47%)'},
        {'br': 0.0055, 'daughters': [221, 211, -211], 'comment': '# eta_c(2S) -> eta pi+pi- (PDG 0.55%)'},
        # Catch-all: remaining ~92%
        {'br': 0.9249, 'daughters': [211, -211, 211, -211, 111, 111], 'comment': '# eta_c(2S) -> catch-all (unmeasured modes)'},
    ],

    # -----------------------------------------------------------------
    # psi(2S):  PDG2025 well-measured transitions.
    # J/psi pi+pi- 34.70%, J/psi pi0pi0 18.25%, J/psi eta 3.37%
    # gamma chi_c0 9.75%, gamma chi_c1 9.75%, gamma chi_c2 9.38%
    # e+e- 0.79%, mu+mu- 0.80%, gamma eta_c 0.36%
    # -----------------------------------------------------------------
    100443: [
        {'br': 0.3470, 'daughters': [443, 211, -211], 'comment': '# psi(2S) -> J/psi pi+pi- (PDG 34.70%)'},
        {'br': 0.1825, 'daughters': [443, 111, 111], 'comment': '# psi(2S) -> J/psi pi0pi0 (PDG 18.25%)'},
        {'br': 0.0975, 'daughters': [22, 10441], 'comment': '# psi(2S) -> gamma chi_c0 (PDG 9.75%)'},
        {'br': 0.0975, 'daughters': [22, 20443], 'comment': '# psi(2S) -> gamma chi_c1 (PDG 9.75%)'},
        {'br': 0.0938, 'daughters': [22, 445], 'comment': '# psi(2S) -> gamma chi_c2 (PDG 9.38%)'},
        {'br': 0.0337, 'daughters': [443, 221], 'comment': '# psi(2S) -> J/psi eta (PDG 3.37%)'},
        {'br': 0.0142, 'daughters': [211, 211, -211, -211, 111, 111, 111], 'comment': '# psi(2S) -> 2(pi+pi-) 3pi0 (PDG 1.42%)'},
        {'br': 0.0080, 'daughters': [-13, 13], 'comment': '# psi(2S) -> mu+mu- (PDG 0.80%)'},
        {'br': 0.0079, 'daughters': [-11, 11], 'comment': '# psi(2S) -> e+e- (PDG 0.79%)'},
        {'br': 0.0036, 'daughters': [22, 441], 'comment': '# psi(2S) -> gamma eta_c (PDG 0.36%)'},
        {'br': 0.0031, 'daughters': [-15, 15], 'comment': '# psi(2S) -> tau+tau- (PDG 0.31%)'},
        {'br': 0.0013, 'daughters': [443, 111], 'comment': '# psi(2S) -> J/psi pi0 (PDG 0.13%)'},
        # Catch-all: remaining ~11%
        {'br': 0.1099, 'daughters': [211, -211, 211, -211, 111], 'comment': '# psi(2S) -> catch-all (unmeasured hadronic modes)'},
    ],

    # -----------------------------------------------------------------
    # psi(3770):  PDG2025 D Dbar ~92.8%, plus radiative transitions.
    # Split D+D- / D0 Dbar0 approximately 56:44 from cross-section data.
    # -----------------------------------------------------------------
    30443: [
        {'br': 0.5200, 'daughters': [421, -421], 'comment': '# psi(3770) -> D0 Dbar0 (PDG ~52%)'},
        {'br': 0.4085, 'daughters': [411, -411], 'comment': '# psi(3770) -> D+ D- (PDG ~41%)'},
        {'br': 0.0069, 'daughters': [22, 10441], 'comment': '# psi(3770) -> gamma chi_c0 (PDG 0.69%)'},
        {'br': 0.0025, 'daughters': [22, 20443], 'comment': '# psi(3770) -> gamma chi_c1 (PDG 0.25%)'},
        {'br': 0.0019, 'daughters': [443, 211, -211], 'comment': '# psi(3770) -> J/psi pi+pi- (PDG 0.19%)'},
        # Catch-all: remaining ~6.0%
        {'br': 0.0602, 'daughters': [211, -211, 211, -211, 111], 'comment': '# psi(3770) -> catch-all non-DDbar'},
    ],

    # -----------------------------------------------------------------
    # chi_c2(3930):  Above DDbar threshold, decays mainly to DDbar.
    # PDG: mainly DDbar, no exclusive BFX data. Use same as before.
    # -----------------------------------------------------------------
    100445: [
        {'br': 0.5200, 'daughters': [421, -421], 'comment': '# chi_c2(3930) -> D0 Dbar0'},
        {'br': 0.4100, 'daughters': [411, -411], 'comment': '# chi_c2(3930) -> D+ D-'},
        {'br': 0.0700, 'daughters': [211, -211, 211, -211], 'comment': '# chi_c2(3930) -> catch-all non-DDbar'},
    ],

    # -----------------------------------------------------------------
    # psi(4040):  Above DD* threshold. Main modes: DDbar, DD*, D*D*.
    # PDG has no exclusive BFX. Use physical estimates.
    # -----------------------------------------------------------------
    9000443: [
        {'br': 0.250, 'daughters': [421, -421], 'comment': '# psi(4040) -> D0 Dbar0'},
        {'br': 0.200, 'daughters': [411, -411], 'comment': '# psi(4040) -> D+ D-'},
        {'br': 0.180, 'daughters': [423, -421], 'comment': '# psi(4040) -> D*0 Dbar0'},
        {'br': 0.180, 'daughters': [413, -411], 'comment': '# psi(4040) -> D*+ D-'},
        {'br': 0.120, 'daughters': [423, -423], 'comment': '# psi(4040) -> D*0 Dbar*0'},
        {'br': 0.070, 'daughters': [211, -211, 211, -211, 111], 'comment': '# psi(4040) -> catch-all non-charm'},
    ],

    # -----------------------------------------------------------------
    # psi(4160):  Above D_s D_s threshold.
    # Main modes: DDbar, DD*, D*D*, Ds Ds.
    # -----------------------------------------------------------------
    9010443: [
        {'br': 0.200, 'daughters': [421, -421], 'comment': '# psi(4160) -> D0 Dbar0'},
        {'br': 0.150, 'daughters': [411, -411], 'comment': '# psi(4160) -> D+ D-'},
        {'br': 0.150, 'daughters': [423, -421], 'comment': '# psi(4160) -> D*0 Dbar0'},
        {'br': 0.150, 'daughters': [413, -411], 'comment': '# psi(4160) -> D*+ D-'},
        {'br': 0.150, 'daughters': [423, -423], 'comment': '# psi(4160) -> D*0 Dbar*0'},
        {'br': 0.050, 'daughters': [431, -431], 'comment': '# psi(4160) -> Ds+ Ds-'},
        {'br': 0.050, 'daughters': [433, -431], 'comment': '# psi(4160) -> Ds*+ Ds-'},
        {'br': 0.100, 'daughters': [211, -211, 211, -211, 111], 'comment': '# psi(4160) -> catch-all non-charm'},
    ],

    # -----------------------------------------------------------------
    # psi(4415):  Highest conventional charmonium in list.
    # PDG: D Dbar2*(2460) -> D0 D- pi+ 10.5%. Rest to generic DDbar modes.
    # -----------------------------------------------------------------
    9020443: [
        {'br': 0.150, 'daughters': [421, -421], 'comment': '# psi(4415) -> D0 Dbar0'},
        {'br': 0.120, 'daughters': [411, -411], 'comment': '# psi(4415) -> D+ D-'},
        {'br': 0.100, 'daughters': [423, -421], 'comment': '# psi(4415) -> D*0 Dbar0'},
        {'br': 0.100, 'daughters': [413, -411], 'comment': '# psi(4415) -> D*+ D-'},
        {'br': 0.100, 'daughters': [423, -423], 'comment': '# psi(4415) -> D*0 Dbar*0'},
        {'br': 0.050, 'daughters': [413, -413], 'comment': '# psi(4415) -> D*+ D*-'},
        {'br': 0.050, 'daughters': [431, -431], 'comment': '# psi(4415) -> Ds+ Ds-'},
        {'br': 0.050, 'daughters': [433, -431], 'comment': '# psi(4415) -> Ds*+ Ds-'},
        {'br': 0.050, 'daughters': [433, -433], 'comment': '# psi(4415) -> Ds*+ Ds*-'},
        {'br': 0.130, 'daughters': [421, -411, 211], 'comment': '# psi(4415) -> D0 D- pi+ (incl. DDbar2*)'},
        {'br': 0.100, 'daughters': [211, -211, 211, -211, 111], 'comment': '# psi(4415) -> catch-all non-charm'},
    ],

    # -----------------------------------------------------------------
    # D*(2007)0:  PDG2025: D0 pi0 61.9%, D0 gamma 38.1%. Already correct.
    # -----------------------------------------------------------------
    423: [
        {'br': 0.619, 'daughters': [421, 111], 'comment': '# D*(2007)0 -> D0 pi0 (PDG 61.9%)'},
        {'br': 0.381, 'daughters': [421, 22], 'comment': '# D*(2007)0 -> D0 gamma (PDG 38.1%)'},
    ],

    # -----------------------------------------------------------------
    # D*(2010)+:  PDG2025: D0 pi+ 67.7%, D+ pi0 30.7%, D+ gamma 1.6%. Already correct.
    # -----------------------------------------------------------------
    413: [
        {'br': 0.677, 'daughters': [421, 211], 'comment': '# D*(2010)+ -> D0 pi+ (PDG 67.7%)'},
        {'br': 0.307, 'daughters': [411, 111], 'comment': '# D*(2010)+ -> D+ pi0 (PDG 30.7%)'},
        {'br': 0.016, 'daughters': [411, 22], 'comment': '# D*(2010)+ -> D+ gamma (PDG 1.6%)'},
    ],

    # -----------------------------------------------------------------
    # D(0)*(2300)0:  Already D+ pi-. Keep as is.
    # -----------------------------------------------------------------
    10421: [
        {'br': 1.0, 'daughters': [411, -211], 'comment': '# D(0)*(2300)0 -> D+ pi-'},
    ],

    # -----------------------------------------------------------------
    # D(0)*(2300)+:  Mirror of neutral. D0 pi+. NEW ENTRY.
    # -----------------------------------------------------------------
    10411: [
        {'br': 1.0, 'daughters': [421, 211], 'comment': '# D(0)*(2300)+ -> D0 pi+'},
    ],

    # -----------------------------------------------------------------
    # D(1)(2420)0:  D*+ pi- (dominant). Already present.
    # -----------------------------------------------------------------
    10423: [
        {'br': 0.50, 'daughters': [413, -211], 'comment': '# D(1)(2420)0 -> D*(2010)+ pi-'},
        {'br': 0.50, 'daughters': [421, 211, -211], 'comment': '# D(1)(2420)0 -> D0 pi+ pi-'},
    ],

    # -----------------------------------------------------------------
    # D(1)(2420)+:  D*0 pi+ + D+ pi+pi-. NEW ENTRY.
    # -----------------------------------------------------------------
    10413: [
        {'br': 0.50, 'daughters': [423, 211], 'comment': '# D(1)(2420)+ -> D*(2007)0 pi+'},
        {'br': 0.50, 'daughters': [411, 211, -211], 'comment': '# D(1)(2420)+ -> D+ pi+ pi-'},
    ],

    # -----------------------------------------------------------------
    # D(1)(2430)0:  Broad state. Mainly D*+ pi-. NEW ENTRY.
    # -----------------------------------------------------------------
    20423: [
        {'br': 1.0, 'daughters': [413, -211], 'comment': '# D(1)(2430)0 -> D*(2010)+ pi-'},
    ],

    # -----------------------------------------------------------------
    # D(3)*(2750)+:  No PDG data. Assume D0 pi+ dominant.
    # -----------------------------------------------------------------
    417: [
        {'br': 0.50, 'daughters': [421, 211], 'comment': '# D(3)*(2750)+ -> D0 pi+'},
        {'br': 0.50, 'daughters': [423, 211], 'comment': '# D(3)*(2750)+ -> D*0 pi+'},
    ],

    # -----------------------------------------------------------------
    # D(3)*(2750)0:  No PDG data. Assume D+ pi- dominant.
    # -----------------------------------------------------------------
    427: [
        {'br': 0.50, 'daughters': [411, -211], 'comment': '# D(3)*(2750)0 -> D+ pi-'},
        {'br': 0.50, 'daughters': [413, -211], 'comment': '# D(3)*(2750)0 -> D*+ pi-'},
    ],

    # -----------------------------------------------------------------
    # D_s3*(2860):  No PDG data. Assume DK dominant.
    # -----------------------------------------------------------------
    437: [
        {'br': 0.50, 'daughters': [421, 321], 'comment': '# D_s3*(2860) -> D0 K+'},
        {'br': 0.50, 'daughters': [411, 311], 'comment': '# D_s3*(2860) -> D+ K0'},
    ],

    # -----------------------------------------------------------------
    # D_s1*(2700):  No PDG data. Assume DK + D*K.
    # -----------------------------------------------------------------
    30433: [
        {'br': 0.50, 'daughters': [421, 321], 'comment': '# D_s1*(2700) -> D0 K+'},
        {'br': 0.50, 'daughters': [423, 321], 'comment': '# D_s1*(2700) -> D*0 K+'},
    ],

    # -----------------------------------------------------------------
    # psi(3)(3842): No PDG data. Assume DDbar dominant.
    # -----------------------------------------------------------------
    447: [
        {'br': 0.50, 'daughters': [421, -421], 'comment': '# psi(3)(3842) -> D0 Dbar0'},
        {'br': 0.50, 'daughters': [411, -411], 'comment': '# psi(3)(3842) -> D+ D-'},
    ],

    # -----------------------------------------------------------------
    # psi(2)(3823): No PDG data. Assume DDbar dominant.
    # -----------------------------------------------------------------
    30445: [
        {'br': 0.50, 'daughters': [421, -421], 'comment': '# psi(2)(3823) -> D0 Dbar0'},
        {'br': 0.50, 'daughters': [411, -411], 'comment': '# psi(2)(3823) -> D+ D-'},
    ],

    # -----------------------------------------------------------------
    # Exotic/unconventional charmonia - catch-all to DDbar or J/psi pi pi
    # -----------------------------------------------------------------
    # chi_c1(3872) aka X(3872)
    9030443: [
        {'br': 0.50, 'daughters': [443, 211, -211], 'comment': '# chi_c1(3872) -> J/psi pi+pi-'},
        {'br': 0.30, 'daughters': [443, 111, 111], 'comment': '# chi_c1(3872) -> J/psi pi0pi0'},
        {'br': 0.20, 'daughters': [421, -421], 'comment': '# chi_c1(3872) -> D0 Dbar0 (near threshold)'},
    ],

    # chi_c0(3915)
    9040441: [
        {'br': 0.50, 'daughters': [443, 223], 'comment': '# chi_c0(3915) -> J/psi omega'},
        {'br': 0.50, 'daughters': [421, -421], 'comment': '# chi_c0(3915) -> D0 Dbar0'},
    ],

    # chi_c1(4140) aka Y(4140)
    9050443: [
        {'br': 0.50, 'daughters': [443, 333], 'comment': '# chi_c1(4140) -> J/psi phi'},
        {'br': 0.50, 'daughters': [431, -431], 'comment': '# chi_c1(4140) -> Ds+ Ds-'},
    ],

    # psi(4230) aka Y(4230)
    9060443: [
        {'br': 0.30, 'daughters': [443, 211, -211], 'comment': '# psi(4230) -> J/psi pi+pi-'},
        {'br': 0.30, 'daughters': [10443, 211, -211], 'comment': '# psi(4230) -> h_c pi+pi-'},
        {'br': 0.40, 'daughters': [421, -421, 111], 'comment': '# psi(4230) -> DDbar pi'},
    ],

    # chi_c1(4274)
    9070443: [
        {'br': 0.50, 'daughters': [443, 333], 'comment': '# chi_c1(4274) -> J/psi phi'},
        {'br': 0.50, 'daughters': [431, -431], 'comment': '# chi_c1(4274) -> Ds+ Ds-'},
    ],

    # psi(4360) aka Y(4360)
    9080443: [
        {'br': 0.40, 'daughters': [443, 211, -211], 'comment': '# psi(4360) -> J/psi pi+pi-'},
        {'br': 0.30, 'daughters': [100443, 211, -211], 'comment': '# psi(4360) -> psi(2S) pi+pi-'},
        {'br': 0.30, 'daughters': [421, -421, 111], 'comment': '# psi(4360) -> DDbar pi'},
    ],

    # psi(4660) aka Y(4660)
    9090443: [
        {'br': 0.40, 'daughters': [100443, 211, -211], 'comment': '# psi(4660) -> psi(2S) pi+pi-'},
        {'br': 0.30, 'daughters': [443, 211, -211], 'comment': '# psi(4660) -> J/psi pi+pi-'},
        {'br': 0.30, 'daughters': [431, -431, 111], 'comment': '# psi(4660) -> Ds Ds pi'},
    ],

    # -----------------------------------------------------------------
    # D0:  PDG2025 exclusive modes. Explicit sum ~95.5%.
    # K0S/K0L -> Kbar0 using K0S+K0L sum where both measured, else x2.
    # Resonances (rho, K*) used where BFX1 shows clear dominance.
    # omega channels listed separately; parent inclusive modes reduced
    # to avoid double-counting (omega -> pi+pi-pi0 overlap subtracted).
    # -----------------------------------------------------------------
    421: [
        # --- Semileptonic modes ---
        {'br': 0.0354, 'daughters': [-321, -11, 12], 'comment': '# D0 -> K- e+ nu_e (PDG 3.54%)'},
        {'br': 0.0342, 'daughters': [-321, -13, 14], 'comment': '# D0 -> K- mu+ nu_mu (PDG 3.42%)'},
        {'br': 0.0216, 'daughters': [-323, -11, 12], 'comment': '# D0 -> K*(892)- e+ nu_e (PDG 2.16%)'},
        {'br': 0.0206, 'daughters': [-323, -13, 14], 'comment': '# D0 -> K*(892)- mu+ nu_mu (PDG 2.06%)'},
        {'br': 0.0160, 'daughters': [-321, 111, -11, 12], 'comment': '# D0 -> K- pi0 e+ nu_e (PDG 1.60%)'},
        {'br': 0.0144, 'daughters': [-311, -211, -11, 12], 'comment': '# D0 -> Kbar0 pi- e+ nu_e (PDG 1.44%)'},
        {'br': 0.0073, 'daughters': [-321, 111, -13, 14], 'comment': '# D0 -> K- pi0 mu+ nu_mu (PDG 0.73%)'},
        {'br': 0.0029, 'daughters': [-211, -11, 12], 'comment': '# D0 -> pi- e+ nu_e (PDG 0.29%)'},
        {'br': 0.0027, 'daughters': [-211, -13, 14], 'comment': '# D0 -> pi- mu+ nu_mu (PDG 0.27%)'},
        {'br': 0.0015, 'daughters': [-213, -11, 12], 'comment': '# D0 -> rho- e+ nu_e (PDG 0.15%)'},
        {'br': 0.0014, 'daughters': [-211, 111, -11, 12], 'comment': '# D0 -> pi- pi0 e+ nu_e (PDG 0.14%)'},
        {'br': 0.0014, 'daughters': [-213, -13, 14], 'comment': '# D0 -> rho- mu+ nu_mu (PDG 0.14%)'},
        # --- Hadronic modes ---
        {'br': 0.0395, 'daughters': [-321, 211], 'comment': '# D0 -> K- pi+ (PDG 3.95%)'},
        {'br': 0.1124, 'daughters': [-321, 213], 'comment': '# D0 -> K- rho+ (PDG BFX1 11.24%)'},
        {'br': 0.0318, 'daughters': [-321, 211, 111], 'comment': '# D0 -> K- pi+ pi0 non-rho (14.42 - 11.24)'},
        {'br': 0.0886, 'daughters': [-321, 211, 111, 111], 'comment': '# D0 -> K- pi+ 2pi0 (PDG 8.86%)'},
        {'br': 0.0686, 'daughters': [-321, 211, 113], 'comment': '# D0 -> K- pi+ rho0 (PDG BFX1 6.86%)'},
        {'br': 0.0136, 'daughters': [-321, 211, 211, -211], 'comment': '# D0 -> K- 2pi+ pi- non-rho (8.22 - 6.86)'},
        {'br': 0.0848, 'daughters': [-311, 211, -211, 111], 'comment': '# D0 -> Kbar0 pi+pi-pi0 non-omega (10.50 - 2.02)'},
        {'br': 0.0128, 'daughters': [-321, 211, 211, -211, 111], 'comment': '# D0 -> K- 2pi+ pi- pi0 non-omega (4.30 - 3.02)'},
        {'br': 0.0339, 'daughters': [-321, 211, 223], 'comment': '# D0 -> K- pi+ omega (PDG 3.39%)'},
        {'br': 0.0128, 'daughters': [-311, 113], 'comment': '# D0 -> Kbar0 rho0 (PDG K0S BFX1 0.64% x2)'},
        {'br': 0.0445, 'daughters': [-311, 211, -211], 'comment': '# D0 -> Kbar0 pi+pi- non-rho (5.73 - 1.28)'},
        {'br': 0.0227, 'daughters': [-311, 223], 'comment': '# D0 -> Kbar0 omega (PDG K0S 1.11% + K0L 1.16%)'},
        {'br': 0.0222, 'daughters': [-311, 111], 'comment': '# D0 -> Kbar0 pi0 (PDG K0S 1.24% + K0L 0.98%)'},
        {'br': 0.0202, 'daughters': [-311, 221, 111], 'comment': '# D0 -> Kbar0 eta pi0 (PDG K0S 1.01% x2)'},
        {'br': 0.0188, 'daughters': [-321, 211, 221], 'comment': '# D0 -> K- pi+ eta (PDG 1.88%) [NEW]'},
        {'br': 0.0176, 'daughters': [-311, 331], 'comment': "# D0 -> Kbar0 eta'(958) (PDG K0S 0.95% + K0L 0.81%)"},
        {'br': 0.0182, 'daughters': [-311, 111, 111], 'comment': '# D0 -> Kbar0 2pi0 (PDG K0S 0.91% x2)'},
        {'br': 0.0170, 'daughters': [-311, 111, 223], 'comment': '# D0 -> Kbar0 pi0 omega (PDG K0S 0.85% x2)'},
        {'br': 0.0152, 'daughters': [-311, 111, 111, 111], 'comment': '# D0 -> Kbar0 3pi0 (PDG K0S 0.76% x2)'},
        {'br': 0.0101, 'daughters': [213, -211], 'comment': '# D0 -> rho+ pi- (PDG BFX1 1.01%)'},
        {'br': 0.0048, 'daughters': [211, -211, 111], 'comment': '# D0 -> pi+pi-pi0 non-rho (1.49 - 1.01)'},
        {'br': 0.0100, 'daughters': [211, -211, 111, 111], 'comment': '# D0 -> pi+pi- 2pi0 (PDG 1.00%)'},
        {'br': 0.0095, 'daughters': [-321, 211, 111, 111, 111], 'comment': '# D0 -> K- pi+ 3pi0 (PDG 0.95%) [NEW]'},
        {'br': 0.0094, 'daughters': [-311, 221], 'comment': '# D0 -> Kbar0 eta (PDG K0S 0.51% + K0L 0.43%)'},
        {'br': 0.0090, 'daughters': [-311, 321, -321], 'comment': '# D0 -> Kbar0 K+K- (PDG K0S 0.45% x2)'},
        {'br': 0.0076, 'daughters': [211, 211, -211, -211], 'comment': '# D0 -> 2pi+2pi- (PDG 0.76%)'},
        {'br': 0.0068, 'daughters': [-311, -321, 211], 'comment': '# D0 -> Kbar0 K- pi+ (PDG K0S 0.34% x2)'},
        {'br': 0.0064, 'daughters': [-321, 211, 331], 'comment': "# D0 -> K- pi+ eta'(958) (PDG 0.64%)"},
        {'br': 0.0048, 'daughters': [211, 211, -211, -211, 111, 111], 'comment': '# D0 -> 2pi+2pi- 2pi0 (PDG 0.48%)'},
        {'br': 0.0045, 'daughters': [-321, 211, 111, 221], 'comment': '# D0 -> K- pi+ pi0 eta (PDG 0.45%)'},
        {'br': 0.0041, 'daughters': [321, -321], 'comment': '# D0 -> K+K- (PDG 0.41%)'},
        {'br': 0.0035, 'daughters': [211, 211, -211, -211, 111], 'comment': '# D0 -> 2pi+2pi- pi0 (PDG 0.35%)'},
        {'br': 0.0034, 'daughters': [321, -321, 111], 'comment': '# D0 -> K+K- pi0 (PDG 0.34%)'},
        {'br': 0.0032, 'daughters': [211, -211, 111, 221], 'comment': '# D0 -> pi+pi- pi0 eta (PDG 0.32%)'},
        {'br': 0.0031, 'daughters': [321, -321, 211, -211, 111], 'comment': '# D0 -> K+K- pi+pi- pi0 (PDG 0.31%)'},
        # Catch-all: remaining 4.52% to generic multi-body
        {'br': 0.0452, 'daughters': [211, -211, 211, -211, 111], 'comment': '# D0 -> catch-all (unmeasured modes)'},
    ],

    # -----------------------------------------------------------------
    # D+:  PDG2025 exclusive sum ~99.6%. Updated from PDG2020.
    # Key changes: added Kbar0 pi+2pi0, Kbar0 2pi+pi-pi0, Kbar0 pi+eta,
    # Kbar0 pi+omega, 2pi+pi-2pi0; fixed K0S->Kbar0 conversions;
    # removed duplicate Kbar0 pi+eta'; resolved omega subfractions.
    # Semileptonic: Kbar*(892)0 l nu includes K-pi+ via K* decay.
    # K0S x2 for Kbar0 unless K0L also measured (then K0S+K0L).
    # -----------------------------------------------------------------
    411: [
        # --- Semileptonic modes ---
        {'br': 0.0881, 'daughters': [-311, -11, 12], 'comment': '# D+ -> Kbar0 e+ nu_e (PDG 8.81%)'},
        {'br': 0.0868, 'daughters': [-311, -13, 14], 'comment': '# D+ -> Kbar0 mu+ nu_mu (PDG 8.68%)'},
        {'br': 0.0540, 'daughters': [-313, -11, 12], 'comment': '# D+ -> Kbar*(892)0 e+ nu_e (PDG 5.40%)'},
        {'br': 0.0527, 'daughters': [-313, -13, 14], 'comment': '# D+ -> Kbar*(892)0 mu+ nu_mu (PDG 5.27%)'},
        {'br': 0.0037, 'daughters': [111, -11, 12], 'comment': '# D+ -> pi0 e+ nu_e (PDG 0.37%)'},
        {'br': 0.0035, 'daughters': [111, -13, 14], 'comment': '# D+ -> pi0 mu+ nu_mu (PDG 0.35%)'},
        {'br': 0.0011, 'daughters': [221, -11, 12], 'comment': '# D+ -> eta e+ nu_e (PDG 0.11%)'},
        {'br': 0.0010, 'daughters': [221, -13, 14], 'comment': '# D+ -> eta mu+ nu_mu (PDG 0.10%)'},
        {'br': 0.0019, 'daughters': [113, -11, 12], 'comment': '# D+ -> rho0 e+ nu_e (PDG 0.19%)'},
        {'br': 0.0016, 'daughters': [113, -13, 14], 'comment': '# D+ -> rho0 mu+ nu_mu (PDG 0.16%)'},
        {'br': 0.0017, 'daughters': [223, -11, 12], 'comment': '# D+ -> omega e+ nu_e (PDG 0.17%)'},
        {'br': 0.0018, 'daughters': [223, -13, 14], 'comment': '# D+ -> omega mu+ nu_mu (PDG 0.18%)'},
        {'br': 0.0012, 'daughters': [-15, 16], 'comment': '# D+ -> tau+ nu_tau (PDG 0.12%)'},
        # --- Hadronic: Cabibbo-favored K+pion modes ---
        {'br': 0.0302, 'daughters': [-311, 211], 'comment': '# D+ -> Kbar0 pi+ (PDG K0S+K0L = 3.02%)'},
        {'br': 0.0834, 'daughters': [-321, 211, 211], 'comment': '# D+ -> K- 2pi+ non-K* (PDG 9.38-1.04%)'},
        {'br': 0.0104, 'daughters': [-313, 211], 'comment': '# D+ -> Kbar*(892)0 pi+ (PDG BFX1 1.04%)'},
        {'br': 0.1228, 'daughters': [-311, 213], 'comment': '# D+ -> Kbar0 rho+ (PDG BFX1 K0S rho+ x2 = 12.28%)'},
        {'br': 0.0244, 'daughters': [-311, 211, 111], 'comment': '# D+ -> Kbar0 pi+ pi0 non-rho (14.72-12.28%)'},
        {'br': 0.0625, 'daughters': [-321, 211, 211, 111], 'comment': '# D+ -> K- 2pi+ pi0 (PDG 6.25%)'},
        {'br': 0.0620, 'daughters': [-311, 211, 211, -211], 'comment': '# D+ -> Kbar0 2pi+ pi- (PDG K0Sx2 = 6.20%)'},
        {'br': 0.0578, 'daughters': [-311, 211, 111, 111], 'comment': '# D+ -> Kbar0 pi+ 2pi0 (PDG K0Sx2 = 5.78%) [NEW]'},
        {'br': 0.0141, 'daughters': [-311, 211, 223], 'comment': '# D+ -> Kbar0 pi+ omega (PDG K0Sx2 = 1.41%) [NEW]'},
        {'br': 0.0165, 'daughters': [-311, 211, 211, -211, 111], 'comment': '# D+ -> Kbar0 2pi+pi-pi0 non-omega (3.06-1.41) [NEW]'},
        {'br': 0.0254, 'daughters': [-311, 211, 221], 'comment': '# D+ -> Kbar0 pi+ eta (PDG K0Sx2 = 2.54%) [NEW]'},
        {'br': 0.0111, 'daughters': [-311, 211, 111, 111, 111], 'comment': '# D+ -> Kbar0 pi+ 3pi0 (PDG K0Sx2 = 1.11%)'},
        {'br': 0.0057, 'daughters': [-321, 211, 211, 211, -211], 'comment': '# D+ -> K- 3pi+ pi- (PDG 0.57%)'},
        {'br': 0.0050, 'daughters': [-321, 211, 211, 111, 111], 'comment': '# D+ -> K- 2pi+ 2pi0 (PDG 0.50%)'},
        {'br': 0.0014, 'daughters': [-321, 211, 211, 221], 'comment': '# D+ -> K- 2pi+ eta (PDG 0.14%)'},
        # --- Hadronic: Cabibbo-suppressed KK/phi modes ---
        {'br': 0.0230, 'daughters': [333, 211, 111], 'comment': '# D+ -> phi pi+ pi0 (PDG 2.30%)'},
        {'br': 0.0150, 'daughters': [321, -321, 211, 111], 'comment': '# D+ -> K+K- pi+ pi0 non-phi (PDG 1.50%)'},
        {'br': 0.0057, 'daughters': [333, 211], 'comment': '# D+ -> phi pi+ (PDG 0.57%)'},
        {'br': 0.0040, 'daughters': [321, -321, 211], 'comment': '# D+ -> K+K- pi+ non-phi (0.97-0.57%)'},
        {'br': 0.0120, 'daughters': [-311, -311, 211], 'comment': '# D+ -> Kbar0 Kbar0 pi+ (PDG K0SK0Sx4 = 1.20%)'},
        {'br': 0.0062, 'daughters': [321, -311], 'comment': '# D+ -> K+ Kbar0 (PDG K0S+K0L = 0.62%)'},
        {'br': 0.0103, 'daughters': [-311, 321, 111], 'comment': '# D+ -> Kbar0 K+ pi0 (PDG K0S+K0L = 1.03%)'},
        {'br': 0.0038, 'daughters': [321, -311, 211, -211], 'comment': '# D+ -> K+ Kbar0 pi+pi- (PDG K0Sx2 = 0.38%)'},
        {'br': 0.0046, 'daughters': [-311, -321, 211, 211], 'comment': '# D+ -> Kbar0 K- 2pi+ (PDG K0Sx2 = 0.46%)'},
        {'br': 0.0051, 'daughters': [321, -311, -311], 'comment': '# D+ -> K+ Kbar0 Kbar0 (PDG ~0.51%)'},
        # --- Hadronic: Pion-only modes ---
        {'br': 0.0012, 'daughters': [211, 111], 'comment': '# D+ -> pi+ pi0 (PDG 0.12%)'},
        {'br': 0.0046, 'daughters': [211, 111, 111], 'comment': '# D+ -> pi+ 2pi0 (PDG 0.46%)'},
        {'br': 0.0033, 'daughters': [211, 211, -211], 'comment': '# D+ -> 2pi+ pi- (PDG 0.33%)'},
        {'br': 0.0042, 'daughters': [211, 111, 111, 111], 'comment': '# D+ -> pi+ 3pi0 (PDG 0.42%)'},
        {'br': 0.0020, 'daughters': [211, 111, 111, 111, 111], 'comment': '# D+ -> pi+ 4pi0 (PDG 0.20%)'},
        {'br': 0.0117, 'daughters': [211, 211, -211, 111], 'comment': '# D+ -> 2pi+pi-pi0 (PDG 1.17%)'},
        {'br': 0.0068, 'daughters': [211, 211, -211, 111, 111], 'comment': '# D+ -> 2pi+pi-2pi0 non-omega (1.07-0.39%)'},
        {'br': 0.0039, 'daughters': [223, 211, 111], 'comment': '# D+ -> omega pi+ pi0 (PDG 0.39%)'},
        {'br': 0.0017, 'daughters': [211, 211, 211, -211, -211], 'comment': '# D+ -> 3pi+2pi- (PDG 0.17%)'},
        {'br': 0.0034, 'daughters': [211, 211, -211, 111, 111, 111], 'comment': '# D+ -> 2pi+pi-3pi0 (PDG 0.34%)'},
        {'br': 0.0023, 'daughters': [211, 211, 211, -211, -211, 111], 'comment': '# D+ -> 3pi+2pi-pi0 (PDG 0.23%)'},
        # --- Hadronic: Eta/eta' modes ---
        {'br': 0.0038, 'daughters': [221, 211], 'comment': '# D+ -> eta pi+ (PDG 0.38%)'},
        {'br': 0.0021, 'daughters': [221, 211, 111], 'comment': '# D+ -> eta pi+ pi0 (PDG 0.21%)'},
        {'br': 0.0034, 'daughters': [221, 211, 211, -211], 'comment': '# D+ -> eta 2pi+pi- (PDG 0.34%)'},
        {'br': 0.0032, 'daughters': [221, 211, 111, 111], 'comment': '# D+ -> eta pi+ 2pi0 (PDG 0.32%)'},
        {'br': 0.0039, 'daughters': [221, 211, 211, -211, 111], 'comment': '# D+ -> eta 2pi+pi-pi0 (PDG 0.39%)'},
        {'br': 0.0030, 'daughters': [221, 221, 211], 'comment': '# D+ -> eta eta pi+ (PDG 0.30%)'},
        {'br': 0.0050, 'daughters': [331, 211], 'comment': "# D+ -> eta' pi+ (PDG 0.50%)"},
        {'br': 0.0038, 'daughters': [-311, 211, 331], 'comment': "# D+ -> Kbar0 pi+ eta' (PDG K0Sx2 = 0.38%)"},
        {'br': 0.0016, 'daughters': [331, 211, 111], 'comment': "# D+ -> eta' pi+ pi0 (PDG 0.16%)"},
        # Catch-all: remaining to generic Cabibbo-favored multi-body
        {'br': 0.0036, 'daughters': [-321, 211, 211, 211, -211, 111], 'comment': '# D+ -> catch-all'},
    ],

    # -----------------------------------------------------------------
    # Ds+:  PDG2025 exclusive BFX sum ~75%.
    # Strategy: use BFX final-state modes directly to avoid double-counting.
    # Uses resonances (rho+, Kbar*(892)0, phi, omega) only when they appear
    # as explicit BFX modes. Total inclusive modes (K+K-pi+, K+K-pi+pi0)
    # absorb sub-fractions (phi pi+, phi rho+) to prevent double-counting.
    # K0S modes converted: Kbar0 with BR x2. Leptonic modes included.
    # -----------------------------------------------------------------
    431: [
        {'br': 0.1267, 'daughters': [321, -313], 'comment': '# Ds+ -> K+ Kbar*(892)0 (PDG 12.67%)'},
        {'br': 0.0890, 'daughters': [221, 213], 'comment': '# Ds+ -> eta rho+ (PDG 8.90%, dominates eta pi+pi0)'},
        {'br': 0.0580, 'daughters': [331, 213], 'comment': "# Ds+ -> eta' rho+ (PDG 5.80%, dominates eta' pi+pi0)"},
        {'br': 0.0553, 'daughters': [321, -321, 211, 111], 'comment': '# Ds+ -> K+K- pi+ pi0 (PDG 5.53%, incl. phi rho+)'},
        {'br': 0.0545, 'daughters': [321, -321, 211], 'comment': '# Ds+ -> K+K- pi+ (PDG 5.45%, incl. phi pi+)'},
        {'br': 0.0539, 'daughters': [-15, 16], 'comment': '# Ds+ -> tau+ nu_tau (PDG 5.39%)'},
        {'br': 0.0490, 'daughters': [211, 211, 211, -211, -211, 111], 'comment': '# Ds+ -> 3pi+ 2pi- pi0 (PDG 4.90%)'},
        {'br': 0.0395, 'daughters': [331, 211], 'comment': "# Ds+ -> eta' pi+ (PDG 3.95%)"},
        {'br': 0.0314, 'daughters': [-311, -321, 211, 211], 'comment': '# Ds+ -> Kbar0 K- 2pi+ (PDG 2x1.57%)'},
        {'br': 0.0308, 'daughters': [211, 211, -211, 221], 'comment': '# Ds+ -> 2pi+ pi- eta (PDG 3.08%) [NEW]'},
        {'br': 0.0295, 'daughters': [321, -311], 'comment': '# Ds+ -> K+ Kbar0 (PDG 2.95%)'},
        {'br': 0.0294, 'daughters': [321, -311, 111], 'comment': '# Ds+ -> K+ Kbar0 pi0 (PDG 2x1.47%)'},
        {'br': 0.0278, 'daughters': [223, 211, 111], 'comment': '# Ds+ -> omega pi+ pi0 (PDG 2.78%)'},
        {'br': 0.0234, 'daughters': [333, -11, 12], 'comment': '# Ds+ -> phi e+ nu_e (PDG 2.34%)'},
        {'br': 0.0227, 'daughters': [221, -11, 12], 'comment': '# Ds+ -> eta e+ nu_e (PDG 2.27%)'},
        {'br': 0.0224, 'daughters': [221, -13, 14], 'comment': '# Ds+ -> eta mu+ nu_mu (PDG 2.24%)'},
        {'br': 0.0224, 'daughters': [333, -13, 14], 'comment': '# Ds+ -> phi mu+ nu_mu (PDG 2.24%)'},
        {'br': 0.0187, 'daughters': [321, -311, 211, -211], 'comment': '# Ds+ -> K+ Kbar0 pi+pi- (PDG 2x0.93%)'},
        {'br': 0.0158, 'daughters': [223, 211, 211, -211], 'comment': '# Ds+ -> omega 2pi+ pi- (PDG 1.58%)'},
        {'br': 0.0146, 'daughters': [311, -311, 211], 'comment': '# Ds+ -> K0 Kbar0 pi+ (PDG 2x0.73%)'},
        {'br': 0.0121, 'daughters': [333, 211, 211, -211], 'comment': '# Ds+ -> phi 2pi+ pi- (PDG 1.21%)'},
        {'br': 0.0109, 'daughters': [211, 211, -211], 'comment': '# Ds+ -> 2pi+ pi- (PDG 1.09%)'},
        {'br': 0.0102, 'daughters': [311, 211, 111], 'comment': '# Ds+ -> K0 pi+ pi0 (PDG 2x0.51%)'},
        {'br': 0.0097, 'daughters': [321, 211, -211, 111], 'comment': '# Ds+ -> K+ pi+pi- pi0 (PDG 0.97%)'},
        {'br': 0.0082, 'daughters': [321, 223, 111], 'comment': '# Ds+ -> K+ omega pi0 (PDG 0.82%)'},
        {'br': 0.0081, 'daughters': [331, -11, 12], 'comment': "# Ds+ -> eta' e+ nu_e (PDG 0.81%)"},
        {'br': 0.0080, 'daughters': [331, -13, 14], 'comment': "# Ds+ -> eta' mu+ nu_mu (PDG 0.80%)"},
        {'br': 0.0080, 'daughters': [211, 211, 211, -211, -211], 'comment': '# Ds+ -> 3pi+ 2pi- (PDG 0.80%)'},
        {'br': 0.0079, 'daughters': [321, 223, 221], 'comment': '# Ds+ -> K+ omega eta (PDG 0.79%)'},
        {'br': 0.0066, 'daughters': [321, -321, 211, 211, -211], 'comment': '# Ds+ -> K+K- 2pi+ pi- (PDG 0.66%)'},
        {'br': 0.0062, 'daughters': [321, 211, -211], 'comment': '# Ds+ -> K+ pi+pi- (PDG 0.62%)'},
        {'br': 0.0054, 'daughters': [223, 221, 211], 'comment': '# Ds+ -> omega eta pi+ (PDG 0.54%)'},
        {'br': 0.0054, 'daughters': [321, 223, 211, -211], 'comment': '# Ds+ -> K+ omega pi+pi- (PDG 0.54%)'},
        {'br': 0.0054, 'daughters': [-13, 14], 'comment': '# Ds+ -> mu+ nu_mu (PDG 0.54%) [NEW]'},
        {'br': 0.0056, 'daughters': [311, 211, 211, -211], 'comment': '# Ds+ -> K0 2pi+ pi- (PDG 2x0.28%)'},
        {'br': 0.0052, 'daughters': [211, 111, 111], 'comment': '# Ds+ -> pi+ 2pi0 (PDG 0.52%)'},
        {'br': 0.0029, 'daughters': [311, -11, 12], 'comment': '# Ds+ -> K0 e+ nu_e (PDG 0.29%)'},
        {'br': 0.0027, 'daughters': [321, 331], 'comment': "# Ds+ -> K+ eta' (PDG 0.27%)"},
        {'br': 0.0024, 'daughters': [311, 211], 'comment': '# Ds+ -> K0 pi+ (PDG 2x0.12%)'},
        {'br': 0.0021, 'daughters': [313, -11, 12], 'comment': '# Ds+ -> K*(892)0 e+ nu_e (PDG 0.21%)'},
        {'br': 0.0020, 'daughters': [223, -11, 12], 'comment': '# Ds+ -> omega e+ nu_e (PDG 0.20%)'},
        {'br': 0.0018, 'daughters': [321, 221], 'comment': '# Ds+ -> K+ eta (PDG 0.18%)'},
        {'br': 0.0012, 'daughters': [2212, -2112], 'comment': '# Ds+ -> p nbar (PDG 0.12%)'},
        {'br': 0.0010, 'daughters': [321, 223], 'comment': '# Ds+ -> K+ omega (PDG 0.10%)'},
        # Catch-all: remaining to generic multi-body
        {'br': 0.0462, 'daughters': [211, 211, -211, 111, 111], 'comment': '# Ds+ -> catch-all (unmeasured modes)'},
    ],

    # -----------------------------------------------------------------
    # Excited charm baryons without decays. Add simple dominant modes.
    # -----------------------------------------------------------------
    # Lambda_c(2860)+ -> D0 p (seen), D+ n
    204124: [
        {'br': 0.50, 'daughters': [421, 2212], 'comment': '# Lambda_c(2860)+ -> D0 p'},
        {'br': 0.50, 'daughters': [411, 2112], 'comment': '# Lambda_c(2860)+ -> D+ n'},
    ],

    # Lambda_c(2940)+ -> D0 p (seen)
    304124: [
        {'br': 0.50, 'daughters': [421, 2212], 'comment': '# Lambda_c(2940)+ -> D0 p'},
        {'br': 0.50, 'daughters': [411, 2112], 'comment': '# Lambda_c(2940)+ -> D+ n'},
    ],

    # Sigma_c(2800)0 -> Lambda_c+ pi-
    9004114: [
        {'br': 1.0, 'daughters': [4122, -211], 'comment': '# Sigma_c(2800)0 -> Lambda_c+ pi-'},
    ],

    # Sigma_c(2800)+ -> Lambda_c+ pi0
    9004214: [
        {'br': 1.0, 'daughters': [4122, 111], 'comment': '# Sigma_c(2800)+ -> Lambda_c+ pi0'},
    ],

    # Sigma_c(2800)++ -> Lambda_c+ pi+
    9004224: [
        {'br': 1.0, 'daughters': [4122, 211], 'comment': '# Sigma_c(2800)++ -> Lambda_c+ pi+'},
    ],

    # Xi_c(2970)0 -> Xi_c0 pi+pi- (or Xi_c'0 pi)
    204312: [
        {'br': 0.50, 'daughters': [4132, 211, -211], 'comment': '# Xi_c(2970)0 -> Xi_c0 pi+pi-'},
        {'br': 0.50, 'daughters': [4312, 111], 'comment': '# Xi_c(2970)0 -> Xi_c\'0 pi0'},
    ],

    # Xi_c(2970)+ -> Xi_c+ pi+pi- (or Xi_c'+ pi)
    204322: [
        {'br': 0.50, 'daughters': [4232, 211, -211], 'comment': '# Xi_c(2970)+ -> Xi_c+ pi+pi-'},
        {'br': 0.50, 'daughters': [4322, 111], 'comment': '# Xi_c(2970)+ -> Xi_c\'+ pi0'},
    ],

    # Xi_c(3055)+ -> Lambda_c+ K- pi+ (or D+ Lambda)
    9104324: [
        {'br': 0.50, 'daughters': [4122, -321, 211], 'comment': '# Xi_c(3055)+ -> Lambda_c+ K- pi+'},
        {'br': 0.50, 'daughters': [411, 3122], 'comment': '# Xi_c(3055)+ -> D+ Lambda'},
    ],

    # Xi_c(3080)0 -> Lambda_c+ K- (seen)
    9204316: [
        {'br': 0.50, 'daughters': [4122, -321], 'comment': '# Xi_c(3080)0 -> Lambda_c+ K-'},
        {'br': 0.50, 'daughters': [421, 3312], 'comment': '# Xi_c(3080)0 -> D0 Xi-'},
    ],

    # Xi_c(3080)+ -> Lambda_c+ K0bar pi+ (seen)
    9204326: [
        {'br': 0.50, 'daughters': [4122, -311, 211], 'comment': '# Xi_c(3080)+ -> Lambda_c+ Kbar0 pi+'},
        {'br': 0.50, 'daughters': [411, 3322], 'comment': '# Xi_c(3080)+ -> D+ Xi0'},
    ],

    # Omega_c excited states -> Omega_c pi or Xi_c K
    9104332: [  # Omega_c(3000)
        {'br': 1.0, 'daughters': [4132, 321], 'comment': '# Omega_c(3000) -> Xi_c0 K+'},
    ],
    9104334: [  # Omega_c(3065)
        {'br': 1.0, 'daughters': [4132, 321], 'comment': '# Omega_c(3065) -> Xi_c0 K+'},
    ],
    9104336: [  # Omega_c(3120)
        {'br': 1.0, 'daughters': [4132, 321], 'comment': '# Omega_c(3120) -> Xi_c0 K+'},
    ],
    9204332: [  # Omega_c(3050)
        {'br': 1.0, 'daughters': [4132, 321], 'comment': '# Omega_c(3050) -> Xi_c0 K+'},
    ],
    9204334: [  # Omega_c(3090)
        {'br': 1.0, 'daughters': [4132, 321], 'comment': '# Omega_c(3090) -> Xi_c0 K+'},
    ],
    9304332: [  # Omega_c(3185)
        {'br': 1.0, 'daughters': [4132, 321], 'comment': '# Omega_c(3185) -> Xi_c0 K+'},
    ],
    9304334: [  # Omega_c(3327)
        {'br': 1.0, 'daughters': [4132, 321], 'comment': '# Omega_c(3327) -> Xi_c0 K+'},
    ],
}

# =====================================================================
# Mass/width updates from PDG2025
# Map: PDG ID -> (new_mass, new_width)
# Only include changes for particles remaining in the list
# =====================================================================

def get_pdg2025_updates():
    """Get mass and width updates from PDG 2025 data"""
    pdg2025 = parse_pdg2025_masswidth(PDG2025_MASSWIDTH)
    return pdg2025

ELEMENTARY_MASSES = {22: 0, 11: 0.000511, -11: 0.000511, 13: 0.10566, -13: 0.10566,
                     14: 0, -14: 0, 12: 0, -12: 0, 15: 1.777, -15: 1.777, 16: 0, -16: 0}

def compute_threshold(pid, decay_channels, particles_by_id):
    """Compute the lightest decay threshold for a particle from its decay channels"""
    min_threshold = float('inf')
    for ch in decay_channels:
        ch_mass = 0
        for d in ch['daughters']:
            if abs(d) in ELEMENTARY_MASSES:
                ch_mass += ELEMENTARY_MASSES.get(d, ELEMENTARY_MASSES.get(abs(d), 0))
            elif d in particles_by_id:
                ch_mass += particles_by_id[d]['mass']
            elif -d in particles_by_id:
                ch_mass += particles_by_id[-d]['mass']
        min_threshold = min(min_threshold, ch_mass)
    return min_threshold if min_threshold < float('inf') else 0

# =====================================================================
# Main generation logic
# =====================================================================

def generate_pdg2025():
    # Parse PDG2020
    particles_2020 = parse_pdg2020_list(PDG2020_LIST)
    decays_2020 = parse_pdg2020_decays(PDG2020_DECAYS)

    # Build a name map from all available lists (includes charm, nuclei)
    # Used for descriptive comments in decays.dat
    pid_name_map = {}
    for extra_list in [PDG2020_LIST_ALL,
                       os.path.join(OUTPUT_DIR, "list-withcharm.dat")]:
        try:
            for p in parse_pdg2020_list(extra_list):
                pid_name_map[p['pdgid']] = p['name']
        except FileNotFoundError:
            pass

    # Parse PDG2025 mass-width
    pdg2025 = get_pdg2025_updates()

    # Build updated particle list
    particles_2025 = []
    removed_names = []

    # Track which PDG IDs from 2020 are kept
    kept_ids = set()

    for p in particles_2020:
        pid = p['pdgid']

        # Check if this particle should be removed
        if pid in REMOVE_IDS:
            removed_names.append(f"  Removed: {p['name']} (ID {pid})")
            continue

        # Update mass and width from PDG2025 if available
        if pid in pdg2025:
            entry = pdg2025[pid]
            if entry['mass'] is not None:
                old_mass = p['mass']
                old_width = p['width']
                p['mass'] = entry['mass']
                if entry['width'] is not None and entry['width'] > 0:
                    p['width'] = entry['width']
                elif old_width > 0 and entry['width'] == 0:
                    # Keep old width if new one is 0 (means stable/long-lived)
                    pass

        # Special mass update for rho(770)0 and rho(770)+
        # PDG2025 gives different masses for neutral and charged rho
        if pid == 113 and 113 in pdg2025:
            p['mass'] = pdg2025[113]['mass']
            p['width'] = pdg2025[113]['width']
        elif pid == 213 and 213 in pdg2025:
            p['mass'] = pdg2025[213]['mass']
            p['width'] = pdg2025[213]['width']

        # Special: K*(892) has different masses for 0 and + in PDG2025
        if pid == 323 and 323 in pdg2025:
            p['mass'] = pdg2025[323]['mass']
            p['width'] = pdg2025[323]['width']
        elif pid == 313 and 313 in pdg2025:
            p['mass'] = pdg2025[313]['mass']
            p['width'] = pdg2025[313]['width']

        # Special: Sigma(1385) has separate masses by charge in PDG2025
        if pid == 3224 and 3224 in pdg2025:
            p['mass'] = pdg2025[3224]['mass']
            p['width'] = pdg2025[3224]['width']
        elif pid == 3214 and 3214 in pdg2025:
            p['mass'] = pdg2025[3214]['mass']
            p['width'] = pdg2025[3214]['width']
        elif pid == 3114 and 3114 in pdg2025:
            p['mass'] = pdg2025[3114]['mass']
            p['width'] = pdg2025[3114]['width']

        # Apply additional mass/width updates for particles not in mass_width_2025.txt
        if pid in ADDITIONAL_UPDATES:
            upd = ADDITIONAL_UPDATES[pid]
            if 'mass' in upd:
                p['mass'] = upd['mass']
            if 'width' in upd:
                p['width'] = upd['width']

        particles_2025.append(p)
        kept_ids.add(pid)

    # Add new mesons (insert at correct position by mass)
    for new_p in NEW_MESONS + NEW_STRANGE_MESONS:
        particles_2025.append(new_p)

    # Update name map with 2025 particles
    for p in particles_2025:
        pid_name_map[p['pdgid']] = p['name']

    # Sort particles: group by sector, then by mass within sector
    # Sectors: light unflavored mesons (B=0, S=0, C=0), strange mesons (B=0, |S|>0),
    # N baryons, Delta baryons, Lambda, Sigma, Xi, Omega

    def particle_sort_key(p):
        pid = abs(p['pdgid'])
        B = p['B']
        S = p['S']
        Q = p['Q']
        mass = p['mass']

        # Light unflavored mesons
        if B == 0 and S == 0 and p['C'] == 0:
            return (0, mass, pid)
        # Strange mesons
        elif B == 0 and (S != 0 or p['absS'] > 0) and p['C'] == 0:
            return (1, mass, pid)
        # Nucleons and N*
        elif B == 1 and S == 0 and Q >= 0 and Q <= 1:
            # Check if it's a Delta (Q can be -1, 0, 1, 2)
            # N baryons have I=1/2, Delta has I=3/2
            # Use PDG ID ranges to distinguish
            if pid in [2212, 2112] or (pid >= 12112 and pid < 100000 and Q in [0,1]):
                # Could be N or Delta - need to check degeneracy
                if p['degeneracy'] in [2, 4, 6, 8, 10, 12]:
                    # Check PDG ID pattern for Delta
                    is_delta = False
                    delta_ids_base = [1114, 2114, 2214, 2224,
                                      31114, 32114, 32214, 32224,
                                      1112, 1212, 2122, 2222,
                                      11114, 12114, 12214, 12224,
                                      11112, 11212, 12122, 12222,
                                      1116, 1216, 2126, 2226,
                                      21112, 21212, 22122, 22222,
                                      21114, 22114, 22214, 22224,
                                      11116, 11216, 12126, 12226,
                                      1118, 2118, 2218, 2228]
                    if pid in delta_ids_base:
                        is_delta = True
                    if is_delta:
                        return (3, mass, pid)  # Delta
                return (2, mass, pid)  # N
            return (2, mass, pid)
        # Delta baryons
        elif B == 1 and S == 0 and (Q == -1 or Q == 2):
            return (3, mass, pid)
        # Lambda (S=-1, Q=0, I=0)
        elif B == 1 and S == -1 and Q == 0 and p['degeneracy'] in [2, 4, 6, 8, 10]:
            return (4, mass, pid)
        # Sigma (S=-1)
        elif B == 1 and S == -1:
            return (5, mass, pid)
        # Xi (S=-2)
        elif B == 1 and S == -2:
            return (6, mass, pid)
        # Omega (S=-3)
        elif B == 1 and S == -3:
            return (7, mass, pid)
        else:
            return (8, mass, pid)

    # Actually, let's preserve the PDG2020 ordering style instead
    # The original file groups: light mesons, strange mesons, N, Delta, Lambda, Sigma, Xi, Omega
    # Let's keep this ordering and just insert new states at appropriate positions

    # Build the output: keep original order for existing particles, insert new ones
    # First, let's separate by sector

    light_mesons = []
    strange_mesons = []
    baryons_other = []

    for p in particles_2025:
        pid = p['pdgid']
        B = p['B']
        S = p['S']

        if B == 0 and S == 0 and p['C'] == 0:
            # Light unflavored meson (including hidden strangeness like eta, phi)
            light_mesons.append(p)
        elif B == 0 and S != 0:
            # Strange meson (open strangeness)
            strange_mesons.append(p)
        else:
            baryons_other.append(p)

    # Also need to handle mesons with hidden strangeness (absS > 0 but S=0)
    # like eta, eta', phi etc. These are in the light unflavored section already
    # Actually looking at the original file, particles like eta (absS=1.33),
    # eta'(958) (absS=0.667), phi(1020) (absS=2) etc are listed in the
    # light unflavored section. The file ordering is by PDG ID groups.

    # Let me just keep the original ordering and sort new particles into position
    # We know the original order: light mesons (up to f(2)(2340)), then kaons, then baryons

    # Sort each sector by mass
    light_mesons.sort(key=lambda p: (p['mass'], p['pdgid']))
    strange_mesons.sort(key=lambda p: (p['mass'], p['pdgid']))
    # Keep baryon order from original file

    all_particles = light_mesons + strange_mesons + baryons_other

    # =====================================================================
    # Compute thresholds for new particles
    # =====================================================================
    particles_by_id = {p['pdgid']: p for p in all_particles}
    for p in all_particles:
        pid = p['pdgid']
        if p['threshold'] == 0 and p['stable'] == 0:
            if pid in NEW_DECAYS:
                thresh = compute_threshold(pid, NEW_DECAYS[pid], particles_by_id)
                p['threshold'] = thresh
                print(f"  Computed threshold for {p['name']}: {thresh:.6f}")

    # =====================================================================
    # Write list.dat
    # =====================================================================

    with open(OUTPUT_LIST, 'w') as f:
        f.write("# PDG2025 list containing strange and non-strange hadrons only\n")
        f.write("#         pdgid                name         stable      mass[GeV]     degeneracy     statistics              B              Q              S              C            |S|            |C|     width[GeV] threshold[GeV]\n")

        for p in all_particles:
            # Format each field to match the original column alignment
            line = f"{p['pdgid']:>15} {p['name']:>20} {p['stable']:>13} {p['mass']:>14g} {p['degeneracy']:>14} {p['statistics']:>14} {p['B']:>14} {p['Q']:>14} {p['S']:>14} {p['C']:>14} {p['absS']:>14g} {p['absC']:>14g} {p['width']:>14g} {p['threshold']:>14g}"
            f.write(line + "\n")

    print(f"Written {len(all_particles)} particles to {OUTPUT_LIST}")

    # =====================================================================
    # Build updated decays
    # =====================================================================

    decays_2025 = OrderedDict()

    # Copy decays from PDG2020 for particles that are kept
    for pid, channels in decays_2020.items():
        if pid in kept_ids:
            # Check if any daughter in any channel is a removed particle
            updated_channels = []
            for ch in channels:
                # Check if all daughters exist
                all_daughters_exist = True
                for d in ch['daughters']:
                    absd = abs(d)
                    # Special daughter IDs: 22 (gamma), 11/-11 (e), 12/-12 (nu_e),
                    # 13/-13 (mu), 14/-14 (nu_mu), 310 (K_S), 130 (K_L)
                    if absd in [22, 11, 12, 13, 14, 310, 130]:
                        continue
                    if absd in REMOVE_IDS:
                        all_daughters_exist = False
                        break
                if all_daughters_exist:
                    updated_channels.append(ch)

            if updated_channels:
                # Renormalize branching ratios
                total_br = sum(ch['br'] for ch in updated_channels)
                if total_br > 0 and abs(total_br - 1.0) > 1e-6:
                    for ch in updated_channels:
                        ch['br'] = ch['br'] / total_br
                decays_2025[pid] = updated_channels
            elif channels:
                # All channels were removed - shouldn't happen for remaining particles
                print(f"  WARNING: All decay channels removed for {pid}")
                decays_2025[pid] = channels  # Keep original

    # Add new decays
    for pid, channels in NEW_DECAYS.items():
        decays_2025[pid] = channels

    # Apply decay updates for existing particles
    n_decay_updates = 0
    for pid, channels in DECAY_UPDATES.items():
        if pid in decays_2025 or pid in kept_ids:
            decays_2025[pid] = channels
            n_decay_updates += 1
            # Find particle name for logging
            pname = str(pid)
            for p in all_particles:
                if p['pdgid'] == pid:
                    pname = p['name']
                    break
            print(f"  Updated decays for {pname} (ID {pid}): {len(channels)} channels")

    # =====================================================================
    # Write decays.dat
    # =====================================================================

    with open(OUTPUT_DECAYS, 'w') as f:
        f.write("# the list of decays\n")
        f.write("# each entry consists of the following:\n")
        f.write("# a line with the pdgid of decaying particle\n")
        f.write("# a line with the number of decay channels\n")
        f.write("# for each channel a line containing whitespace-separated values of the channel branching ratio and pdg ids of the daughter products\n")
        f.write("# everything after the # symbol is treated as a comment and ignored\n")
        f.write("# decays of antiparticles are not listed but generated from the listed decays of particles\n")
        f.write("\n")

        # Write decays in particle list order
        pid_order = [p['pdgid'] for p in all_particles]
        written_pids = set()

        for pid in pid_order:
            if pid in decays_2025 and pid not in written_pids:
                channels = decays_2025[pid]
                # Find particle name
                name = ""
                for p in all_particles:
                    if p['pdgid'] == pid:
                        name = p['name']
                        break

                f.write(f"{pid:<20}                     # {name}\n")
                f.write(f"{len(channels):<20}                     # {len(channels)} decay channel{'s' if len(channels) > 1 else ''}\n")
                for ch in channels:
                    daughters_str = " ".join(str(d) for d in ch['daughters'])
                    comment = ch.get('comment', '')
                    # Format branching ratio with limited precision
                    br_val = ch['br']
                    if abs(br_val - 1.0) < 1e-10:
                        br_str = "1.0"
                    elif abs(br_val - round(br_val, 6)) < 1e-10:
                        br_str = f"{br_val:.6g}"
                    else:
                        br_str = f"{br_val:.6g}"
                    f.write(f"{br_str:<16}{daughters_str:<20} {comment}\n")
                f.write("\n")
                written_pids.add(pid)

        # Also write extra decay entries (charm, nuclei) from PDG2020
        # Apply CHARM_DECAY_UPDATES where available
        extra_count = 0
        charm_updated_count = 0
        for pid, channels in decays_2020.items():
            if pid not in written_pids and pid not in REMOVE_IDS:
                # Check if this charm particle has an update
                pname = pid_name_map.get(pid, str(pid))
                if pid in CHARM_DECAY_UPDATES:
                    channels = CHARM_DECAY_UPDATES[pid]
                    charm_updated_count += 1
                    label = f"{pname} (PDG2025 updated)"
                else:
                    label = f"{pname} (from PDG2020)"
                f.write(f"{pid:<20}                     # {label}\n")
                f.write(f"{len(channels):<20}                     # {len(channels)} decay channel{'s' if len(channels) > 1 else ''}\n")
                for ch in channels:
                    daughters_str = " ".join(str(d) for d in ch['daughters'])
                    comment = ch.get('comment', '')
                    br_val = ch['br']
                    if abs(br_val - 1.0) < 1e-10:
                        br_str = "1.0"
                    elif abs(br_val - round(br_val, 6)) < 1e-10:
                        br_str = f"{br_val:.6g}"
                    else:
                        br_str = f"{br_val:.6g}"
                    f.write(f"{br_str:<16}{daughters_str:<20} {comment}\n")
                f.write("\n")
                written_pids.add(pid)
                extra_count += 1

        # Write new charm decay entries not in PDG2020
        for pid, channels in CHARM_DECAY_UPDATES.items():
            if pid not in written_pids:
                pname = pid_name_map.get(pid, str(pid))
                f.write(f"{pid:<20}                     # {pname} (PDG2025 new)\n")
                f.write(f"{len(channels):<20}                     # {len(channels)} decay channel{'s' if len(channels) > 1 else ''}\n")
                for ch in channels:
                    daughters_str = " ".join(str(d) for d in ch['daughters'])
                    comment = ch.get('comment', '')
                    br_val = ch['br']
                    if abs(br_val - 1.0) < 1e-10:
                        br_str = "1.0"
                    elif abs(br_val - round(br_val, 6)) < 1e-10:
                        br_str = f"{br_val:.6g}"
                    else:
                        br_str = f"{br_val:.6g}"
                    f.write(f"{br_str:<16}{daughters_str:<20} {comment}\n")
                f.write("\n")
                written_pids.add(pid)
                extra_count += 1
                charm_updated_count += 1

    print(f"Written {len(decays_2025)} hadron decay entries + {extra_count} extra entries ({charm_updated_count} charm updated/new) to {OUTPUT_DECAYS}")

    # =====================================================================
    # Summary of changes
    # =====================================================================

    print("\n" + "="*70)
    print("SUMMARY OF CHANGES FROM PDG2020 TO PDG2025")
    print("="*70)

    print(f"\nTotal particles: {len(particles_2020)} (PDG2020) -> {len(all_particles)} (PDG2025)")

    print(f"\nRemoved {len(removed_names)} particles:")
    for name in removed_names:
        print(name)

    print(f"\nAdded {len(NEW_MESONS) + len(NEW_STRANGE_MESONS)} new particles:")
    for p in NEW_MESONS + NEW_STRANGE_MESONS:
        print(f"  Added: {p['name']} (ID {p['pdgid']}, M={p['mass']:.3f} GeV, W={p['width']:.3f} GeV)")

    # Count mass/width updates
    n_mass_updates = 0
    n_width_updates = 0
    for p_new in all_particles:
        pid = p_new['pdgid']
        for p_old in particles_2020:
            if p_old['pdgid'] == pid:
                if abs(p_new['mass'] - p_old['mass']) > 1e-6:
                    n_mass_updates += 1
                if abs(p_new['width'] - p_old['width']) > 1e-10 and p_old['width'] > 0:
                    n_width_updates += 1
                break

    print(f"\nMass updates: {n_mass_updates} particles")
    print(f"Width updates: {n_width_updates} particles")
    print(f"Decay updates: {n_decay_updates} particles (branching fractions updated from PDG 2025)")

if __name__ == '__main__':
    generate_pdg2025()
