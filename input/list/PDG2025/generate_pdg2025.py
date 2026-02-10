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
BASE_DIR = os.path.dirname(SCRIPT_DIR)  # input/list/
PDG2020_LIST = os.path.join(BASE_DIR, "PDG2020", "list.dat")
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
    # PDG 2025: pi+pi-pi0 (89.2%), pi0 gamma (8.33%), pi+pi- (1.53%)
    # Remaining ~0.9% is eta gamma + other rare modes - distribute among main channels
    223: [
        {'br': 0.892, 'daughters': [211, -211, 111], 'comment': '# omega(782) -> pi+ + pi- + pi0'},
        {'br': 0.0833, 'daughters': [111, 22], 'comment': '# omega(782) -> pi0 + gamma'},
        {'br': 0.0153, 'daughters': [211, -211], 'comment': '# omega(782) -> pi+ + pi-'},
        {'br': 0.0045, 'daughters': [221, 22], 'comment': '# omega(782) -> eta + gamma'},
        {'br': 0.0049, 'daughters': [111, -11, 11], 'comment': '# omega(782) -> pi0 + e+ + e-'},
    ],

    # eta'(958): Updated branching fractions from PDG 2025
    # PDG 2025: pi+pi-eta (42.5%), rho0 gamma (29.48%), pi0pi0eta (22.4%), omega gamma (2.52%), gamma gamma (2.307%)
    # Very minor changes vs PDG2020
    331: [
        {'br': 0.425, 'daughters': [211, -211, 221], 'comment': "# eta'(958) -> pi+ + pi- + eta"},
        {'br': 0.2948, 'daughters': [113, 22], 'comment': "# eta'(958) -> rho(770)0 + gamma"},
        {'br': 0.224, 'daughters': [111, 111, 221], 'comment': "# eta'(958) -> pi0 + pi0 + eta"},
        {'br': 0.0252, 'daughters': [223, 22], 'comment': "# eta'(958) -> omega(782) + gamma"},
        {'br': 0.023, 'daughters': [22, 22], 'comment': "# eta'(958) -> gamma + gamma"},
        {'br': 0.008, 'daughters': [211, -211, 111], 'comment': "# eta'(958) -> pi+ + pi- + pi0"},
    ],

    # phi(1020): Updated branching fractions from PDG 2025
    # PDG 2025: K+K- (49.9%), KL KS (33.6%), rho pi + pi+pi-pi0 (14.9%), eta gamma (1.306%)
    333: [
        {'br': 0.492, 'daughters': [321, -321], 'comment': '# phi(1020) -> K+ + K-'},
        {'br': 0.340, 'daughters': [311, -311], 'comment': '# phi(1020) -> K0 + anti-K0'},
        {'br': 0.153, 'daughters': [211, -211, 111], 'comment': '# phi(1020) -> pi+ + pi- + pi0'},
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

    # Delta(1600)-: PDG 2025 ranges: Npi (8-24%), Delta(1232)pi (58-82%), N(1440)pi (17-27%)
    # Current: Npi 15%, Deltapi 75%, N(1440)pi 10% -> N(1440)pi too low
    # Updated: Npi 15%, Deltapi 65%, N(1440)pi 20%
    # Delta(1600) is P33 (I=3/2)
    # Isospin C-G for I=3/2(Q=-1) -> I=1/2 N + I=1 pi: only n+pi-
    # Isospin C-G for I=3/2(Q=-1) -> I=3/2 Delta + I=1 pi: Delta0 pi- (2/5), Delta- pi0 (3/5)
    31114: [  # Delta(1600)-
        {'br': 0.15, 'daughters': [2112, -211], 'comment': '# Delta(1600)- -> n + pi-'},
        {'br': 0.26, 'daughters': [2114, -211], 'comment': '# Delta(1600)- -> Delta(1232)0 + pi-'},
        {'br': 0.39, 'daughters': [1114, 111], 'comment': '# Delta(1600)- -> Delta(1232)- + pi0'},
        {'br': 0.20, 'daughters': [12112, -211], 'comment': '# Delta(1600)- -> N(1440)0 + pi-'},
    ],
    32114: [  # Delta(1600)0
        {'br': 0.10, 'daughters': [2112, 111], 'comment': '# Delta(1600)0 -> n + pi0'},
        {'br': 0.05, 'daughters': [2212, -211], 'comment': '# Delta(1600)0 -> p + pi-'},
        {'br': 0.347, 'daughters': [2214, -211], 'comment': '# Delta(1600)0 -> Delta(1232)+ + pi-'},
        {'br': 0.043, 'daughters': [2114, 111], 'comment': '# Delta(1600)0 -> Delta(1232)0 + pi0'},
        {'br': 0.26, 'daughters': [1114, 211], 'comment': '# Delta(1600)0 -> Delta(1232)- + pi+'},
        {'br': 0.133, 'daughters': [12112, 111], 'comment': '# Delta(1600)0 -> N(1440)0 + pi0'},
        {'br': 0.067, 'daughters': [12212, -211], 'comment': '# Delta(1600)0 -> N(1440)+ + pi-'},
    ],
    32214: [  # Delta(1600)+
        {'br': 0.05, 'daughters': [2112, 211], 'comment': '# Delta(1600)+ -> n + pi+'},
        {'br': 0.10, 'daughters': [2212, 111], 'comment': '# Delta(1600)+ -> p + pi0'},
        {'br': 0.26, 'daughters': [2224, -211], 'comment': '# Delta(1600)+ -> Delta(1232)++ + pi-'},
        {'br': 0.043, 'daughters': [2214, 111], 'comment': '# Delta(1600)+ -> Delta(1232)+ + pi0'},
        {'br': 0.347, 'daughters': [2114, 211], 'comment': '# Delta(1600)+ -> Delta(1232)0 + pi+'},
        {'br': 0.067, 'daughters': [12112, 211], 'comment': '# Delta(1600)+ -> N(1440)0 + pi+'},
        {'br': 0.133, 'daughters': [12212, 111], 'comment': '# Delta(1600)+ -> N(1440)+ + pi0'},
    ],
    32224: [  # Delta(1600)++
        {'br': 0.15, 'daughters': [2212, 211], 'comment': '# Delta(1600)++ -> p + pi+'},
        {'br': 0.39, 'daughters': [2224, 111], 'comment': '# Delta(1600)++ -> Delta(1232)++ + pi0'},
        {'br': 0.26, 'daughters': [2214, 211], 'comment': '# Delta(1600)++ -> Delta(1232)+ + pi+'},
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
        # These are not in the particle list but used by the event generator
        extra_count = 0
        for pid, channels in decays_2020.items():
            if pid not in written_pids and pid not in REMOVE_IDS:
                # This is an extra entry (charm particle or nucleus)
                f.write(f"{pid:<20}                     # extra (charm/nuclear)\n")
                f.write(f"{len(channels):<20}                     # {len(channels)} decay channel{'s' if len(channels) > 1 else ''}\n")
                for ch in channels:
                    daughters_str = " ".join(str(d) for d in ch['daughters'])
                    comment = ch.get('comment', '')
                    br_val = ch['br']
                    if abs(br_val - 1.0) < 1e-10:
                        br_str = "1.0"
                    else:
                        br_str = f"{br_val:.6g}"
                    f.write(f"{br_str:<16}{daughters_str:<20} {comment}\n")
                f.write("\n")
                written_pids.add(pid)
                extra_count += 1

    print(f"Written {len(decays_2025)} hadron decay entries + {extra_count} extra (charm/nuclear) entries to {OUTPUT_DECAYS}")

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
