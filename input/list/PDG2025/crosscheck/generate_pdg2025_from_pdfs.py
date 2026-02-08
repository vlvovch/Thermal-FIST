#!/usr/bin/env python3
"""
Generate PDG2025 particle list and decay list for Thermal-FIST
Independent generation from PDG 2025 Summary Table PDFs
All particle data hardcoded from manual reading of the PDFs
"""

import os
from collections import OrderedDict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_LIST = os.path.join(SCRIPT_DIR, "list_from_pdfs.dat")
OUTPUT_DECAYS = os.path.join(SCRIPT_DIR, "decays_from_pdfs.dat")

# =====================================================================
# PARTICLE DATA: Hardcoded from PDG 2025 Summary Table PDFs
# Each entry: pdgid -> dict with properties
# Mass and width in GeV
# =====================================================================

PARTICLES = OrderedDict()

# =====================================================================
# LIGHT UNFLAVORED MESONS (from rpp2025-tab-mesons-light.pdf)
# =====================================================================

# pi0: I^G(J^PC) = 1-(0-+), mass = 134.9768 MeV, width = 7.81 eV
PARTICLES[111] = {
    'name': 'pi0', 'mass': 0.134977, 'width': 7.81e-09,
    'J': 0, 'I': 1, 'G': -1, 'P': -1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 1, 'statistics': -1,
}

# pi+: I^G(J^PC) = 1-(0-+), mass = 139.570 MeV, width = 2.5284e-8 eV (tau = 2.6e-8 s)
PARTICLES[211] = {
    'name': 'pi+', 'mass': 0.13957, 'width': 2.5284e-17,
    'J': 0, 'I': 1, 'G': -1, 'P': -1, 'C': 1,
    'B': 0, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 1, 'statistics': -1,
}

# eta: I^G(J^PC) = 0+(0-+), mass = 547.862 MeV, width = 1.31 keV
PARTICLES[221] = {
    'name': 'eta', 'mass': 0.547862, 'width': 1.31e-06,
    'J': 0, 'I': 0, 'G': 1, 'P': -1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 1.33333, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# f0(500)/sigma: I^G(J^PC) = 0+(0++), mass ~ 400-800 MeV (use 600), width ~ 400-800 (use 450)
PARTICLES[9000221] = {
    'name': 'f(0)(500)', 'mass': 0.6, 'width': 0.45,
    'J': 0, 'I': 0, 'G': 1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1, 'degeneracy_override': 0,
}

# rho(770)+: I^G(J^PC) = 1+(1--), mass = 775.11 MeV, width = 149.1 MeV
PARTICLES[213] = {
    'name': 'rho(770)+', 'mass': 0.77511, 'width': 0.1491,
    'J': 1, 'I': 1, 'G': 1, 'P': -1, 'C': -1,
    'B': 0, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# rho(770)0: I^G(J^PC) = 1+(1--), mass = 775.26 MeV, width = 147.4 MeV
PARTICLES[113] = {
    'name': 'rho(770)0', 'mass': 0.77526, 'width': 0.1474,
    'J': 1, 'I': 1, 'G': 1, 'P': -1, 'C': -1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# omega(782): I^G(J^PC) = 0-(1--), mass = 782.66 MeV, width = 8.68 MeV
PARTICLES[223] = {
    'name': 'omega(782)', 'mass': 0.78266, 'width': 0.00868,
    'J': 1, 'I': 0, 'G': -1, 'P': -1, 'C': -1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# eta'(958): I^G(J^PC) = 0+(0-+), mass = 957.78 MeV, width = 0.188 MeV
PARTICLES[331] = {
    'name': "eta'(958)", 'mass': 0.95778, 'width': 0.000188,
    'J': 0, 'I': 0, 'G': 1, 'P': -1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0.66667, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# a0(980)0: I^G(J^PC) = 1-(0++), mass = 980 MeV, width = 50-100 MeV (use 75)
PARTICLES[9000111] = {
    'name': 'a(0)(980)0', 'mass': 0.98, 'width': 0.075,
    'J': 0, 'I': 1, 'G': -1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# a0(980)+: same as neutral
PARTICLES[9000211] = {
    'name': 'a(0)(980)+', 'mass': 0.98, 'width': 0.075,
    'J': 0, 'I': 1, 'G': -1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# f0(980): I^G(J^PC) = 0+(0++), mass = 990 MeV, width = 10-100 MeV (use 60)
PARTICLES[9010221] = {
    'name': 'f(0)(980)', 'mass': 0.99, 'width': 0.06,
    'J': 0, 'I': 0, 'G': 1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# phi(1020): I^G(J^PC) = 0-(1--), mass = 1019.461 MeV, width = 4.249 MeV
PARTICLES[333] = {
    'name': 'phi(1020)', 'mass': 1.01946, 'width': 0.004249,
    'J': 1, 'I': 0, 'G': -1, 'P': -1, 'C': -1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 2, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# h1(1170): I^G(J^PC) = 0-(1+-), mass = 1166 MeV, width = 375 MeV
PARTICLES[10223] = {
    'name': 'h(1)(1170)', 'mass': 1.166, 'width': 0.375,
    'J': 1, 'I': 0, 'G': -1, 'P': 1, 'C': -1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# b1(1235)0: I^G(J^PC) = 1+(1+-), mass = 1229.5 MeV, width = 142 MeV
PARTICLES[10113] = {
    'name': 'b(1)(1235)0', 'mass': 1.2295, 'width': 0.142,
    'J': 1, 'I': 1, 'G': 1, 'P': 1, 'C': -1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[10213] = {
    'name': 'b(1)(1235)+', 'mass': 1.2295, 'width': 0.142,
    'J': 1, 'I': 1, 'G': 1, 'P': 1, 'C': -1,
    'B': 0, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# a1(1260)0: I^G(J^PC) = 1-(1++), mass = 1230 MeV, width = 420 MeV
PARTICLES[20113] = {
    'name': 'a(1)(1260)0', 'mass': 1.23, 'width': 0.42,
    'J': 1, 'I': 1, 'G': -1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[20213] = {
    'name': 'a(1)(1260)+', 'mass': 1.23, 'width': 0.42,
    'J': 1, 'I': 1, 'G': -1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# f2(1270): I^G(J^PC) = 0+(2++), mass = 1275.4 MeV, width = 186.6 MeV
PARTICLES[225] = {
    'name': 'f(2)(1270)', 'mass': 1.2754, 'width': 0.1866,
    'J': 2, 'I': 0, 'G': 1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# f1(1285): I^G(J^PC) = 0+(1++), mass = 1281.8 MeV, width = 23.0 MeV
PARTICLES[20223] = {
    'name': 'f(1)(1285)', 'mass': 1.2818, 'width': 0.023,
    'J': 1, 'I': 0, 'G': 1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# eta(1295): I^G(J^PC) = 0+(0-+), mass = 1294 MeV, width = 55 MeV
PARTICLES[100221] = {
    'name': 'eta(1295)', 'mass': 1.294, 'width': 0.055,
    'J': 0, 'I': 0, 'G': 1, 'P': -1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# pi(1300)0: I^G(J^PC) = 1-(0-+), mass = 1300 MeV, width = 200-600 (use 400) MeV
PARTICLES[100111] = {
    'name': 'pi(1300)0', 'mass': 1.3, 'width': 0.4,
    'J': 0, 'I': 1, 'G': -1, 'P': -1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[100211] = {
    'name': 'pi(1300)+', 'mass': 1.3, 'width': 0.4,
    'J': 0, 'I': 1, 'G': -1, 'P': -1, 'C': 1,
    'B': 0, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# a2(1320)0: I^G(J^PC) = 1-(2++), mass = 1318.2 MeV, width = 107 MeV
PARTICLES[115] = {
    'name': 'a(2)(1320)0', 'mass': 1.3182, 'width': 0.107,
    'J': 2, 'I': 1, 'G': -1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[215] = {
    'name': 'a(2)(1320)+', 'mass': 1.3182, 'width': 0.107,
    'J': 2, 'I': 1, 'G': -1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# f0(1370): I^G(J^PC) = 0+(0++), mass = 1200-1500 (use 1350), width = 200-500 (use 350)
PARTICLES[10221] = {
    'name': 'f(0)(1370)', 'mass': 1.35, 'width': 0.35,
    'J': 0, 'I': 0, 'G': 1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# eta(1405): I^G(J^PC) = 0+(0-+), mass = 1408.7 MeV, width = 50.3 MeV
PARTICLES[9020221] = {
    'name': 'eta(1405)', 'mass': 1.4087, 'width': 0.0503,
    'J': 0, 'I': 0, 'G': 1, 'P': -1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0.8, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# h1(1415): I^G(J^PC) = 0-(1+-), mass = 1409 MeV, width = 78 MeV (Note: called h1(1415) in PDG, contains ss-bar)
PARTICLES[10333] = {
    'name': 'h(1)(1415)', 'mass': 1.409, 'width': 0.078,
    'J': 1, 'I': 0, 'G': -1, 'P': 1, 'C': -1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# omega(1420): I^G(J^PC) = 0-(1--), mass = 1410 MeV, width = 290 MeV
PARTICLES[100223] = {
    'name': 'omega(1420)', 'mass': 1.41, 'width': 0.29,
    'J': 1, 'I': 0, 'G': -1, 'P': -1, 'C': -1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# f1(1420): I^G(J^PC) = 0+(1++), mass = 1428.4 MeV, width = 56.7 MeV
PARTICLES[20333] = {
    'name': 'f(1)(1420)', 'mass': 1.4284, 'width': 0.0567,
    'J': 1, 'I': 0, 'G': 1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 2, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# a0(1450)0: I^G(J^PC) = 1-(0++), mass = 1439 MeV, width = 258 MeV
PARTICLES[10111] = {
    'name': 'a(0)(1450)0', 'mass': 1.439, 'width': 0.258,
    'J': 0, 'I': 1, 'G': -1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[10211] = {
    'name': 'a(0)(1450)+', 'mass': 1.439, 'width': 0.258,
    'J': 0, 'I': 1, 'G': -1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# rho(1450)0: I^G(J^PC) = 1+(1--), mass = 1465 MeV, width = 400 MeV
PARTICLES[100113] = {
    'name': 'rho(1450)0', 'mass': 1.465, 'width': 0.4,
    'J': 1, 'I': 1, 'G': 1, 'P': -1, 'C': -1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[100213] = {
    'name': 'rho(1450)+', 'mass': 1.465, 'width': 0.4,
    'J': 1, 'I': 1, 'G': 1, 'P': -1, 'C': -1,
    'B': 0, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# eta(1475) [called eta(1475) in PDG, listed under eta(1405) section in some editions]
# I^G(J^PC) = 0+(0-+), mass = 1476 MeV, width = 96 MeV
PARTICLES[100331] = {
    'name': 'eta(1475)', 'mass': 1.476, 'width': 0.096,
    'J': 0, 'I': 0, 'G': 1, 'P': -1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 2, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# f2'(1525): I^G(J^PC) = 0+(2++), mass = 1517.3 MeV, width = 72 MeV
PARTICLES[335] = {
    'name': "f(2)'(1525)", 'mass': 1.5173, 'width': 0.072,
    'J': 2, 'I': 0, 'G': 1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 2, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# f0(1500): I^G(J^PC) = 0+(0++), mass = 1522 MeV, width = 108 MeV
PARTICLES[9030221] = {
    'name': 'f(0)(1500)', 'mass': 1.522, 'width': 0.108,
    'J': 0, 'I': 0, 'G': 1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# f2(1565): I^G(J^PC) = 0+(2++), mass = 1571 MeV, width = 133 MeV  [NEW in PDG2025]
PARTICLES[9010225] = {
    'name': 'f(2)(1565)', 'mass': 1.571, 'width': 0.133,
    'J': 2, 'I': 0, 'G': 1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# eta2(1645): I^G(J^PC) = 0+(2-+), mass = 1617 MeV, width = 181 MeV
PARTICLES[10225] = {
    'name': 'eta(2)(1645)', 'mass': 1.617, 'width': 0.181,
    'J': 2, 'I': 0, 'G': 1, 'P': -1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# pi1(1600)0: I^G(J^PC) = 1-(1-+), mass = 1645 MeV (Breit-Wigner), width = 370 MeV
PARTICLES[9010113] = {
    'name': 'pi(1)(1600)0', 'mass': 1.645, 'width': 0.37,
    'J': 1, 'I': 1, 'G': -1, 'P': -1, 'C': -1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[9010213] = {
    'name': 'pi(1)(1600)+', 'mass': 1.645, 'width': 0.37,
    'J': 1, 'I': 1, 'G': -1, 'P': -1, 'C': -1,
    'B': 0, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# a1(1640)0: I^G(J^PC) = 1-(1++), mass = 1655 MeV, width = 250 MeV
PARTICLES[9020113] = {
    'name': 'a(1)(1640)0', 'mass': 1.655, 'width': 0.25,
    'J': 1, 'I': 1, 'G': -1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[9020213] = {
    'name': 'a(1)(1640)+', 'mass': 1.655, 'width': 0.25,
    'J': 1, 'I': 1, 'G': -1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# omega3(1670): I^G(J^PC) = 0-(3--), mass = 1667 MeV, width = 168 MeV
PARTICLES[227] = {
    'name': 'omega(3)(1670)', 'mass': 1.667, 'width': 0.168,
    'J': 3, 'I': 0, 'G': -1, 'P': -1, 'C': -1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# omega(1650): I^G(J^PC) = 0-(1--), mass = 1670 MeV, width = 315 MeV
PARTICLES[30223] = {
    'name': 'omega(1650)', 'mass': 1.67, 'width': 0.315,
    'J': 1, 'I': 0, 'G': -1, 'P': -1, 'C': -1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# pi2(1670)0: I^G(J^PC) = 1-(2-+), mass = 1670.6 MeV, width = 258 MeV
PARTICLES[10115] = {
    'name': 'pi(2)(1670)0', 'mass': 1.6706, 'width': 0.258,
    'J': 2, 'I': 1, 'G': -1, 'P': -1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[10215] = {
    'name': 'pi(2)(1670)+', 'mass': 1.6706, 'width': 0.258,
    'J': 2, 'I': 1, 'G': -1, 'P': -1, 'C': 1,
    'B': 0, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# phi(1680): I^G(J^PC) = 0-(1--), mass = 1680 MeV, width = 150 MeV
PARTICLES[100333] = {
    'name': 'phi(1680)', 'mass': 1.68, 'width': 0.15,
    'J': 1, 'I': 0, 'G': -1, 'P': -1, 'C': -1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 2, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# rho3(1690)0: I^G(J^PC) = 1+(3--), mass = 1688.8 MeV, width = 161 MeV
PARTICLES[117] = {
    'name': 'rho(3)(1690)0', 'mass': 1.6888, 'width': 0.161,
    'J': 3, 'I': 1, 'G': 1, 'P': -1, 'C': -1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[217] = {
    'name': 'rho(3)(1690)+', 'mass': 1.6888, 'width': 0.161,
    'J': 3, 'I': 1, 'G': 1, 'P': -1, 'C': -1,
    'B': 0, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# a2(1700)0: I^G(J^PC) = 1-(2++), mass = 1706 MeV, width = 380 MeV
PARTICLES[9000115] = {
    'name': 'a(2)(1700)0', 'mass': 1.706, 'width': 0.38,
    'J': 2, 'I': 1, 'G': -1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[9000215] = {
    'name': 'a(2)(1700)+', 'mass': 1.706, 'width': 0.38,
    'J': 2, 'I': 1, 'G': -1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# rho(1700)0: I^G(J^PC) = 1+(1--), mass = 1720 MeV, width = 250 MeV
PARTICLES[30113] = {
    'name': 'rho(1700)0', 'mass': 1.72, 'width': 0.25,
    'J': 1, 'I': 1, 'G': 1, 'P': -1, 'C': -1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[30213] = {
    'name': 'rho(1700)+', 'mass': 1.72, 'width': 0.25,
    'J': 1, 'I': 1, 'G': 1, 'P': -1, 'C': -1,
    'B': 0, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# f0(1710): I^G(J^PC) = 0+(0++), mass = 1733 MeV, width = 150 MeV
PARTICLES[10331] = {
    'name': 'f(0)(1710)', 'mass': 1.733, 'width': 0.15,
    'J': 0, 'I': 0, 'G': 1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# pi(1800)0: I^G(J^PC) = 1-(0-+), mass = 1810 MeV, width = 215 MeV
PARTICLES[9010111] = {
    'name': 'pi(1800)0', 'mass': 1.81, 'width': 0.215,
    'J': 0, 'I': 1, 'G': -1, 'P': -1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[9010211] = {
    'name': 'pi(1800)+', 'mass': 1.81, 'width': 0.215,
    'J': 0, 'I': 1, 'G': -1, 'P': -1, 'C': 1,
    'B': 0, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# eta2(1870): I^G(J^PC) = 0+(2-+), mass = 1842 MeV, width = 225 MeV
PARTICLES[10335] = {
    'name': 'eta2(1870)', 'mass': 1.842, 'width': 0.225,
    'J': 2, 'I': 0, 'G': 1, 'P': -1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# phi3(1850): I^G(J^PC) = 0-(3--), mass = 1854 MeV, width = 87 MeV
PARTICLES[337] = {
    'name': 'phi(3)(1850)', 'mass': 1.854, 'width': 0.087,
    'J': 3, 'I': 0, 'G': -1, 'P': -1, 'C': -1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 2, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# pi2(1880)0: I^G(J^PC) = 1-(2-+), mass = 1895 MeV, width = 237 MeV
PARTICLES[20115] = {
    'name': 'pi2(1880)0', 'mass': 1.895, 'width': 0.237,
    'J': 2, 'I': 1, 'G': -1, 'P': -1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[20215] = {
    'name': 'pi2(1880)+', 'mass': 1.895, 'width': 0.237,
    'J': 2, 'I': 1, 'G': -1, 'P': -1, 'C': 1,
    'B': 0, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# f2(1950): I^G(J^PC) = 0+(2++), mass = 1936 MeV, width = 464 MeV
PARTICLES[9050225] = {
    'name': 'f(2)(1950)', 'mass': 1.936, 'width': 0.464,
    'J': 2, 'I': 0, 'G': 1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# a4(1970)0: I^G(J^PC) = 1-(4++), mass = 1967 MeV, width = 324 MeV
PARTICLES[119] = {
    'name': 'a(4)(1970)0', 'mass': 1.967, 'width': 0.324,
    'J': 4, 'I': 1, 'G': -1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[219] = {
    'name': 'a(4)(1970)+', 'mass': 1.967, 'width': 0.324,
    'J': 4, 'I': 1, 'G': -1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# f0(2020): I^G(J^PC) = 0+(0++), mass = 1982 MeV, width = 440 MeV [NEW in PDG2025]
PARTICLES[9050221] = {
    'name': 'f(0)(2020)', 'mass': 1.982, 'width': 0.44,
    'J': 0, 'I': 0, 'G': 1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# f2(2010): I^G(J^PC) = 0+(2++), mass = 2010 MeV, width = 200 MeV (Note: called f2(2010) in PDG listings)
PARTICLES[9060225] = {
    'name': 'f(2)(2010)', 'mass': 2.01, 'width': 0.2,
    'J': 2, 'I': 0, 'G': 1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 4, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# f4(2050): I^G(J^PC) = 0+(4++), mass = 2018 MeV, width = 237 MeV
PARTICLES[229] = {
    'name': 'f(4)(2050)', 'mass': 2.018, 'width': 0.237,
    'J': 4, 'I': 0, 'G': 1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# f2(2150): I^G(J^PC) = 0+(2++), mass = 2157 MeV, width = 152 MeV [NEW in PDG2025]
PARTICLES[9070225] = {
    'name': 'f(2)(2150)', 'mass': 2.157, 'width': 0.152,
    'J': 2, 'I': 0, 'G': 1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# phi(2170): I^G(J^PC) = 0-(1--), mass = 2164 MeV, width = 88 MeV
PARTICLES[200333] = {
    'name': 'phi(2170)', 'mass': 2.164, 'width': 0.088,
    'J': 1, 'I': 0, 'G': -1, 'P': -1, 'C': -1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 2, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# f2(2300): I^G(J^PC) = 0+(2++), mass = 2297 MeV, width = 150 MeV
PARTICLES[9080225] = {
    'name': 'f(2)(2300)', 'mass': 2.297, 'width': 0.15,
    'J': 2, 'I': 0, 'G': 1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 4, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# f2(2340): I^G(J^PC) = 0+(2++), mass = 2346 MeV, width = 331 MeV
PARTICLES[9090225] = {
    'name': 'f(2)(2340)', 'mass': 2.346, 'width': 0.331,
    'J': 2, 'I': 0, 'G': 1, 'P': 1, 'C': 1,
    'B': 0, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 4, 'absC': 0,
    'stable': 0, 'statistics': -1,
}


# =====================================================================
# STRANGE MESONS (from rpp2025-tab-mesons-strange.pdf)
# =====================================================================

# K+: I(J^P) = 1/2(0-), mass = 493.677 MeV, width (tau = 1.238e-8 s)
PARTICLES[321] = {
    'name': 'K+', 'mass': 0.493677, 'width': 5.317e-17,
    'J': 0, 'I': 0.5, 'P': -1,
    'B': 0, 'Q': 1, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 1, 'statistics': -1,
}

# K0: I(J^P) = 1/2(0-), mass = 497.611 MeV
PARTICLES[311] = {
    'name': 'K0', 'mass': 0.497611, 'width': 0,
    'J': 0, 'I': 0.5, 'P': -1,
    'B': 0, 'Q': 0, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 1, 'statistics': -1,
}

# K0*(700)0 (kappa): I(J^P) = 1/2(0+), BW mass = 838 MeV, width = 463 MeV [NEW in PDG2025]
# degeneracy_override=0: cancellation in I=3/2 Kpi repulsive channel, analogous to f(0)(500)/sigma
PARTICLES[9000311] = {
    'name': 'K(0)*(700)0', 'mass': 0.838, 'width': 0.463,
    'J': 0, 'I': 0.5, 'P': 1,
    'B': 0, 'Q': 0, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1, 'degeneracy_override': 0,
}
PARTICLES[9000321] = {
    'name': 'K(0)*(700)+', 'mass': 0.838, 'width': 0.463,
    'J': 0, 'I': 0.5, 'P': 1,
    'B': 0, 'Q': 1, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1, 'degeneracy_override': 0,
}

# K*(892)+: I(J^P) = 1/2(1-), mass = 891.67 MeV, width = 51.4 MeV
PARTICLES[323] = {
    'name': 'K*(892)+', 'mass': 0.89167, 'width': 0.0514,
    'J': 1, 'I': 0.5, 'P': -1,
    'B': 0, 'Q': 1, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# K*(892)0: I(J^P) = 1/2(1-), mass = 895.55 MeV, width = 47.3 MeV
PARTICLES[313] = {
    'name': 'K*(892)0', 'mass': 0.89555, 'width': 0.0473,
    'J': 1, 'I': 0.5, 'P': -1,
    'B': 0, 'Q': 0, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# K1(1270)0: I(J^P) = 1/2(1+), mass = 1253 MeV, width = 90 MeV
PARTICLES[10313] = {
    'name': 'K(1)(1270)0', 'mass': 1.253, 'width': 0.09,
    'J': 1, 'I': 0.5, 'P': 1,
    'B': 0, 'Q': 0, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[10323] = {
    'name': 'K(1)(1270)+', 'mass': 1.253, 'width': 0.09,
    'J': 1, 'I': 0.5, 'P': 1,
    'B': 0, 'Q': 1, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# K1(1400)0: I(J^P) = 1/2(1+), mass = 1403 MeV, width = 174 MeV
PARTICLES[20313] = {
    'name': 'K(1)(1400)0', 'mass': 1.403, 'width': 0.174,
    'J': 1, 'I': 0.5, 'P': 1,
    'B': 0, 'Q': 0, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[20323] = {
    'name': 'K(1)(1400)+', 'mass': 1.403, 'width': 0.174,
    'J': 1, 'I': 0.5, 'P': 1,
    'B': 0, 'Q': 1, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# K*(1410)0: I(J^P) = 1/2(1-), mass = 1414 MeV, width = 232 MeV
PARTICLES[100313] = {
    'name': 'K*(1410)0', 'mass': 1.414, 'width': 0.232,
    'J': 1, 'I': 0.5, 'P': -1,
    'B': 0, 'Q': 0, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[100323] = {
    'name': 'K*(1410)+', 'mass': 1.414, 'width': 0.232,
    'J': 1, 'I': 0.5, 'P': -1,
    'B': 0, 'Q': 1, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# K2*(1430)+: I(J^P) = 1/2(2+), mass = 1427.3 MeV, width = 100 MeV
PARTICLES[325] = {
    'name': 'K(2)*(1430)+', 'mass': 1.4273, 'width': 0.1,
    'J': 2, 'I': 0.5, 'P': 1,
    'B': 0, 'Q': 1, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# K0*(1430)0: I(J^P) = 1/2(0+), mass = 1430 MeV, width = 270 MeV
PARTICLES[10311] = {
    'name': 'K(0)*(1430)0', 'mass': 1.43, 'width': 0.27,
    'J': 0, 'I': 0.5, 'P': 1,
    'B': 0, 'Q': 0, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[10321] = {
    'name': 'K(0)*(1430)+', 'mass': 1.43, 'width': 0.27,
    'J': 0, 'I': 0.5, 'P': 1,
    'B': 0, 'Q': 1, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# K2*(1430)0: I(J^P) = 1/2(2+), mass = 1432.4 MeV, width = 109 MeV
PARTICLES[315] = {
    'name': 'K(2)*(1430)0', 'mass': 1.4324, 'width': 0.109,
    'J': 2, 'I': 0.5, 'P': 1,
    'B': 0, 'Q': 0, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# K1(1650)0: I(J^P) = 1/2(1+), mass = 1650 MeV, width = 150 MeV [NEW in PDG2025]
PARTICLES[9000313] = {
    'name': 'K(1)(1650)0', 'mass': 1.65, 'width': 0.15,
    'J': 1, 'I': 0.5, 'P': 1,
    'B': 0, 'Q': 0, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[9000323] = {
    'name': 'K(1)(1650)+', 'mass': 1.65, 'width': 0.15,
    'J': 1, 'I': 0.5, 'P': 1,
    'B': 0, 'Q': 1, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# K*(1680)0: I(J^P) = 1/2(1-), mass = 1718 MeV, width = 320 MeV
PARTICLES[30313] = {
    'name': 'K*(1680)0', 'mass': 1.718, 'width': 0.32,
    'J': 1, 'I': 0.5, 'P': -1,
    'B': 0, 'Q': 0, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[30323] = {
    'name': 'K*(1680)+', 'mass': 1.718, 'width': 0.32,
    'J': 1, 'I': 0.5, 'P': -1,
    'B': 0, 'Q': 1, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# K2(1770)0: I(J^P) = 1/2(2-), mass = 1773 MeV, width = 186 MeV
PARTICLES[10315] = {
    'name': 'K(2)(1770)0', 'mass': 1.773, 'width': 0.186,
    'J': 2, 'I': 0.5, 'P': -1,
    'B': 0, 'Q': 0, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[10325] = {
    'name': 'K(2)(1770)+', 'mass': 1.773, 'width': 0.186,
    'J': 2, 'I': 0.5, 'P': -1,
    'B': 0, 'Q': 1, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# K3*(1780)0: I(J^P) = 1/2(3-), mass = 1779 MeV, width = 161 MeV
PARTICLES[317] = {
    'name': 'K(3)*(1780)0', 'mass': 1.779, 'width': 0.161,
    'J': 3, 'I': 0.5, 'P': -1,
    'B': 0, 'Q': 0, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[327] = {
    'name': 'K(3)*(1780)+', 'mass': 1.779, 'width': 0.161,
    'J': 3, 'I': 0.5, 'P': -1,
    'B': 0, 'Q': 1, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# K2(1820)0: I(J^P) = 1/2(2-), mass = 1819 MeV, width = 264 MeV
PARTICLES[20315] = {
    'name': 'K(2)(1820)0', 'mass': 1.819, 'width': 0.264,
    'J': 2, 'I': 0.5, 'P': -1,
    'B': 0, 'Q': 0, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[20325] = {
    'name': 'K(2)(1820)+', 'mass': 1.819, 'width': 0.264,
    'J': 2, 'I': 0.5, 'P': -1,
    'B': 0, 'Q': 1, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# K0*(1950)0: I(J^P) = 1/2(0+), mass = 1957 MeV, width = 170 MeV [NEW in PDG2025]
PARTICLES[9020311] = {
    'name': 'K(0)*(1950)0', 'mass': 1.957, 'width': 0.17,
    'J': 0, 'I': 0.5, 'P': 1,
    'B': 0, 'Q': 0, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[9020321] = {
    'name': 'K(0)*(1950)+', 'mass': 1.957, 'width': 0.17,
    'J': 0, 'I': 0.5, 'P': 1,
    'B': 0, 'Q': 1, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# K2*(1980)0: I(J^P) = 1/2(2+), mass = 1990 MeV, width = 348 MeV [NEW in PDG2025]
PARTICLES[9010315] = {
    'name': 'K(2)*(1980)0', 'mass': 1.99, 'width': 0.348,
    'J': 2, 'I': 0.5, 'P': 1,
    'B': 0, 'Q': 0, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[9010325] = {
    'name': 'K(2)*(1980)+', 'mass': 1.99, 'width': 0.348,
    'J': 2, 'I': 0.5, 'P': 1,
    'B': 0, 'Q': 1, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}

# K4*(2045)0: I(J^P) = 1/2(4+), mass = 2048 MeV, width = 199 MeV
PARTICLES[319] = {
    'name': 'K(4)*(2045)0', 'mass': 2.048, 'width': 0.199,
    'J': 4, 'I': 0.5, 'P': 1,
    'B': 0, 'Q': 0, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}
PARTICLES[329] = {
    'name': 'K(4)*(2045)+', 'mass': 2.048, 'width': 0.199,
    'J': 4, 'I': 0.5, 'P': 1,
    'B': 0, 'Q': 1, 'S': 1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': -1,
}


# =====================================================================
# N BARYONS (from rpp2025-tab-baryons-N.pdf)
# =====================================================================

# p: I(J^P) = 1/2(1/2+), mass = 938.272 MeV
PARTICLES[2212] = {
    'name': 'p', 'mass': 0.938272, 'width': 0,
    'J': 0.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 1, 'statistics': 1,
}

# n: I(J^P) = 1/2(1/2+), mass = 939.565 MeV, width = 7.493e-28 GeV (tau = 878.4 s)
PARTICLES[2112] = {
    'name': 'n', 'mass': 0.939565, 'width': 7.493e-28,
    'J': 0.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 1, 'statistics': 1,
}

# N(1440)0 (Roper): I(J^P) = 1/2(1/2+), mass = 1440 MeV, width = 350 MeV
PARTICLES[12112] = {
    'name': 'N(1440)0', 'mass': 1.44, 'width': 0.35,
    'J': 0.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[12212] = {
    'name': 'N(1440)+', 'mass': 1.44, 'width': 0.35,
    'J': 0.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# N(1520)0: I(J^P) = 1/2(3/2-), mass = 1515 MeV, width = 110 MeV
PARTICLES[1214] = {
    'name': 'N(1520)0', 'mass': 1.515, 'width': 0.11,
    'J': 1.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[2124] = {
    'name': 'N(1520)+', 'mass': 1.515, 'width': 0.11,
    'J': 1.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# N(1535)0: I(J^P) = 1/2(1/2-), mass = 1530 MeV, width = 150 MeV
PARTICLES[22112] = {
    'name': 'N(1535)0', 'mass': 1.53, 'width': 0.15,
    'J': 0.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[22212] = {
    'name': 'N(1535)+', 'mass': 1.53, 'width': 0.15,
    'J': 0.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# N(1650)0: I(J^P) = 1/2(1/2-), mass = 1650 MeV, width = 125 MeV
PARTICLES[32112] = {
    'name': 'N(1650)0', 'mass': 1.65, 'width': 0.125,
    'J': 0.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[32212] = {
    'name': 'N(1650)+', 'mass': 1.65, 'width': 0.125,
    'J': 0.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# N(1675)0: I(J^P) = 1/2(5/2-), mass = 1675 MeV, width = 145 MeV
PARTICLES[2116] = {
    'name': 'N(1675)0', 'mass': 1.675, 'width': 0.145,
    'J': 2.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[2216] = {
    'name': 'N(1675)+', 'mass': 1.675, 'width': 0.145,
    'J': 2.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# N(1680)0: I(J^P) = 1/2(5/2+), mass = 1685 MeV, width = 120 MeV
PARTICLES[12116] = {
    'name': 'N(1680)0', 'mass': 1.685, 'width': 0.12,
    'J': 2.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[12216] = {
    'name': 'N(1680)+', 'mass': 1.685, 'width': 0.12,
    'J': 2.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# N(1700)0: I(J^P) = 1/2(3/2-), mass = 1720 MeV, width = 200 MeV
PARTICLES[21214] = {
    'name': 'N(1700)0', 'mass': 1.72, 'width': 0.2,
    'J': 1.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[22124] = {
    'name': 'N(1700)+', 'mass': 1.72, 'width': 0.2,
    'J': 1.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# N(1710)0: I(J^P) = 1/2(1/2+), mass = 1710 MeV, width = 140 MeV
PARTICLES[42112] = {
    'name': 'N(1710)0', 'mass': 1.71, 'width': 0.14,
    'J': 0.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[42212] = {
    'name': 'N(1710)+', 'mass': 1.71, 'width': 0.14,
    'J': 0.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# N(1720)0: I(J^P) = 1/2(3/2+), mass = 1720 MeV, width = 250 MeV
PARTICLES[31214] = {
    'name': 'N(1720)0', 'mass': 1.72, 'width': 0.25,
    'J': 1.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[32124] = {
    'name': 'N(1720)+', 'mass': 1.72, 'width': 0.25,
    'J': 1.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# N(1875)0: I(J^P) = 1/2(3/2-), mass = 1875 MeV, width = 200 MeV
PARTICLES[9902114] = {
    'name': 'N(1875)0', 'mass': 1.875, 'width': 0.2,
    'J': 1.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[9902214] = {
    'name': 'N(1875)+', 'mass': 1.875, 'width': 0.2,
    'J': 1.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# N(1880)0: I(J^P) = 1/2(1/2+), mass = 1880 MeV, width = 300 MeV
PARTICLES[9902112] = {
    'name': 'N(1880)0', 'mass': 1.88, 'width': 0.3,
    'J': 0.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[9902212] = {
    'name': 'N(1880)+', 'mass': 1.88, 'width': 0.3,
    'J': 0.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# N(1895)0: I(J^P) = 1/2(1/2-), mass = 1895 MeV, width = 120 MeV
PARTICLES[9912112] = {
    'name': 'N(1895)0', 'mass': 1.895, 'width': 0.12,
    'J': 0.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[9912212] = {
    'name': 'N(1895)+', 'mass': 1.895, 'width': 0.12,
    'J': 0.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# N(1900)0: I(J^P) = 1/2(3/2+), mass = 1920 MeV, width = 200 MeV
PARTICLES[9912114] = {
    'name': 'N(1900)0', 'mass': 1.92, 'width': 0.2,
    'J': 1.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[9912214] = {
    'name': 'N(1900)+', 'mass': 1.92, 'width': 0.2,
    'J': 1.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# N(2060)0: I(J^P) = 1/2(5/2-), mass = 2100 MeV, width = 400 MeV
PARTICLES[9902116] = {
    'name': 'N(2060)0', 'mass': 2.1, 'width': 0.4,
    'J': 2.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[9902216] = {
    'name': 'N(2060)+', 'mass': 2.1, 'width': 0.4,
    'J': 2.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# N(2100)0: I(J^P) = 1/2(1/2+), mass = 2100 MeV, width = 260 MeV
PARTICLES[9922112] = {
    'name': 'N(2100)0', 'mass': 2.1, 'width': 0.26,
    'J': 0.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[9922212] = {
    'name': 'N(2100)+', 'mass': 2.1, 'width': 0.26,
    'J': 0.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# N(2120)0: I(J^P) = 1/2(3/2-), mass = 2120 MeV, width = 300 MeV
PARTICLES[9922114] = {
    'name': 'N(2120)0', 'mass': 2.12, 'width': 0.3,
    'J': 1.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[9922214] = {
    'name': 'N(2120)+', 'mass': 2.12, 'width': 0.3,
    'J': 1.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# N(2190)0: I(J^P) = 1/2(7/2-), mass = 2180 MeV, width = 400 MeV
PARTICLES[1218] = {
    'name': 'N(2190)0', 'mass': 2.18, 'width': 0.4,
    'J': 3.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[2128] = {
    'name': 'N(2190)+', 'mass': 2.18, 'width': 0.4,
    'J': 3.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# N(2220)0: I(J^P) = 1/2(9/2+), mass = 2225 MeV, width = 400 MeV
PARTICLES[99021110] = {
    'name': 'N(2220)0', 'mass': 2.225, 'width': 0.4,
    'J': 4.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[99022110] = {
    'name': 'N(2220)+', 'mass': 2.225, 'width': 0.4,
    'J': 4.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# N(2250)0: I(J^P) = 1/2(9/2-), mass = 2280 MeV, width = 500 MeV
PARTICLES[99121110] = {
    'name': 'N(2250)0', 'mass': 2.28, 'width': 0.5,
    'J': 4.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[99122110] = {
    'name': 'N(2250)+', 'mass': 2.28, 'width': 0.5,
    'J': 4.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# N(2600)0: I(J^P) = 1/2(11/2-), mass = 2600 MeV, width = 650 MeV
PARTICLES[99021112] = {
    'name': 'N(2600)0', 'mass': 2.6, 'width': 0.65,
    'J': 5.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[99022112] = {
    'name': 'N(2600)+', 'mass': 2.6, 'width': 0.65,
    'J': 5.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# =====================================================================
# DELTA BARYONS (from rpp2025-tab-baryons-Delta.pdf)
# =====================================================================

# Delta(1232): I(J^P) = 3/2(3/2+), mass = 1232 MeV, width = 117 MeV
PARTICLES[1114] = {
    'name': 'Delta(1232)-', 'mass': 1.232, 'width': 0.117,
    'J': 1.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': -1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[2114] = {
    'name': 'Delta(1232)0', 'mass': 1.232, 'width': 0.117,
    'J': 1.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[2214] = {
    'name': 'Delta(1232)+', 'mass': 1.232, 'width': 0.117,
    'J': 1.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[2224] = {
    'name': 'Delta(1232)++', 'mass': 1.232, 'width': 0.117,
    'J': 1.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': 2, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Delta(1600): I(J^P) = 3/2(3/2+), mass = 1570 MeV, width = 250 MeV
PARTICLES[31114] = {
    'name': 'Delta(1600)-', 'mass': 1.57, 'width': 0.25,
    'J': 1.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': -1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[32114] = {
    'name': 'Delta(1600)0', 'mass': 1.57, 'width': 0.25,
    'J': 1.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[32214] = {
    'name': 'Delta(1600)+', 'mass': 1.57, 'width': 0.25,
    'J': 1.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[32224] = {
    'name': 'Delta(1600)++', 'mass': 1.57, 'width': 0.25,
    'J': 1.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': 2, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Delta(1620): I(J^P) = 3/2(1/2-), mass = 1610 MeV, width = 130 MeV
PARTICLES[1112] = {
    'name': 'Delta(1620)-', 'mass': 1.61, 'width': 0.13,
    'J': 0.5, 'I': 1.5, 'P': -1,
    'B': 1, 'Q': -1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[1212] = {
    'name': 'Delta(1620)0', 'mass': 1.61, 'width': 0.13,
    'J': 0.5, 'I': 1.5, 'P': -1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[2122] = {
    'name': 'Delta(1620)+', 'mass': 1.61, 'width': 0.13,
    'J': 0.5, 'I': 1.5, 'P': -1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[2222] = {
    'name': 'Delta(1620)++', 'mass': 1.61, 'width': 0.13,
    'J': 0.5, 'I': 1.5, 'P': -1,
    'B': 1, 'Q': 2, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Delta(1700): I(J^P) = 3/2(3/2-), mass = 1710 MeV, width = 300 MeV
PARTICLES[11114] = {
    'name': 'Delta(1700)-', 'mass': 1.71, 'width': 0.3,
    'J': 1.5, 'I': 1.5, 'P': -1,
    'B': 1, 'Q': -1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[12114] = {
    'name': 'Delta(1700)0', 'mass': 1.71, 'width': 0.3,
    'J': 1.5, 'I': 1.5, 'P': -1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[12214] = {
    'name': 'Delta(1700)+', 'mass': 1.71, 'width': 0.3,
    'J': 1.5, 'I': 1.5, 'P': -1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[12224] = {
    'name': 'Delta(1700)++', 'mass': 1.71, 'width': 0.3,
    'J': 1.5, 'I': 1.5, 'P': -1,
    'B': 1, 'Q': 2, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Delta(1900): I(J^P) = 3/2(1/2-), mass = 1860 MeV, width = 250 MeV
PARTICLES[11112] = {
    'name': 'Delta(1900)-', 'mass': 1.86, 'width': 0.25,
    'J': 0.5, 'I': 1.5, 'P': -1,
    'B': 1, 'Q': -1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[11212] = {
    'name': 'Delta(1900)0', 'mass': 1.86, 'width': 0.25,
    'J': 0.5, 'I': 1.5, 'P': -1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[12122] = {
    'name': 'Delta(1900)+', 'mass': 1.86, 'width': 0.25,
    'J': 0.5, 'I': 1.5, 'P': -1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[12222] = {
    'name': 'Delta(1900)++', 'mass': 1.86, 'width': 0.25,
    'J': 0.5, 'I': 1.5, 'P': -1,
    'B': 1, 'Q': 2, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Delta(1905): I(J^P) = 3/2(5/2+), mass = 1880 MeV, width = 330 MeV
PARTICLES[1116] = {
    'name': 'Delta(1905)-', 'mass': 1.88, 'width': 0.33,
    'J': 2.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': -1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[1216] = {
    'name': 'Delta(1905)0', 'mass': 1.88, 'width': 0.33,
    'J': 2.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[2126] = {
    'name': 'Delta(1905)+', 'mass': 1.88, 'width': 0.33,
    'J': 2.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[2226] = {
    'name': 'Delta(1905)++', 'mass': 1.88, 'width': 0.33,
    'J': 2.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': 2, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Delta(1910): I(J^P) = 3/2(1/2+), mass = 1900 MeV, width = 300 MeV
PARTICLES[21112] = {
    'name': 'Delta(1910)-', 'mass': 1.9, 'width': 0.3,
    'J': 0.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': -1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[21212] = {
    'name': 'Delta(1910)0', 'mass': 1.9, 'width': 0.3,
    'J': 0.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[22122] = {
    'name': 'Delta(1910)+', 'mass': 1.9, 'width': 0.3,
    'J': 0.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[22222] = {
    'name': 'Delta(1910)++', 'mass': 1.9, 'width': 0.3,
    'J': 0.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': 2, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Delta(1920): I(J^P) = 3/2(3/2+), mass = 1920 MeV, width = 300 MeV
PARTICLES[21114] = {
    'name': 'Delta(1920)-', 'mass': 1.92, 'width': 0.3,
    'J': 1.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': -1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[22114] = {
    'name': 'Delta(1920)0', 'mass': 1.92, 'width': 0.3,
    'J': 1.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[22214] = {
    'name': 'Delta(1920)+', 'mass': 1.92, 'width': 0.3,
    'J': 1.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[22224] = {
    'name': 'Delta(1920)++', 'mass': 1.92, 'width': 0.3,
    'J': 1.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': 2, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Delta(1930): I(J^P) = 3/2(5/2-), mass = 1950 MeV, width = 300 MeV
PARTICLES[11116] = {
    'name': 'Delta(1930)-', 'mass': 1.95, 'width': 0.3,
    'J': 2.5, 'I': 1.5, 'P': -1,
    'B': 1, 'Q': -1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[11216] = {
    'name': 'Delta(1930)0', 'mass': 1.95, 'width': 0.3,
    'J': 2.5, 'I': 1.5, 'P': -1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[12126] = {
    'name': 'Delta(1930)+', 'mass': 1.95, 'width': 0.3,
    'J': 2.5, 'I': 1.5, 'P': -1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[12226] = {
    'name': 'Delta(1930)++', 'mass': 1.95, 'width': 0.3,
    'J': 2.5, 'I': 1.5, 'P': -1,
    'B': 1, 'Q': 2, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Delta(1950): I(J^P) = 3/2(7/2+), mass = 1930 MeV, width = 280 MeV
PARTICLES[1118] = {
    'name': 'Delta(1950)-', 'mass': 1.93, 'width': 0.28,
    'J': 3.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': -1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[2118] = {
    'name': 'Delta(1950)0', 'mass': 1.93, 'width': 0.28,
    'J': 3.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[2218] = {
    'name': 'Delta(1950)+', 'mass': 1.93, 'width': 0.28,
    'J': 3.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[2228] = {
    'name': 'Delta(1950)++', 'mass': 1.93, 'width': 0.28,
    'J': 3.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': 2, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Delta(2200): I(J^P) = 3/2(7/2-), mass = 2200 MeV, width = 350 MeV
PARTICLES[9901118] = {
    'name': 'Delta(2200)-', 'mass': 2.2, 'width': 0.35,
    'J': 3.5, 'I': 1.5, 'P': -1,
    'B': 1, 'Q': -1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[9901218] = {
    'name': 'Delta(2200)0', 'mass': 2.2, 'width': 0.35,
    'J': 3.5, 'I': 1.5, 'P': -1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[9902128] = {
    'name': 'Delta(2200)+', 'mass': 2.2, 'width': 0.35,
    'J': 3.5, 'I': 1.5, 'P': -1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[9902228] = {
    'name': 'Delta(2200)++', 'mass': 2.2, 'width': 0.35,
    'J': 3.5, 'I': 1.5, 'P': -1,
    'B': 1, 'Q': 2, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Delta(2420): I(J^P) = 3/2(11/2+), mass = 2420 MeV, width = 400 MeV
PARTICLES[99011112] = {
    'name': 'Delta(2420)-', 'mass': 2.42, 'width': 0.4,
    'J': 5.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': -1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[99012112] = {
    'name': 'Delta(2420)0', 'mass': 2.42, 'width': 0.4,
    'J': 5.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': 0, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[99021212] = {
    'name': 'Delta(2420)+', 'mass': 2.42, 'width': 0.4,
    'J': 5.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': 1, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}
PARTICLES[99022212] = {
    'name': 'Delta(2420)++', 'mass': 2.42, 'width': 0.4,
    'J': 5.5, 'I': 1.5, 'P': 1,
    'B': 1, 'Q': 2, 'S': 0, 'charm': 0, 'absS': 0, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# ==========================================
# Lambda Baryons
# ==========================================

# Lambda: I(J^P) = 0(1/2+), mass = 1115.68 MeV, width = 2.515e-15 GeV
PARTICLES[3122] = {
    'name': 'Lambda', 'mass': 1.11568, 'width': 2.515e-15,
    'J': 0.5, 'I': 0, 'P': 1,
    'B': 1, 'Q': 0, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 1, 'statistics': 1,
}

# Lambda(1405): I(J^P) = 0(1/2-), mass = 1405.1 MeV, width = 50.5 MeV
PARTICLES[13122] = {
    'name': 'Lambda(1405)', 'mass': 1.4051, 'width': 0.0505,
    'J': 0.5, 'I': 0, 'P': -1,
    'B': 1, 'Q': 0, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Lambda(1520): I(J^P) = 0(3/2-), mass = 1519 MeV, width = 16 MeV
PARTICLES[3124] = {
    'name': 'Lambda(1520)', 'mass': 1.519, 'width': 0.016,
    'J': 1.5, 'I': 0, 'P': -1,
    'B': 1, 'Q': 0, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Lambda(1600): I(J^P) = 0(1/2+), mass = 1600 MeV, width = 200 MeV
PARTICLES[23122] = {
    'name': 'Lambda(1600)', 'mass': 1.6, 'width': 0.2,
    'J': 0.5, 'I': 0, 'P': 1,
    'B': 1, 'Q': 0, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Lambda(1670): I(J^P) = 0(1/2-), mass = 1674 MeV, width = 30 MeV
PARTICLES[33122] = {
    'name': 'Lambda(1670)', 'mass': 1.674, 'width': 0.03,
    'J': 0.5, 'I': 0, 'P': -1,
    'B': 1, 'Q': 0, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Lambda(1690): I(J^P) = 0(3/2-), mass = 1690 MeV, width = 70 MeV
PARTICLES[13124] = {
    'name': 'Lambda(1690)', 'mass': 1.69, 'width': 0.07,
    'J': 1.5, 'I': 0, 'P': -1,
    'B': 1, 'Q': 0, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Lambda(1800): I(J^P) = 0(1/2-), mass = 1800 MeV, width = 200 MeV
PARTICLES[43122] = {
    'name': 'Lambda(1800)', 'mass': 1.8, 'width': 0.2,
    'J': 0.5, 'I': 0, 'P': -1,
    'B': 1, 'Q': 0, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Lambda(1810): I(J^P) = 0(1/2+), mass = 1790 MeV, width = 110 MeV
PARTICLES[53122] = {
    'name': 'Lambda(1810)', 'mass': 1.79, 'width': 0.11,
    'J': 0.5, 'I': 0, 'P': 1,
    'B': 1, 'Q': 0, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Lambda(1820): I(J^P) = 0(5/2+), mass = 1820 MeV, width = 80 MeV
PARTICLES[3126] = {
    'name': 'Lambda(1820)', 'mass': 1.82, 'width': 0.08,
    'J': 2.5, 'I': 0, 'P': 1,
    'B': 1, 'Q': 0, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Lambda(1830): I(J^P) = 0(5/2-), mass = 1825 MeV, width = 90 MeV
PARTICLES[13126] = {
    'name': 'Lambda(1830)', 'mass': 1.825, 'width': 0.09,
    'J': 2.5, 'I': 0, 'P': -1,
    'B': 1, 'Q': 0, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Lambda(1890): I(J^P) = 0(3/2+), mass = 1890 MeV, width = 120 MeV
PARTICLES[23124] = {
    'name': 'Lambda(1890)', 'mass': 1.89, 'width': 0.12,
    'J': 1.5, 'I': 0, 'P': 1,
    'B': 1, 'Q': 0, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Lambda(2110): I(J^P) = 0(5/2+), mass = 2090 MeV, width = 250 MeV
PARTICLES[23126] = {
    'name': 'Lambda(2110)', 'mass': 2.09, 'width': 0.25,
    'J': 2.5, 'I': 0, 'P': 1,
    'B': 1, 'Q': 0, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Lambda(2100): I(J^P) = 0(7/2-), mass = 2100 MeV, width = 200 MeV
PARTICLES[3128] = {
    'name': 'Lambda(2100)', 'mass': 2.1, 'width': 0.2,
    'J': 3.5, 'I': 0, 'P': -1,
    'B': 1, 'Q': 0, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Lambda(2350): I(J^P) = 0(9/2+), mass = 2350 MeV, width = 150 MeV
PARTICLES[99031210] = {
    'name': 'Lambda(2350)', 'mass': 2.35, 'width': 0.15,
    'J': 4.5, 'I': 0, 'P': 1,
    'B': 1, 'Q': 0, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# ==========================================
# Sigma Baryons
# ==========================================

# Sigma+: I(J^P) = 1(1/2+), mass = 1189.37 MeV, width = 8.209e-15 GeV
PARTICLES[3222] = {
    'name': 'Sigma+', 'mass': 1.18937, 'width': 8.209e-15,
    'J': 0.5, 'I': 1, 'P': 1,
    'B': 1, 'Q': 1, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 1, 'statistics': 1,
}

# Sigma0: I(J^P) = 1(1/2+), mass = 1192.64 MeV, width = 8.9e-06 GeV
PARTICLES[3212] = {
    'name': 'Sigma0', 'mass': 1.19264, 'width': 8.9e-06,
    'J': 0.5, 'I': 1, 'P': 1,
    'B': 1, 'Q': 0, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma-: I(J^P) = 1(1/2+), mass = 1197.45 MeV, width = 4.45e-15 GeV
PARTICLES[3112] = {
    'name': 'Sigma-', 'mass': 1.19745, 'width': 4.45e-15,
    'J': 0.5, 'I': 1, 'P': 1,
    'B': 1, 'Q': -1, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 1, 'statistics': 1,
}

# Sigma(1385)+: I(J^P) = 1(3/2+), mass = 1382.83 MeV, width = 36.2 MeV
PARTICLES[3224] = {
    'name': 'Sigma(1385)+', 'mass': 1.38283, 'width': 0.0362,
    'J': 1.5, 'I': 1, 'P': 1,
    'B': 1, 'Q': 1, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma(1385)0: I(J^P) = 1(3/2+), mass = 1383.7 MeV, width = 36 MeV
PARTICLES[3214] = {
    'name': 'Sigma(1385)0', 'mass': 1.3837, 'width': 0.036,
    'J': 1.5, 'I': 1, 'P': 1,
    'B': 1, 'Q': 0, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma(1385)-: I(J^P) = 1(3/2+), mass = 1387.2 MeV, width = 39.4 MeV
PARTICLES[3114] = {
    'name': 'Sigma(1385)-', 'mass': 1.3872, 'width': 0.0394,
    'J': 1.5, 'I': 1, 'P': 1,
    'B': 1, 'Q': -1, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma(1660)+: I(J^P) = 1(1/2+), mass = 1660 MeV, width = 200 MeV
PARTICLES[13222] = {
    'name': 'Sigma(1660)+', 'mass': 1.66, 'width': 0.2,
    'J': 0.5, 'I': 1, 'P': 1,
    'B': 1, 'Q': 1, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma(1660)0: I(J^P) = 1(1/2+), mass = 1660 MeV, width = 200 MeV
PARTICLES[13212] = {
    'name': 'Sigma(1660)0', 'mass': 1.66, 'width': 0.2,
    'J': 0.5, 'I': 1, 'P': 1,
    'B': 1, 'Q': 0, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma(1660)-: I(J^P) = 1(1/2+), mass = 1660 MeV, width = 200 MeV
PARTICLES[13112] = {
    'name': 'Sigma(1660)-', 'mass': 1.66, 'width': 0.2,
    'J': 0.5, 'I': 1, 'P': 1,
    'B': 1, 'Q': -1, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma(1670)+: I(J^P) = 1(3/2-), mass = 1675 MeV, width = 70 MeV
PARTICLES[13224] = {
    'name': 'Sigma(1670)+', 'mass': 1.675, 'width': 0.07,
    'J': 1.5, 'I': 1, 'P': -1,
    'B': 1, 'Q': 1, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma(1670)0: I(J^P) = 1(3/2-), mass = 1675 MeV, width = 70 MeV
PARTICLES[13214] = {
    'name': 'Sigma(1670)0', 'mass': 1.675, 'width': 0.07,
    'J': 1.5, 'I': 1, 'P': -1,
    'B': 1, 'Q': 0, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma(1670)-: I(J^P) = 1(3/2-), mass = 1675 MeV, width = 70 MeV
PARTICLES[13114] = {
    'name': 'Sigma(1670)-', 'mass': 1.675, 'width': 0.07,
    'J': 1.5, 'I': 1, 'P': -1,
    'B': 1, 'Q': -1, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma(1750)+: I(J^P) = 1(1/2-), mass = 1750 MeV, width = 150 MeV
PARTICLES[23222] = {
    'name': 'Sigma(1750)+', 'mass': 1.75, 'width': 0.15,
    'J': 0.5, 'I': 1, 'P': -1,
    'B': 1, 'Q': 1, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma(1750)0: I(J^P) = 1(1/2-), mass = 1750 MeV, width = 150 MeV
PARTICLES[23212] = {
    'name': 'Sigma(1750)0', 'mass': 1.75, 'width': 0.15,
    'J': 0.5, 'I': 1, 'P': -1,
    'B': 1, 'Q': 0, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma(1750)-: I(J^P) = 1(1/2-), mass = 1750 MeV, width = 150 MeV
PARTICLES[23112] = {
    'name': 'Sigma(1750)-', 'mass': 1.75, 'width': 0.15,
    'J': 0.5, 'I': 1, 'P': -1,
    'B': 1, 'Q': -1, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma(1775)+: I(J^P) = 1(5/2-), mass = 1775 MeV, width = 120 MeV
PARTICLES[3226] = {
    'name': 'Sigma(1775)+', 'mass': 1.775, 'width': 0.12,
    'J': 2.5, 'I': 1, 'P': -1,
    'B': 1, 'Q': 1, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma(1775)0: I(J^P) = 1(5/2-), mass = 1775 MeV, width = 120 MeV
PARTICLES[3216] = {
    'name': 'Sigma(1775)0', 'mass': 1.775, 'width': 0.12,
    'J': 2.5, 'I': 1, 'P': -1,
    'B': 1, 'Q': 0, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma(1775)-: I(J^P) = 1(5/2-), mass = 1775 MeV, width = 120 MeV
PARTICLES[3116] = {
    'name': 'Sigma(1775)-', 'mass': 1.775, 'width': 0.12,
    'J': 2.5, 'I': 1, 'P': -1,
    'B': 1, 'Q': -1, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma(1910)+: I(J^P) = 1(3/2+), mass = 1910 MeV, width = 220 MeV
PARTICLES[23224] = {
    'name': 'Sigma(1910)+', 'mass': 1.91, 'width': 0.22,
    'J': 1.5, 'I': 1, 'P': 1,
    'B': 1, 'Q': 1, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma(1910)0: I(J^P) = 1(3/2+), mass = 1910 MeV, width = 220 MeV
PARTICLES[23214] = {
    'name': 'Sigma(1910)0', 'mass': 1.91, 'width': 0.22,
    'J': 1.5, 'I': 1, 'P': 1,
    'B': 1, 'Q': 0, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma(1910)-: I(J^P) = 1(3/2+), mass = 1910 MeV, width = 220 MeV
PARTICLES[23114] = {
    'name': 'Sigma(1910)-', 'mass': 1.91, 'width': 0.22,
    'J': 1.5, 'I': 1, 'P': 1,
    'B': 1, 'Q': -1, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma(1915)+: I(J^P) = 1(5/2+), mass = 1915 MeV, width = 120 MeV
PARTICLES[13226] = {
    'name': 'Sigma(1915)+', 'mass': 1.915, 'width': 0.12,
    'J': 2.5, 'I': 1, 'P': 1,
    'B': 1, 'Q': 1, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma(1915)0: I(J^P) = 1(5/2+), mass = 1915 MeV, width = 120 MeV
PARTICLES[13216] = {
    'name': 'Sigma(1915)0', 'mass': 1.915, 'width': 0.12,
    'J': 2.5, 'I': 1, 'P': 1,
    'B': 1, 'Q': 0, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma(1915)-: I(J^P) = 1(5/2+), mass = 1915 MeV, width = 120 MeV
PARTICLES[13116] = {
    'name': 'Sigma(1915)-', 'mass': 1.915, 'width': 0.12,
    'J': 2.5, 'I': 1, 'P': 1,
    'B': 1, 'Q': -1, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma(2030)+: I(J^P) = 1(7/2-), mass = 2030 MeV, width = 180 MeV
PARTICLES[3228] = {
    'name': 'Sigma(2030)+', 'mass': 2.03, 'width': 0.18,
    'J': 3.5, 'I': 1, 'P': -1,
    'B': 1, 'Q': 1, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma(2030)0: I(J^P) = 1(7/2-), mass = 2030 MeV, width = 180 MeV
PARTICLES[3218] = {
    'name': 'Sigma(2030)0', 'mass': 2.03, 'width': 0.18,
    'J': 3.5, 'I': 1, 'P': -1,
    'B': 1, 'Q': 0, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Sigma(2030)-: I(J^P) = 1(7/2-), mass = 2030 MeV, width = 180 MeV
PARTICLES[3118] = {
    'name': 'Sigma(2030)-', 'mass': 2.03, 'width': 0.18,
    'J': 3.5, 'I': 1, 'P': -1,
    'B': 1, 'Q': -1, 'S': -1, 'charm': 0, 'absS': 1, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# ==========================================
# Xi Baryons
# ==========================================

# Xi0: I(J^P) = 1/2(1/2+), mass = 1314.86 MeV, width = 2.27e-15 GeV
PARTICLES[3322] = {
    'name': 'Xi0', 'mass': 1.31486, 'width': 2.27e-15,
    'J': 0.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': 0, 'S': -2, 'charm': 0, 'absS': 2, 'absC': 0,
    'stable': 1, 'statistics': 1,
}

# Xi-: I(J^P) = 1/2(1/2+), mass = 1321.71 MeV, width = 4.02e-15 GeV
PARTICLES[3312] = {
    'name': 'Xi-', 'mass': 1.32171, 'width': 4.02e-15,
    'J': 0.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': -1, 'S': -2, 'charm': 0, 'absS': 2, 'absC': 0,
    'stable': 1, 'statistics': 1,
}

# Xi(1530)0: I(J^P) = 1/2(3/2+), mass = 1531.8 MeV, width = 9.1 MeV
PARTICLES[3324] = {
    'name': 'Xi(1530)0', 'mass': 1.5318, 'width': 0.0091,
    'J': 1.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': 0, 'S': -2, 'charm': 0, 'absS': 2, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Xi(1530)-: I(J^P) = 1/2(3/2+), mass = 1535 MeV, width = 9.9 MeV
PARTICLES[3314] = {
    'name': 'Xi(1530)-', 'mass': 1.535, 'width': 0.0099,
    'J': 1.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': -1, 'S': -2, 'charm': 0, 'absS': 2, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Xi(1690)0: I(J^P) = 1/2(1/2-), mass = 1690 MeV, width = 20 MeV
PARTICLES[203322] = {
    'name': 'Xi(1690)0', 'mass': 1.69, 'width': 0.02,
    'J': 0.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 0, 'S': -2, 'charm': 0, 'absS': 2, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Xi(1690)-: I(J^P) = 1/2(1/2-), mass = 1690 MeV, width = 20 MeV
PARTICLES[203312] = {
    'name': 'Xi(1690)-', 'mass': 1.69, 'width': 0.02,
    'J': 0.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': -1, 'S': -2, 'charm': 0, 'absS': 2, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Xi(1820)0: I(J^P) = 1/2(3/2-), mass = 1823 MeV, width = 24 MeV
PARTICLES[13324] = {
    'name': 'Xi(1820)0', 'mass': 1.823, 'width': 0.024,
    'J': 1.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 0, 'S': -2, 'charm': 0, 'absS': 2, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Xi(1820)-: I(J^P) = 1/2(3/2-), mass = 1823 MeV, width = 24 MeV
PARTICLES[13314] = {
    'name': 'Xi(1820)-', 'mass': 1.823, 'width': 0.024,
    'J': 1.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': -1, 'S': -2, 'charm': 0, 'absS': 2, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Xi(1950)0: I(J^P) = 1/2(5/2-), mass = 1950 MeV, width = 60 MeV
PARTICLES[103326] = {
    'name': 'Xi(1950)0', 'mass': 1.95, 'width': 0.06,
    'J': 2.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': 0, 'S': -2, 'charm': 0, 'absS': 2, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Xi(1950)-: I(J^P) = 1/2(5/2-), mass = 1950 MeV, width = 60 MeV
PARTICLES[103316] = {
    'name': 'Xi(1950)-', 'mass': 1.95, 'width': 0.06,
    'J': 2.5, 'I': 0.5, 'P': -1,
    'B': 1, 'Q': -1, 'S': -2, 'charm': 0, 'absS': 2, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Xi(2030)0: I(J^P) = 1/2(5/2+), mass = 2025 MeV, width = 20 MeV
PARTICLES[203326] = {
    'name': 'Xi(2030)0', 'mass': 2.025, 'width': 0.02,
    'J': 2.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': 0, 'S': -2, 'charm': 0, 'absS': 2, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Xi(2030)-: I(J^P) = 1/2(5/2+), mass = 2025 MeV, width = 20 MeV
PARTICLES[203316] = {
    'name': 'Xi(2030)-', 'mass': 2.025, 'width': 0.02,
    'J': 2.5, 'I': 0.5, 'P': 1,
    'B': 1, 'Q': -1, 'S': -2, 'charm': 0, 'absS': 2, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# ==========================================
# Omega Baryons
# ==========================================

# Omega-: I(J^P) = 0(3/2+), mass = 1672.45 MeV, width = 8.02e-15 GeV
PARTICLES[3334] = {
    'name': 'Omega', 'mass': 1.67245, 'width': 8.02e-15,
    'J': 1.5, 'I': 0, 'P': 1,
    'B': 1, 'Q': -1, 'S': -3, 'charm': 0, 'absS': 3, 'absC': 0,
    'stable': 1, 'statistics': 1,
}

# Omega(2012)-: I(J^P) = 0(3/2-), mass = 2012.4 MeV, width = 6.4 MeV
PARTICLES[9903334] = {
    'name': 'Omega(2012)', 'mass': 2.0124, 'width': 0.0064,
    'J': 1.5, 'I': 0, 'P': -1,
    'B': 1, 'Q': -1, 'S': -3, 'charm': 0, 'absS': 3, 'absC': 0,
    'stable': 0, 'statistics': 1,
}

# Omega(2250)-: I(J^P) = 0(7/2-), mass = 2252 MeV, width = 55 MeV
PARTICLES[203338] = {
    'name': 'Omega(2250)', 'mass': 2.252, 'width': 0.055,
    'J': 3.5, 'I': 0, 'P': -1,
    'B': 1, 'Q': -1, 'S': -3, 'charm': 0, 'absS': 3, 'absC': 0,
    'stable': 0, 'statistics': 1,
}


# =====================================================================
# Elementary particle masses (for threshold computation)
# =====================================================================
ELEMENTARY_MASSES = {
    22: 0, 11: 0.000511, -11: 0.000511, 13: 0.10566, -13: 0.10566,
    14: 0, -14: 0, 12: 0, -12: 0, 15: 1.777, -15: 1.777, 16: 0, -16: 0
}

# =====================================================================
# Parse existing decays.dat
# =====================================================================
def parse_decays_file(filepath):
    """Parse decays.dat and return dict of pdgid -> list of {br, daughters}"""
    decays = OrderedDict()
    with open(filepath, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        # Skip empty lines and comments
        if not line or line.startswith('#'):
            i += 1
            continue

        # Try to parse PDG ID line
        parts = line.split('#')[0].strip().split()
        if len(parts) == 1:
            try:
                pid = int(parts[0])
            except ValueError:
                i += 1
                continue

            i += 1
            # Next line should be number of channels
            while i < len(lines):
                cline = lines[i].strip()
                if cline and not cline.startswith('#'):
                    break
                i += 1
            if i >= len(lines):
                break

            n_channels = int(cline.split('#')[0].strip().split()[0])
            i += 1

            channels = []
            for _ in range(n_channels):
                while i < len(lines):
                    dline = lines[i].strip()
                    if dline and not dline.startswith('#'):
                        break
                    i += 1
                if i >= len(lines):
                    break

                dparts = dline.split('#')[0].strip().split()
                br = float(dparts[0])
                daughters = [int(d) for d in dparts[1:]]
                comment = ''
                if '#' in dline:
                    comment = '#' + dline.split('#', 1)[1]
                channels.append({'br': br, 'daughters': daughters, 'comment': comment})
                i += 1

            decays[pid] = channels
        else:
            i += 1

    return decays


# =====================================================================
# Compute decay threshold
# =====================================================================
def compute_threshold(pid, decay_channels, particles_by_id):
    """Compute the lightest decay threshold for a particle"""
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
# Sorting function for particle list output
# =====================================================================
def particle_sort_key(pid, pdata):
    """Sort key: light mesons, strange mesons, N, Delta, Lambda, Sigma, Xi, Omega"""
    B = pdata['B']
    S = pdata['S']
    Q = pdata['Q']
    mass = pdata['mass']
    I = pdata['I']

    if B == 0 and S == 0:
        return (0, mass, pid)  # Light unflavored mesons
    elif B == 0 and S != 0:
        return (1, mass, pid)  # Strange mesons
    elif B == 1 and S == 0:
        return (2, mass, pid)  # N and Delta baryons (interleaved by mass)
    elif B == 1 and S == -1 and I == 0:
        return (3, mass, pid)  # Lambda
    elif B == 1 and S == -1 and I > 0:
        return (4, mass, pid)  # Sigma
    elif B == 1 and S == -2:
        return (5, mass, pid)  # Xi
    elif B == 1 and S == -3:
        return (6, mass, pid)  # Omega
    else:
        return (7, mass, pid)


# =====================================================================
# Main generation function
# =====================================================================
def generate():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(script_dir)

    # Input: existing decays file for decay channels
    existing_decays_file = os.path.join(parent_dir, 'decays.dat')

    # Output files
    output_list = os.path.join(script_dir, 'list_from_pdfs.dat')
    output_decays = os.path.join(script_dir, 'decays_from_pdfs.dat')

    # Build particles_by_id lookup
    particles_by_id = dict(PARTICLES)

    # Parse existing decay channels
    print("Parsing existing decays.dat...")
    decays = parse_decays_file(existing_decays_file)
    print(f"  Found {len(decays)} decay entries")

    # Sort particles for output
    sorted_pids = sorted(PARTICLES.keys(), key=lambda pid: particle_sort_key(pid, PARTICLES[pid]))

    # Compute thresholds
    print("Computing thresholds...")
    thresholds = {}
    for pid in sorted_pids:
        pdata = PARTICLES[pid]
        if pdata['stable'] == 1:
            # For stable particles, threshold = mass of lightest decay
            if pid in decays:
                thresholds[pid] = compute_threshold(pid, decays[pid], particles_by_id)
            else:
                thresholds[pid] = 0
        else:
            if pid in decays:
                thresholds[pid] = compute_threshold(pid, decays[pid], particles_by_id)
            else:
                thresholds[pid] = 0

    # =====================================================================
    # Write list_from_pdfs.dat
    # =====================================================================
    print(f"Writing {output_list}...")
    with open(output_list, 'w') as f:
        f.write("# PDG2025 list containing strange and non-strange hadrons only\n")
        f.write("#         pdgid                name         stable      mass[GeV]"
                "     degeneracy     statistics              B              Q"
                "              S              C            |S|            |C|"
                "     width[GeV] threshold[GeV]\n")

        for pid in sorted_pids:
            p = PARTICLES[pid]
            # Compute degeneracy
            J = p['J']
            degeneracy = int(2 * J + 1)
            # Special override for f0(500)
            if 'degeneracy_override' in p:
                degeneracy = p['degeneracy_override']

            threshold = thresholds.get(pid, 0)

            line = (f"{pid:>15} {p['name']:>20} {p['stable']:>13} "
                    f"{p['mass']:>14g} {degeneracy:>14} {p['statistics']:>14} "
                    f"{p['B']:>14} {p['Q']:>14} {p['S']:>14} "
                    f"{p.get('charm', 0):>14} {p['absS']:>14g} {p['absC']:>14g} "
                    f"{p['width']:>14g} {threshold:>14g}")
            f.write(line + "\n")

    print(f"  Written {len(sorted_pids)} particles")

    # =====================================================================
    # Write decays_from_pdfs.dat
    # =====================================================================
    print(f"Writing {output_decays}...")
    with open(output_decays, 'w') as f:
        f.write("# the list of decays\n")
        f.write("# each entry consists of the following:\n")
        f.write("# a line with the pdgid of decaying particle\n")
        f.write("# a line with the number of decay channels\n")
        f.write("# for each channel a line containing whitespace-separated values of the channel branching ratio and pdg ids of the daughter products\n")
        f.write("# everything after the # symbol is treated as a comment and ignored\n")
        f.write("# decays of antiparticles are not listed but generated from the listed decays of particles\n")
        f.write("\n")

        written_pids = set()
        n_written = 0

        # Write decays in particle list order
        for pid in sorted_pids:
            if pid in decays and pid not in written_pids:
                channels = decays[pid]
                p = PARTICLES[pid]

                f.write(f"{pid:<20}                     # {p['name']}\n")
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
                n_written += 1

        # Also write extra decay entries not in our particle list but present in existing decays
        extra_count = 0
        for pid, channels in decays.items():
            if pid not in written_pids:
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

    print(f"  Written {n_written} hadron decay entries + {extra_count} extra (charm/nuclear) entries")

    # =====================================================================
    # Cross-check against existing list.dat
    # =====================================================================
    existing_list = os.path.join(parent_dir, 'list.dat')
    if os.path.exists(existing_list):
        print("\n" + "=" * 70)
        print("CROSS-CHECK AGAINST EXISTING list.dat")
        print("=" * 70)

        # Parse existing list.dat
        existing_particles = {}
        with open(existing_list, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 14:
                    epid = int(parts[0])
                    existing_particles[epid] = {
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

        print(f"\nExisting list: {len(existing_particles)} particles")
        print(f"New list:      {len(PARTICLES)} particles")

        # Check for missing/extra particles
        existing_ids = set(existing_particles.keys())
        new_ids = set(PARTICLES.keys())

        missing = existing_ids - new_ids
        extra = new_ids - existing_ids

        if missing:
            print(f"\nMISSING from new list ({len(missing)}):")
            for pid in sorted(missing):
                print(f"  {pid}: {existing_particles[pid]['name']}")

        if extra:
            print(f"\nEXTRA in new list ({len(extra)}):")
            for pid in sorted(extra):
                print(f"  {pid}: {PARTICLES[pid]['name']}")

        # Compare properties for common particles
        common = existing_ids & new_ids
        mass_diffs = []
        width_diffs = []
        quantum_diffs = []

        for pid in sorted(common):
            ep = existing_particles[pid]
            np_ = PARTICLES[pid]

            # Check mass
            if abs(ep['mass'] - np_['mass']) > 1e-6:
                mass_diffs.append((pid, ep['name'], ep['mass'], np_['mass']))

            # Check width
            if ep['width'] > 0 and np_['width'] > 0:
                if abs(ep['width'] - np_['width']) / max(ep['width'], np_['width']) > 0.01:
                    width_diffs.append((pid, ep['name'], ep['width'], np_['width']))
            elif abs(ep['width'] - np_['width']) > 1e-18:
                width_diffs.append((pid, ep['name'], ep['width'], np_['width']))

            # Check quantum numbers
            J_new = np_['J']
            deg_new = int(2 * J_new + 1)
            if 'degeneracy_override' in np_:
                deg_new = np_['degeneracy_override']
            if ep['degeneracy'] != deg_new:
                quantum_diffs.append((pid, ep['name'], 'degeneracy', ep['degeneracy'], deg_new))
            if ep['B'] != np_['B']:
                quantum_diffs.append((pid, ep['name'], 'B', ep['B'], np_['B']))
            if ep['Q'] != np_['Q']:
                quantum_diffs.append((pid, ep['name'], 'Q', ep['Q'], np_['Q']))
            if ep['S'] != np_['S']:
                quantum_diffs.append((pid, ep['name'], 'S', ep['S'], np_['S']))

        if mass_diffs:
            print(f"\nMass differences ({len(mass_diffs)}):")
            for pid, name, old_m, new_m in mass_diffs[:20]:
                print(f"  {name:>20} (ID {pid}): {old_m:.6f} -> {new_m:.6f}  (diff={new_m-old_m:+.6f})")
            if len(mass_diffs) > 20:
                print(f"  ... and {len(mass_diffs)-20} more")
        else:
            print("\nNo mass differences found!")

        if width_diffs:
            print(f"\nWidth differences ({len(width_diffs)}):")
            for pid, name, old_w, new_w in width_diffs[:20]:
                print(f"  {name:>20} (ID {pid}): {old_w:.6g} -> {new_w:.6g}")
            if len(width_diffs) > 20:
                print(f"  ... and {len(width_diffs)-20} more")
        else:
            print("\nNo width differences found!")

        if quantum_diffs:
            print(f"\nQuantum number differences ({len(quantum_diffs)}):")
            for pid, name, qn, old_v, new_v in quantum_diffs[:20]:
                print(f"  {name:>20} (ID {pid}): {qn} {old_v} -> {new_v}")
        else:
            print("\nNo quantum number differences found!")

        print("\n" + "=" * 70)
        print("Cross-check complete.")
        print("=" * 70)


if __name__ == '__main__':
    generate()
