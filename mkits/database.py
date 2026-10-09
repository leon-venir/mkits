# -*- coding: utf-8 -*-
"""Shared constant tables used by mkits.structure and other mkits modules
(mkits.builder's adsorption molecule geometries included).

DFT-code-specific templates (INCAR/QE keyword tables, gnuplot scripts, ...)
live in their respective modules (vasp.py, qe.py, ...) instead.
"""

import numpy as np


# ============================================================== #
# unit conversion
# ============================================================== #
uc_bohr2ang = 0.529177  # bohr to angstrom
uc_ry2ev = 13.605693  # rydberg to eV
uc_ha2ev = 27.211386  # hartree to eV


# ============================================================== #
# periodic table
# ------------------------------------------------------------ #
# atom_data[Z] = [Z, symbol, name, atomic_mass, magmom]
#   atomic_mass: standard atomic weight in amu, None if unknown/unstable.
#   magmom: recommended initial magnetic moment (mu_B) used to seed
#           spin-polarized DFT calculations; 0.0 for elements without a
#           commonly used non-zero starting guess.
# Z = 0 is reserved as a placeholder for an unrecognized element ("X").
# ============================================================== #
atom_data = [
    [0, "X", "Unknown", None, 0.0],
    [1, "H", "Hydrogen", 1.00794, 0.0],
    [2, "He", "Helium", 4.002602, 0.0],
    [3, "Li", "Lithium", 6.941, 0.0],
    [4, "Be", "Beryllium", 9.012182, 0.0],
    [5, "B", "Boron", 10.811, 0.0],
    [6, "C", "Carbon", 12.0107, 1.0],
    [7, "N", "Nitrogen", 14.0067, 1.0],
    [8, "O", "Oxygen", 15.9994, 0.6],
    [9, "F", "Fluorine", 18.9984032, 0.6],
    [10, "Ne", "Neon", 20.1797, 0.0],
    [11, "Na", "Sodium", 22.98976928, 0.0],
    [12, "Mg", "Magnesium", 24.3050, 0.0],
    [13, "Al", "Aluminium", 26.9815386, 0.0],
    [14, "Si", "Silicon", 28.0855, 0.0],
    [15, "P", "Phosphorus", 30.973762, 0.0],
    [16, "S", "Sulfur", 32.065, 0.0],
    [17, "Cl", "Chlorine", 35.453, 0.0],
    [18, "Ar", "Argon", 39.948, 0.0],
    [19, "K", "Potassium", 39.0983, 0.0],
    [20, "Ca", "Calcium", 40.078, 0.0],
    [21, "Sc", "Scandium", 44.955912, 0.0],
    [22, "Ti", "Titanium", 47.867, 0.6],
    [23, "V", "Vanadium", 50.9415, 1.0],
    [24, "Cr", "Chromium", 51.9961, -2.0],
    [25, "Mn", "Manganese", 54.938045, 3.0],
    [26, "Fe", "Iron", 55.845, 3.0],
    [27, "Co", "Cobalt", 58.933195, 3.0],
    [28, "Ni", "Nickel", 58.6934, 1.0],
    [29, "Cu", "Copper", 63.546, 0.0],
    [30, "Zn", "Zinc", 65.38, 0.0],
    [31, "Ga", "Gallium", 69.723, 0.0],
    [32, "Ge", "Germanium", 72.64, 0.0],
    [33, "As", "Arsenic", 74.92160, 0.0],
    [34, "Se", "Selenium", 78.96, 0.0],
    [35, "Br", "Bromine", 79.904, 0.0],
    [36, "Kr", "Krypton", 83.798, 0.0],
    [37, "Rb", "Rubidium", 85.4678, 0.0],
    [38, "Sr", "Strontium", 87.62, 0.0],
    [39, "Y", "Yttrium", 88.90585, 0.0],
    [40, "Zr", "Zirconium", 91.224, 0.0],
    [41, "Nb", "Niobium", 92.90638, 0.0],
    [42, "Mo", "Molybdenum", 95.96, 2.0],
    [43, "Tc", "Technetium", None, 0.0],
    [44, "Ru", "Ruthenium", 101.07, 2.0],
    [45, "Rh", "Rhodium", 102.90550, 0.0],
    [46, "Pd", "Palladium", 106.42, 0.0],
    [47, "Ag", "Silver", 107.8682, 0.0],
    [48, "Cd", "Cadmium", 112.411, 0.0],
    [49, "In", "Indium", 114.818, 0.0],
    [50, "Sn", "Tin", 118.710, 0.0],
    [51, "Sb", "Antimony", 121.760, 0.0],
    [52, "Te", "Tellurium", 127.60, 0.0],
    [53, "I", "Iodine", 126.90447, 0.0],
    [54, "Xe", "Xenon", 131.293, 0.0],
    [55, "Cs", "Caesium", 132.9054519, 0.0],
    [56, "Ba", "Barium", 137.327, 0.0],
    [57, "La", "Lanthanum", 138.90547, 0.0],
    [58, "Ce", "Cerium", 140.116, 0.0],
    [59, "Pr", "Praseodymium", 140.90765, 0.0],
    [60, "Nd", "Neodymium", 144.242, 0.0],
    [61, "Pm", "Promethium", None, 0.0],
    [62, "Sm", "Samarium", 150.36, 0.0],
    [63, "Eu", "Europium", 151.964, 0.0],
    [64, "Gd", "Gadolinium", 157.25, 0.0],
    [65, "Tb", "Terbium", 158.92535, 0.0],
    [66, "Dy", "Dysprosium", 162.500, 0.0],
    [67, "Ho", "Holmium", 164.93032, 0.0],
    [68, "Er", "Erbium", 167.259, 0.0],
    [69, "Tm", "Thulium", 168.93421, 0.0],
    [70, "Yb", "Ytterbium", 173.054, 0.0],
    [71, "Lu", "Lutetium", 174.9668, 0.0],
    [72, "Hf", "Hafnium", 178.49, 0.0],
    [73, "Ta", "Tantalum", 180.94788, 0.0],
    [74, "W", "Tungsten", 183.84, 2.0],
    [75, "Re", "Rhenium", 186.207, 0.0],
    [76, "Os", "Osmium", 190.23, 0.0],
    [77, "Ir", "Iridium", 192.217, 0.0],
    [78, "Pt", "Platinum", 195.084, 0.0],
    [79, "Au", "Gold", 196.966569, 0.0],
    [80, "Hg", "Mercury", 200.59, 0.0],
    [81, "Tl", "Thallium", 204.3833, 0.0],
    [82, "Pb", "Lead", 207.2, 0.0],
    [83, "Bi", "Bismuth", 208.98040, 0.0],
    [84, "Po", "Polonium", None, 0.0],
    [85, "At", "Astatine", None, 0.0],
    [86, "Rn", "Radon", None, 0.0],
    [87, "Fr", "Francium", None, 0.0],
    [88, "Ra", "Radium", None, 0.0],
    [89, "Ac", "Actinium", None, 0.0],
    [90, "Th", "Thorium", 232.03806, 0.0],
    [91, "Pa", "Protactinium", 231.03588, 0.0],
    [92, "U", "Uranium", 238.02891, 0.0],
    [93, "Np", "Neptunium", None, 0.0],
    [94, "Pu", "Plutonium", None, 0.0],
    [95, "Am", "Americium", None, 0.0],
    [96, "Cm", "Curium", None, 0.0],
    [97, "Bk", "Berkelium", None, 0.0],
    [98, "Cf", "Californium", None, 0.0],
    [99, "Es", "Einsteinium", None, 0.0],
    [100, "Fm", "Fermium", None, 0.0],
    [101, "Md", "Mendelevium", None, 0.0],
    [102, "No", "Nobelium", None, 0.0],
    [103, "Lr", "Lawrencium", None, 0.0],
    [104, "Rf", "Rutherfordium", None, 0.0],
    [105, "Db", "Dubnium", None, 0.0],
    [106, "Sg", "Seaborgium", None, 0.0],
    [107, "Bh", "Bohrium", None, 0.0],
    [108, "Hs", "Hassium", None, 0.0],
    [109, "Mt", "Meitnerium", None, 0.0],
    [110, "Ds", "Darmstadtium", None, 0.0],
    [111, "Rg", "Roentgenium", None, 0.0],
    [112, "Cn", "Copernicium", None, 0.0],
    [113, "Uut", "Ununtrium", None, 0.0],
    [114, "Uuq", "Ununquadium", None, 0.0],
    [115, "Uup", "Ununpentium", None, 0.0],
    [116, "Uuh", "Ununhexium", None, 0.0],
    [117, "Uus", "Ununseptium", None, 0.0],
    [118, "Uuo", "Ununoctium", None, 0.0],
]

symbol_map = {row[1]: row[0] for row in atom_data if row[0] != 0}


# ============================================================== #
# mkits.builder.adsorpt_molecue: adsorbate molecule geometries
# ------------------------------------------------------------ #
# Row 0 of each table is a placeholder "central position" (unused by
# adsorpt_molecue, kept for layout parity with atom_data); real atoms are
# rows 1: [atomic_number, x, y, z] in angstrom, already centered so that
# distance/orientation offsets can be applied directly.
# ============================================================== #
mol_h2o = np.array([
    [0, 0, 0, 0],
    [8, 0, 0, 0],
    [1, 0, -0.76913, 0.59479],
    [1, 0, 0.76913, 0.59479],
])

mol_h2 = np.array([
    [0, 0, 0, 0],
    [1, 0, -0.35913, 0],
    [1, 0, 0.35913, 0],
])

mol_oh = np.array([
    [0, 0, 0, 0],
    [8, 0, 0, 0],
    [1, 0, 0.76913, 0.59479],
])

mol_ooh = np.array([
    [0, 0, 0, 0],
    [8, 0, 0, 0],
    [8, 0.86900, 1.07800, 0.52270],
    [1, 1.67520, 1.39700, 1.21580],
])

mol_h = np.array([
    [0, 0, 0, 0],
    [1, 0, 0, 0],
])

mol_o = np.array([
    [0, 0, 0, 0],
    [8, 0, 0, 0],
])

mol_o2 = np.array([
    [0, 0, 0, 0],
    [8, 0, -0.61913, 0],
    [8, 0, 0.61913, 0],
])


# ============================================================== #
# WIEN2k templates
# ------------------------------------------------------------ #
# Identity local rotation matrix, used when writing case.struct files.
# Symmetry-adapted local rotation matrices are left to WIEN2k's own
# symmetry utilities (x nn / symmetso) run after the file is generated.
# ============================================================== #
local_rot_matrix = """LOCAL ROT MATRIX:    1.0000000 0.0000000 0.0000000
                     0.0000000 1.0000000 0.0000000
                     0.0000000 0.0000000 1.0000000
"""


# ============================================================== #
# high-symmetry k-path placeholders
# ------------------------------------------------------------ #
# NOTE: these tables are a first, self-built approximation intended to
# unblock band-structure workflows. They only distinguish crystal systems
# (3D) / 2D Bravais lattice types (2D) and do NOT yet account for lattice
# centering (P/I/F/C). Replace with literature-verified paths once
# available (see mkits.structure.seekpath).
# ============================================================== #
kpath_3d = {
    "cubic": {
        "points": {
            "GAMMA": [0.0, 0.0, 0.0],
            "X": [0.5, 0.0, 0.0],
            "M": [0.5, 0.5, 0.0],
            "R": [0.5, 0.5, 0.5],
        },
        "path": [["GAMMA", "X", "M", "GAMMA", "R", "X"]],
    }, # DOI: 
    "tetragonal": {
        "points": {
            "GAMMA": [0.0, 0.0, 0.0],
            "X": [0.5, 0.0, 0.0],
            "M": [0.5, 0.5, 0.0],
            "Z": [0.0, 0.0, 0.5],
            "R": [0.5, 0.0, 0.5],
            "A": [0.5, 0.5, 0.5],
        },
        "path": [["GAMMA", "X", "M", "GAMMA", "Z", "R", "A", "Z"]],
    },
    "orthorhombic": {
        "points": {
            "GAMMA": [0.0, 0.0, 0.0],
            "X": [0.5, 0.0, 0.0],
            "Y": [0.0, 0.5, 0.0],
            "Z": [0.0, 0.0, 0.5],
            "S": [0.5, 0.5, 0.0],
            "T": [0.0, 0.5, 0.5],
            "U": [0.5, 0.0, 0.5],
            "R": [0.5, 0.5, 0.5],
        },
        "path": [["GAMMA", "X", "S", "Y", "GAMMA", "Z", "U", "R", "T", "Z"]],
    },
    "hexagonal": {
        "points": {
            "GAMMA": [0.0, 0.0, 0.0],
            "M": [0.5, 0.0, 0.0],
            "K": [1.0 / 3.0, 1.0 / 3.0, 0.0],
            "A": [0.0, 0.0, 0.5],
            "L": [0.5, 0.0, 0.5],
            "H": [1.0 / 3.0, 1.0 / 3.0, 0.5],
        },
        "path": [["GAMMA", "M", "K", "GAMMA", "A", "L", "H", "A"]],
    },
    "trigonal": {
        "points": {
            "GAMMA": [0.0, 0.0, 0.0],
            "L": [0.5, 0.0, 0.0],
            "X": [0.5, 0.5, 0.0],
            "Z": [0.5, 0.5, 0.5],
        },
        "path": [["GAMMA", "L", "X", "GAMMA", "Z"]],
    },
    "monoclinic": {
        "points": {
            "GAMMA": [0.0, 0.0, 0.0],
            "Y": [0.0, 0.5, 0.0],
            "H": [0.0, 0.5, 0.5],
            "C": [0.0, 0.0, 0.5],
        },
        "path": [["GAMMA", "Y", "H", "C", "GAMMA"]],
    },
    "triclinic": {
        "points": {
            "GAMMA": [0.0, 0.0, 0.0],
            "X": [0.5, 0.0, 0.0],
            "Y": [0.0, 0.5, 0.0],
            "Z": [0.0, 0.0, 0.5],
        },
        "path": [["X", "GAMMA", "Y", "GAMMA", "Z"]],
    },
}

kpath_2d = {
    "square": {
        "points": {
            "GAMMA": [0.0, 0.0, 0.0],
            "X": [0.5, 0.0, 0.0],
            "M": [0.5, 0.5, 0.0],
        },
        "path": [["GAMMA", "X", "M", "GAMMA"]],
    },
    "rectangular": {
        "points": {
            "GAMMA": [0.0, 0.0, 0.0],
            "X": [0.5, 0.0, 0.0],
            "S": [0.5, 0.5, 0.0],
            "Y": [0.0, 0.5, 0.0],
        },
        "path": [["GAMMA", "X", "S", "Y", "GAMMA"]],
    },
    "hexagonal": {
        "points": {
            "GAMMA": [0.0, 0.0, 0.0],
            "K": [1.0 / 3.0, 1.0 / 3.0, 0.0],
            "M": [0.5, 0.0, 0.0],
        },
        "path": [["GAMMA", "M", "K", "GAMMA"]],
    },
    "oblique": {
        "points": {
            "GAMMA": [0.0, 0.0, 0.0],
            "X": [0.5, 0.0, 0.0],
            "Y": [0.0, 0.5, 0.0],
        },
        "path": [["X", "GAMMA", "Y"]],
    },
}
