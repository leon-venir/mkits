import numpy as np

"""

"""


"""
:constants
"""
cons_hbar = 1.0545718e-34  # m^2kg/s
cons_planck = 6.62607015e-34 # J/Hz
cons_kb = 1.38064852e-23  # m^2kg/(Ks^2)
cons_pi = 3.14159265359
cons_emass = 9.10938356e-31  # electron mass kg
cons_echarge = 1.60217662e-19  # electron charge



"""
:units convertor
"""
uc_ev2j = 1.602176634e-19  # electonic volt to Joule
uc_ry2j = 2.179874099E-18  # Ry
uc_gpa2nm2 = 1e9  # GPa to N/m^2
uc_bohr2ang = 0.529177 # bohr to angstrom
uc_ang2m = 1e-10
uc_ry2ev = 13.6056980659 # Rydberg to electonic volt
uc_d2a = np.pi/180 # degree to radian
uc_atomsphere2nm2 = 101325 


"""
:database 
"""


# :database bond length
bond_length = {
    "Sr-O": 2.6,
    "Ti-O": 1.93
}
rev_bond = {}
for key in bond_length.keys():
    atom1, atom2 = key.split("-")
    rev_pair = "%s-%s" %(atom2, atom1)
    if rev_pair not in bond_length.keys():
        rev_bond[rev_pair] = bond_length[key]
bond_length.update(rev_bond)
del(rev_bond)


# :database WIEN2K tempelates =================================== #
local_rot_matrix = """LOCAL ROT MATRIX:    1.0000000 0.0000000 0.0000000                             
                     0.0000000 1.0000000 0.0000000                             
                     0.0000000 0.0000000 1.0000000                             
"""


# :database quantum: VASP input tempelates ====================== #
qe_control_block = {
    "qeblock": "control",
    "restart_mode": '"from_scratch"',
    "nstep": "100",
    "outdir": '"./outdir"',
    "max_seconds": "4.d+4",
    "etot_conv_thr": "1.d-4",
    "pseudo_dir": '"./pseudo"'
}


qe_system_block = {
    "qeblock": "system",
    "ecutwfc":  "80",
}


qe_electrons_block = {
    "qeblock": "electrons",
    "electron_maxstep": "200",
    "conv_thr": "1.d-8",
    "mixing_beta": "0.7",
    "mixing_mode": '"TF" ! not converging? try "TF" for highly homogeneous, "local-TF" for highly inhomogeneous system'
}

qe_ions_block = {
    "qeblock": "ions",
    "ion_dynamics": '"bfgs"'
}

qe_cell_block = {
    "qeblock": "cell",
    "cell_dynamics": '"bfgs"',
    "press": "0.d0",
    "press_conv_thr": "0.01d0",
    "cell_factor": "2.0",
    "cell_dofree": '"all"'
}

qe_dos_block = {
    "qeblock": "dos",
    "outdir": './outdir',
    "fildos": 'dos.dat',
    "emin": "-10",
    "emax": "35"
}

qe_control_key = [
    "calculation", "title", "verbosity", "restart_mode", 
    "wf_collect", "nstep", "iprint", "tstress", "tprnfor", "dt", 
    "outdir", "wfcdir", "prefix", "lkpoint_dir", "max_seconds", 
    "etot_conv_thr", "forc_conv_thr", "disk_io", "pseudo_dir", 
    "tefield", "dipfield", "lelfield", "nberrycyc", "lorbm", 
    "lberry", "gdir", "nppstr", "gate", "lfcp", "trism"
]

qe_system_key = [
    "ibrav", "celldm", "A", "B", "C", "cosAB", "cosAC", "cosBC", 
    "nat", "ntyp", "nbnd", "tot_charge", "starting_charge", 
    "tot_magnetization", "starting_magnetization", "ecutwfc", 
    "ecutrho", "ecutfock", "nr1", "nr2", "nr3", "nr1s", "nr2s", 
    "nr3s", "nosym", "nosym_evc", "noinv", "no_t_rev", 
    "force_symmorphic", "use_all_frac", "occupations", 
    "one_atom_occupations", "starting_spin_angle", "degauss", 
    "smearing", "nspin", "noncolin", "ecfixed", "qcutz", 
    "q2sigma", "input_dft", "ace", "exx_fraction", 
    "screening_parameter", "exxdiv_treatment", 
    "x_gamma_extrapolation", "ecutvcut", "nqx1", "nqx2", "nqx3", 
    "localization_thr", "Hubbard_occ", "Hubbard_alpha", 
    "Hubbard_beta", "starting_ns_eigenvalue", "dmft", 
    "dmft_prefix", "ensemble_energies", "edir", "emaxpos", 
    "eopreg", "eamp", "angle1", "angle2", "lforcet", 
    "constrained_magnetization", "fixed_magnetization", "lambda", "report", "lspinorb", "assume_isolated", "esm_bc", "esm_w", "esm_efield", "esm_nfit", "lgcscf", "gcscf_mu", "gcscf_conv_thr", "gcscf_beta", "vdw_corr", "london", "london_s6", "london_c6", "london_rvdw", "london_rcut", "dftd3_version", "dftd3_threebody", "ts_vdw_econv_thr", "ts_vdw_isolated", "xdm", "xdm_a1", "xdm_a2", "space_group", "uniqueb", "origin_choice", "rhombohedral", "zgate", "relaxz", "block", "block_1", "block_2", "block_height"
]

qe_electrons_key = [
    "electron_maxstep", "scf_must_converge", "conv_thr", 
    "adaptive_thr", "conv_thr_init", "conv_thr_multi", 
    "mixing_mode", "mixing_beta", "mixing_ndim", 
    "mixing_fixed_ns", "diagonalization", "diago_thr_init", 
    "diago_cg_maxiter", "diago_ppcg_maxiter", 
    "diago_david_ndim", "diago_rmm_ndim", "diago_rmm_conv", 
    "diago_gs_nblock", "diago_full_acc", "efield", 
    "efield_cart", "efield_phase", "startingpot", 
    "startingwfc", "tqr", "real_space"
]

qe_ions_key = ["ion_positions", "ion_velocities", "ion_dynamics", 
               "pot_extrapolation", "wfc_extrapolation", "remove_rigid_rot", 
               "ion_temperature", "tempw", "tolp", "delta_t", "nraise", 
               "refold_pos", "upscale", "bfgs_ndim", "trust_radius_max", 
               "trust_radius_min", "trust_radius_ini", "w_1", "w_2", 
               "fire_alpha_init", "fire_falpha", "fire_nmin", "fire_f_inc", 
               "fire_f_dec", "fire_dtmax"]

qe_cell_key = [
    "cell_dynamics", "press", "wmass", "cell_factor", "press_conv_thr", 
    "cell_dofree"
]

qe_rism_key = [
    "nsolv", "closure", "tempv", "ecutsolv", "solute_lj", "solute_epsilon", 
    "solute_sigma", "starting1d", "starting3d", "smear1d", "smear3d", 
    "rism1d_maxstep", "rism3d_maxstep", "rism1d_conv_thr", "rism3d_conv_thr", 
    "mdiis1d_size", "mdiis3d_size", "mdiis1d_step", "mdiis3d_step", 
    "rism1d_bond_width", "rism1d_dielectric", "rism1d_molesize", 
    "rism1d_nproc", "rism3d_conv_level", "rism3d_planar_average", 
    "laue_nfit", "laue_expand_right", "laue_expand_left", 
    "laue_starting_right", "laue_starting_left", "laue_buffer_right", 
    "laue_buffer_left", "laue_both_hands", "laue_wall", "laue_wall_z", 
    "laue_wall_rho", "laue_wall_epsilon", "laue_wall_sigma", "laue_wall_lj6"
]


# :database vasp: VASP input tag ==============================================
incar_tag = [
    "ENCUT", "PREC", "ALGO", "ISMEAR", "SIGMA", "ISTART", "NELM", "LWAVE", 
    "LREAL", "AMIX", "NFREE", "BMIX", "NCORE", "NELECT", "IOPTCELL", "LAECHG",
    "LELF", "TIME", "LCHARG", "AMIN",
    # dos band
    "ICHARG", "EMIN", "EMAX", "NEDOS", "LORBIT",
    # DFT+U
    "LDAU", "LDAUTYPE", "LDAUL", "LDAUU", "LDAUJ", "LDAUPRINT", 
    # spin
    "ISPIN", "MAGMOM", 
    # so
    "LNONCOLLINEAR", "VOSKOWN", "LSORBIT", "SAXIS", "NBANDS",
    # opt
    "EDIFF", "EDIFFG", "NSW", "ISIF", "IBRION",
    # md
    "TEBEG", "TEEND", "SMASS",
    # sym
    "ISYM",
    # dipole
    "LDIPOL", "IDIPOL",
    # work function
    "LVTOT", "LVHAR",
    # vdw
    "IVDW",
    # optic
    "LOPTICS", "CSHIFT"
]


# :database vasp: VASP input tempelates =======================================
incar_glob = {
    "PREC": "Normal",
    "EDIFF": "1e-6",
    "NCORE": "16",
    "NELM": "1000",
    "NELMIN": "5",
    "LMAXMIX": "2",
    "AMIX": "0.4 # for metals smaller 0.02, or try other value 0.1",
    "BMIX": "1.0 # try 0.5"
}
incar_opt = {
    "NELM": "60",
    "IBRION": "2",
    "ISIF": "3",
    "NSW": "500",
    "EDIFFG": "-0.05",
    "POTIM": "0.5",
    "LWAVE": ".FALSE.",
    "LCHARG": ".FALSE."
}
incar_scf = {
    "IBRION": "-1",
    "NSW": "0"
}
incar_metal = {
    "ISMEAR": "2",
    "SIGMA" : "0.2"
}
incar_semi = {
    "ISMEAR": "0",
    "SIGMA" : "0.05"
}
incar_dos = {
    "ICHARG": "11", 
    "LORBIT": "11", 
    "IBRION": "-1", 
    "NSW": "0", 
    "ISMEAR": "-5", 
    "SIGMA": "0.05", 
    "EMIN": "-6",
    "EMAX": "6",
    "NEDOS": "2500", 
} # "NBANDS": "set_your_own_value"
incar_band = {
    "ICHARG": "11", 
    "IBRION": "-1", 
    "NSW": "0", 
    "ISMEAR": "0", 
    "SIGMA": "0.05", 
} # "NBANDS": "set_your_own_value"
incar_optic = {
    "ICHARG": "11", 
    "LORBIT": "11", 
    "IBRION": "-1", 
    "NSW": "0", 
    "ISMEAR": "-5", 
    "SIGMA": "0.05", 
    "EMIN": "-6",
    "EMAX": "6",
    "NEDOS": "2500", 
}
incar_freq = {
    "IBRION": "5", 
    "NSW": "500", 
    "ISMEAR": "0", 
    "SIGMA": "0.05",
    "ISIF": "2",
    "POTIM": "0.02",
    "NFREE": "2",
    "ALGO": "Fast"
} # "NBANDS": "set_your_own_value"
incar_nvt = {
    "IBRION": "0",
    "MDALGO": "2 # 1 Andersen 2 Nose-Hoover 3 Langevin 4 Multiple Andersen",
    "ISIF": "2",
    "SMASS": "1.0",
    "ISYM": "0",
    "ALGO": "FAST",
    "TEBEG": "300",
    "TEEND": "300",
    "POTIM": "1.0",
    "LREAL": "Auto",
    "NSW": "2000",
    "NELMIN": "4",
    "LWAVE": ".FALSE.",
    "LCHARG" : ".FALSE."
}
incar_nve = {
    # To calculate an NVE_ensemble we recommend the user to use MDALGO=1 and ANDERSEN_PROB=0.0.
}
incar_ml_heat = {
    "ML_LMLFF": ".TRUE.",
    " ": ""
}

incar_prop = {
    "xanes": {
        "#XANES": "",
        "ALGO": "FAST",
        "ISMEAR": "0",
        "SIGMA": "0.05",
        "ICORELEVEL": "2",
        "CLNT": "xxx",
        "CLN": "1",
        "CLL": "0",
        "CLZ": "1.0",
        "CH_LSPEC": ".TRUE.",
        "CH_SIGMA": "0.3",
        "CH_NEDOS": "1000",
        "LREAL": "Auto"
    },
    "dfpt": {
        "#DFPT": "input for DFPT",
        "PREC": "Accurate",
        "ENCUT": "450",
        "EDIFF": "1e-8",
        "IBRION": "-1",
        "NSW": "0",
        "NELM": "200",
        "ISMEAR": "0",
        "SIGMA": "0.01",
        "IALGO": "38",
        "LREAL": ".FALSE." 
    },
    "born": {
        "#BORN": "effective charge with DFPT",
        "PREC": "Accurate",
        "IBRION": "-1",
        "ISMEAR": "0",
        "SIGMA": "0.01",
        "EDIFF": "1.E-8",
        "LREAL": ".FALSE.",
        "LEPSILON": ".TRUE."
    }
}


incar_functionals = {
    "pbelda": {
        "#GGA_or_LDA": " "
    },
    "pbesol": {
        "#pbesol": " ",
        "GGA": "PS"
    },
    "rev-vdW-DF2": {
        "#vdw-rev-vdW-DF2": " ",
        "GGA": "MK",
        "LUSE_VDW": ".TRUE.",
        "AGGAC": "0.0",
        "LASPH": ".TRUE.",
        "PARAM1": "0.1234",
        "PARAM2": "0.711357",
        "ZAB_VDW": "-1.8867"
    },
    "optB88": {
        "#vdw-optB88": " ",
        "GGA": "BO",
        "PARAM1": "0.1833333333",
        "PARAM2": "0.22",
        "LUSE_VDW": ".TRUE.",
        "AGGAC": "0.0",
        "LASPH": ".TRUE."
    },
    "optPBE": {
        "#vdw-optPBE": " ",
        "GGA": "OR",
        "LUSE_VDW": ".TRUE.",
        "AGGAC": "0.0",
        "LASPH": ".TRUE."
    },
    "hsesol": {
        "#HSEsol": " ",
        "LHFCALC": ".TRUE.",
        "GGA": "PS",
        "HFSCREEN": "0.2",
        "ALGO": "Damped",
        "AEXX": "0.25",
        "AGGAX": "0.75",
        "AGGAC": "1.0",
        "ALDAC": "1.0"
    },
    "hse06" : {
        "#HSE06": " ",
        "LHFCALC": ".TRUE.",
        "GGA": "PE",
        "HFSCREEN": "0.2",
        "ALGO": "All",
        "AEXX": "0.25",
        "AGGAX": "0.75",
        "AGGAC": "1.0",
        "ALDAC": "1.0"
    }
}



"""
molecule data
atomic number, x, y, z in angstrom
"""
mol_h2o = np.array(
    [
        [8, 0, 0, 0],
        [1, 0, 0, 0],
        [1, 0, 0, 0]
    ]
)




"""
: Recommend PBE
"""
vasp_upf_recommend=[
    "H", "He", "Li_sv", "Be", "B", "C", "N", "OF", "Ne", "Na_pv", "Mg", "Al", 
    "Si", "P", "S", "Cl", "Ar", "K_sv", "Ca_sv", "Sc_sv", "Ti_sv", 
    "V_sv", "Cr_pv", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga_d", 
    "Ge_d", "As", "Se", "Br", "Kr", "Rb_sv", "Sr_sv", "Y_sv", 
    "Zr_sv", "Nb_sv", "Mo_sv", "Tc_pv", "Ru_pv", "Rh_pv", "Pd", 
    "Ag", "Cd", "In_d", "Sn_d", "Sb", "Te", "I", "Xe", "Cs_sv", 
    "Ba_sv", "La", "Ce", "Pr_3", "Nd_3", "Pm_3", "Sm_3", "Eu_2", 
    "Gd_3", "Tb_3", "Dy_3", "Ho_3", "Er_3", "Tm_3", "Yb_2", "Lu_3", 
    "Hf_pv", "Ta_pv", "W_sv", "Re", "Os", "Ir", "Pt", "Au", "Hg", 
    "Tl_d", "Pb_d", "Bi_d", "Po_d", "At", "Rn", "Fr_sv", "Ra_sv", 
    "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm"
    ]

atom_data = [ 
    #  0   1     2            3             4
    [  0, "N",   "Name",      "atomic mass", "magnetic", "mag_diret"], # 0
    [  1, "H",   "Hydrogen",  1.00794, 0], # 1
    [  2, "He",  "Helium",    4.002602, 0], # 2
    [  3, "Li",  "Lithium",   6.941, 0], # 3
    [  4, "Be",  "Beryllium", 9.012182], # 4
    [  5, "B",   "Boron",     10.811], # 5
    [  6, "C",   "Carbon",    12.0107, 1], # 6
    [  7, "N",   "Nitrogen",  14.0067, 1], # 7
    [  8, "O",   "Oxygen",    15.9994, 0.6], # 8
    [  9, "F",   "Fluorine",  18.9984032, 0.6], # 9
    [ 10, "Ne",  "Neon",      20.1797], # 10
    [ 11, "Na",  "Sodium",    22.98976928], # 11
    [ 12, "Mg",  "Magnesium", 24.3050], # 12
    [ 13, "Al",  "Aluminium", 26.9815386], # 13
    [ 14, "Si",  "Silicon",   28.0855], # 14
    [ 15, "P",   "Phosphorus", 30.973762], # 15
    [ 16, "S",   "Sulfur",    32.065, 0], # 16
    [ 17, "Cl",  "Chlorine", 35.453, 0], # 17
    [ 18, "Ar",  "Argon", 39.948], # 18
    [ 19, "K",   "Potassium", 39.0983], # 19
    [ 20, "Ca",  "Calcium", 40.078], # 20
    [ 21, "Sc",  "Scandium", 44.955912], # 21
    [ 22, "Ti",  "Titanium", 47.867, 0.6], # 22
    [ 23, "V",   "Vanadium", 50.9415, 1], # 23
    [ 24, "Cr",  "Chromium", 51.9961, 3], # 24
    [ 25, "Mn",  "Manganese", 54.938045, 3], # 25
    [ 26, "Fe",  "Iron", 55.845, 5], # 26
    [ 27, "Co",  "Cobalt", 58.933195, 3], # 27
    [ 28, "Ni",  "Nickel", 58.6934, 2], # 28
    [ 29, "Cu",  "Copper", 63.546, 1], # 29
    [ 30, "Zn",  "Zinc", 65.38], # 30
    [ 31, "Ga",  "Gallium", 69.723], # 31
    [ 32, "Ge",  "Germanium", 72.64], # 32
    [ 33, "As",  "Arsenic", 74.92160], # 33
    [ 34, "Se",  "Selenium", 78.96, 0], # 34
    [ 35, "Br",  "Bromine", 79.904], # 35
    [ 36, "Kr",  "Krypton", 83.798], # 36
    [ 37, "Rb",  "Rubidium", 85.4678], # 37
    [ 38, "Sr",  "Strontium", 87.62], # 38
    [ 39, "Y",   "Yttrium", 88.90585], # 39
    [ 40, "Zr",  "Zirconium", 91.224], # 40
    [ 41, "Nb",  "Niobium", 92.90638], # 41
    [ 42, "Mo",  "Molybdenum", 95.96, 2], # 42
    [ 43, "Tc",  "Technetium", None], # 43
    [ 44, "Ru",  "Ruthenium", 101.07, 2], # 44
    [ 45, "Rh",  "Rhodium", 102.90550], # 45
    [ 46, "Pd",  "Palladium", 106.42], # 46
    [ 47, "Ag",  "Silver", 107.8682], # 47
    [ 48, "Cd",  "Cadmium", 112.411], # 48
    [ 49, "In",  "Indium", 114.818], # 49
    [ 50, "Sn",  "Tin", 118.710], # 50
    [ 51, "Sb",  "Antimony", 121.760], # 51
    [ 52, "Te",  "Tellurium", 127.60, 0], # 52
    [ 53, "I",   "Iodine", 126.90447], # 53
    [ 54, "Xe",  "Xenon", 131.293], # 54
    [ 55, "Cs",  "Caesium", 132.9054519], # 55
    [ 56, "Ba",  "Barium", 137.327], # 56
    [ 57, "La",  "Lanthanum", 138.90547], # 57
    [ 58, "Ce",  "Cerium", 140.116], # 58
    [ 59, "Pr",  "Praseodymium", 140.90765], # 59
    [ 60, "Nd",  "Neodymium", 144.242], # 60
    [ 61, "Pm",  "Promethium", None], # 61
    [ 62, "Sm",  "Samarium", 150.36], # 62
    [ 63, "Eu",  "Europium", 151.964], # 63
    [ 64, "Gd",  "Gadolinium", 157.25], # 64
    [ 65, "Tb",  "Terbium", 158.92535], # 65
    [ 66, "Dy",  "Dysprosium", 162.500], # 66
    [ 67, "Ho",  "Holmium", 164.93032], # 67
    [ 68, "Er",  "Erbium", 167.259], # 68
    [ 69, "Tm",  "Thulium", 168.93421], # 69
    [ 70, "Yb",  "Ytterbium", 173.054], # 70
    [ 71, "Lu",  "Lutetium", 174.9668], # 71
    [ 72, "Hf",  "Hafnium", 178.49], # 72
    [ 73, "Ta",  "Tantalum", 180.94788], # 73
    [ 74, "W",   "Tungsten", 183.84, 2], # 74
    [ 75, "Re",  "Rhenium", 186.207], # 75
    [ 76, "Os",  "Osmium", 190.23], # 76
    [ 77, "Ir",  "Iridium", 192.217], # 77
    [ 78, "Pt",  "Platinum", 195.084], # 78
    [ 79, "Au",  "Gold", 196.966569], # 79
    [ 80, "Hg",  "Mercury", 200.59], # 80
    [ 81, "Tl",  "Thallium", 204.3833], # 81
    [ 82, "Pb",  "Lead", 207.2], # 82
    [ 83, "Bi",  "Bismuth", 208.98040], # 83
    [ 84, "Po",  "Polonium", None], # 84
    [ 85, "At",  "Astatine", None], # 85
    [ 86, "Rn",  "Radon", None], # 86
    [ 87, "Fr",  "Francium", None], # 87
    [ 88, "Ra",  "Radium", None], # 88
    [ 89, "Ac",  "Actinium", None], # 89
    [ 90, "Th",  "Thorium", 232.03806], # 90
    [ 91, "Pa",  "Protactinium", 231.03588], # 91
    [ 92, "U",   "Uranium", 238.02891], # 92
    [ 93, "Np",  "Neptunium", None], # 93
    [ 94, "Pu",  "Plutonium", None], # 94
    [ 95, "Am",  "Americium", None], # 95
    [ 96, "Cm",  "Curium", None], # 96
    [ 97, "Bk",  "Berkelium", None], # 97
    [ 98, "Cf",  "Californium", None], # 98
    [ 99, "Es",  "Einsteinium", None], # 99
    [100, "Fm",  "Fermium", None], # 100
    [101, "Md",  "Mendelevium", None], # 101
    [102, "No",  "Nobelium", None], # 102
    [103, "Lr",  "Lawrencium", None], # 103
    [104, "Rf",  "Rutherfordium", None], # 104
    [105, "Db",  "Dubnium", None], # 105
    [106, "Sg",  "Seaborgium", None], # 106
    [107, "Bh",  "Bohrium", None], # 107
    [108, "Hs",  "Hassium", None], # 108
    [109, "Mt",  "Meitnerium", None], # 109
    [110, "Ds",  "Darmstadtium", None], # 110
    [111, "Rg",  "Roentgenium", None], # 111
    [112, "Cn",  "Copernicium", None], # 112
    [113, "Uut", "Ununtrium", None], # 113
    [114, "Uuq", "Ununquadium", None], # 114
    [115, "Uup", "Ununpentium", None], # 115
    [116, "Uuh", "Ununhexium", None], # 116
    [117, "Uus", "Ununseptium", None], # 117
    [118, "Uuo", "Ununoctium", None], # 118
    ]

symbol_map = {
    "H":1,
    "He":2,
    "Li":3,
    "Be":4,
    "B":5,
    "C":6,
    "N":7,
    "O":8,
    "F":9,
    "Ne":10,
    "Na":11,
    "Mg":12,
    "Al":13,
    "Si":14,
    "P":15,
    "S":16,
    "Cl":17,
    "Ar":18,
    "K":19,
    "Ca":20,
    "Sc":21,
    "Ti":22,
    "V":23,
    "Cr":24,
    "Mn":25,
    "Fe":26,
    "Co":27,
    "Ni":28,
    "Cu":29,
    "Zn":30,
    "Ga":31,
    "Ge":32,
    "As":33,
    "Se":34,
    "Br":35,
    "Kr":36,
    "Rb":37,
    "Sr":38,
    "Y":39,
    "Zr":40,
    "Nb":41,
    "Mo":42,
    "Tc":43,
    "Ru":44,
    "Rh":45,
    "Pd":46,
    "Ag":47,
    "Cd":48,
    "In":49,
    "Sn":50,
    "Sb":51,
    "Te":52,
    "I":53,
    "Xe":54,
    "Cs":55,
    "Ba":56,
    "La":57,
    "Ce":58,
    "Pr":59,
    "Nd":60,
    "Pm":61,
    "Sm":62,
    "Eu":63,
    "Gd":64,
    "Tb":65,
    "Dy":66,
    "Ho":67,
    "Er":68,
    "Tm":69,
    "Yb":70,
    "Lu":71,
    "Hf":72,
    "Ta":73,
    "W":74,
    "Re":75,
    "Os":76,
    "Ir":77,
    "Pt":78,
    "Au":79,
    "Hg":80,
    "Tl":81,
    "Pb":82,
    "Bi":83,
    "Po":84,
    "At":85,
    "Rn":86,
    "Fr":87,
    "Ra":88,
    "Ac":89,
    "Th":90,
    "Pa":91,
    "U":92,
    "Np":93,
    "Pu":94,
    "Am":95,
    "Cm":96,
    "Bk":97,
    "Cf":98,
    "Es":99,
    "Fm":100,
    "Md":101,
    "No":102,
    "Lr":103,
    "Rf":104,
    "Db":105,
    "Sg":106,
    "Bh":107,
    "Hs":108,
    "Mt":109,
    "Ds":110,
    "Rg":111,
    "Cn":112,
    "Uut":113,
    "Uuq":114,
    "Uup":115,
    "Uuh":116,
    "Uus":117,
    "Uuo":118,
    }

""" :database encut: recommended kinetic cut off of each element. """
kinetic_cutoff = [
    [  0, "vasp",   "qe", "", ""], # 0
    [ 14,    450], # 14 Si
    [ 22,    450], # 22 Ti
    [ 23,    500], # 
    [ 32,    450], # Ge
    [ 51,    500], # Sb
    [ 52,    500]  # Te
]

""" :database smearing: recommended smearing, ie sigma in VASP, """
smearing = [
    [  0, "vasp",   "qe", ""], # 0
    [ 14,   0.05,      0], # Si
    [ 22,   0.05], # Ti
    [ 23,   0.05], # 
    [ 32,   0.05], # Ge
    [ 51,   0.05], # Sb
    [ 52,   0.05]  # Te
]

""" :database VASP recommended PAW """
vasp_paw = [
    [  0, "PBE", "val_num",  "LDA"],
    [ 14, "Si", 6, "Si", 6],
    [ 32, "Ge_h", 14, "Ge_h", 14],
    [ 51, "Sb", 5, "Sb", 5],
    [ 52, "Te", 6, "Te", 6]
]


""" :database kpath_sg: recommended kpath for each space group. """
kpath_sg = {
    164: np.array([[0,0,0], [1,2,3]]),
    164+500: ["Gamma", ""]
}





""" :script """
""" :script gnu2dlines: """
gnu2dline = """set terminal pngcairo enhanced font "Times,24" lw 5 size 1000, 900
set output '%s'
unset key
set title '%s'
set xlabel "%s" font "Times,28"
set ylabel "%s" font "Times,28"
set grid
set colorsequence default # [default, classic]

set style line 1 dt 1 lt 1 lw 0.8 pt 6 ps 3
set style line 2 dt 1 lt 2 lw 0.8 pt 7 ps 3
set style line 3 dt 1 lt 3 lw 0.8 pt 8 ps 3
set style line 4 dt 1 lt 4 lw 0.8 pt 9 ps 3
set style line 5 dt 1 lt 5 lw 0.8 pt 10 ps 3
set style line 6 dt 1 lt 6 lw 0.8 pt 11 ps 3
set style line 7 dt 1 lt 7 lw 0.8 pt 12 ps 3

set style line 11 dt 3 lt 1 lw 0.3 pt 6 ps 3
set style line 12 dt 3 lt 2 lw 0.3 pt 7 ps 3
set style line 13 dt 3 lt 3 lw 0.3 pt 8 ps 3
set style line 14 dt 3 lt 4 lw 0.3 pt 9 ps 3
set style line 15 dt 3 lt 5 lw 0.3 pt 10 ps 3
set style line 16 dt 3 lt 6 lw 0.3 pt 11 ps 3
set style line 17 dt 3 lt 7 lw 0.3 pt 12 ps 3

plot "%s" with linespoints ls 1 notitle
"""


""" :gnuplot band structures plottig """
gnubandcar = """set terminal pngcairo enhanced font "Times,20" lw 5 size 900, 900
set output '%s.png'
unset key

x_max = %.8f
x_min = 0
y_max = 6
y_min = -6

set title 'Band structures'

set ylabel "Energy (eV)" font "Times,24"
# replace with your own labels
set xtics (%s) 

set ytics # font "Times-Roman, 16"
set xrange [x_min:x_max]
set yrange [y_min:y_max]

set key left top # font "Times-Roman,16" 

set colorsequence default
set style line 3 dt 1 lt 3 lw 0.8 pt 8 ps 3
set style line 16 dt 3 lt 6 lw 0.3 pt 11 ps 3
set style line 17 dt 3 lt 7 lw 0.3 pt 12 ps 3

high_points = "%s"

set for [x_coord in high_points ] arrow from x_coord,y_min to x_coord,y_max nohead ls 17 front
set arrow from x_min,0 to x_max,0 nohead ls 16 front

plot "BANDCAR" using ($1):($2 - %.8f) with line ls 3 notitle"""



""" :script extract_vasp_xanes_outcar: """
extract_xanes_vasp = """#!/bin/bash

parallel=-1
normal=-1
all=-1
tauc=-1
trace=-1
while [[ $# -gt 0 ]]
do
   key="$1"
   case $key in
      -parallel) parallel=0
      ;;
      -normal) normal=0
      ;;
      -trace) trace=0
      ;;
      -tauc) tauc=0
      ;;
   esac
   shift
done


cat > helpscript.perl  <<EOF
#!/bin/perl

use strict;
use warnings;
my \$mode=shift;

while(<>)
{
   chomp;
   if(\$_ =~ /frequency dependent IMAGINARY DIELECTRIC FUNCTION/)
   {
      \$_=<>;
      \$_=<>;
      while (<>)
      {
         my \$sum=0;
         if (\$_ !~ /[0-9]/) {last;}
         chomp;
         \$_=~s/^/ /;
         my @help=split(/[\t,\s]+/);
         if (\$help[2]=~/NaN/||\$help[3]=~/NaN/||\$help[4]=~/NaN/) {next;}
         if (\$help[5]=~/NaN/||\$help[6]=~/NaN/||\$help[4]=~/NaN/) {next;}
         if (\$mode==0) {\$sum=\$help[2]+\$help[3]+\$help[4]+\$help[5]+\$help[6]+\$help[7];}
         if (\$mode==1) {\$sum=\$help[4];}
         if (\$mode==2) {\$sum=\$help[2]+\$help[3];}
         if (\$mode==3) {\$sum=\$help[2]+\$help[3]+\$help[4];}
         if (\$mode==4) {\$sum=(\$help[1]*\$help[1]*(\$help[2]+\$help[3]+\$help[4]+\$help[5]+\$help[6]+\$help[7]))**0.5;}
         if (\$mode==5) {\$sum=(\$help[1]*\$help[1]*(\$help[4]))**0.5;}
         if (\$mode==6) {\$sum=(\$help[1]*\$help[1]*(\$help[2]+\$help[3]))**0.5;}
         if (\$mode==7) {\$sum=(\$help[1]*\$help[1]*(\$help[2]+\$help[3]+\$help[4]))**0.5;}
         print \$help[1]," ",\$sum,"\n";
      }
   }
   last if eof;
}
EOF

if [[ $normal -eq 0 ]]; then
   if [[ $tauc -eq 0 ]]; then
      perl helpscript.perl 4 OUTCAR > xanes.dat
   else
      perl helpscript.perl 1 OUTCAR > xanes.dat
   fi
else
   if [[ $parallel -eq 0 ]]; then
      if [[ $tauc -eq 0 ]]; then
         perl helpscript.perl 5 OUTCAR > xanes.dat
      else
         perl helpscript.perl 2 OUTCAR > xanes.dat
      fi
   else
      if [[ $trace -eq 0 ]]; then
         if [[ $tauc -eq 0 ]]; then
            perl helpscript.perl 6 OUTCAR > xanes.dat
         else
            perl helpscript.perl 3 OUTCAR > xanes.dat
         fi
      else
         if [[ $tauc -eq 0 ]]; then
            perl helpscript.perl 7 OUTCAR > xanes.dat
         else
            perl helpscript.perl 0 OUTCAR > xanes.dat
         fi
      fi
   fi
fi
rm helpscript.perl
"""


""" :scipt monitor_cluster """
monitor_cluster = """import subprocess
import time
from datetime import datetime

"""


""" :script py_cif2vasp """
py_cif2vasp = """%s
import sys
import ase.io
import ase.io.vasp

ciffile = sys.argv[1]

if ciffile[-4:] != ".cif":
    print("ERROR")
else:
    atoms = ase.io.read(ciffile)
    #atoms.write(ciffile[:-4]+'.vasp', format = 'vasp')
    ase.io.vasp.write_vasp(ciffile[:-4]+'.vasp', atoms, direct=True)"""



"""
: Useful equations
"""
def o_chemical_potential(t, p, e_tot_o2, e_zpe_o2, ):
    """
    This is the equation to calculation the chemical potential for 
    oxygen, 

    Parameters
    ----------

    
    
    
    
    The $\mu_{\mathrm{O}}$ is defined as:

    $$
    \mu_{\mathrm{O}} (T,P) = \frac{1}{2} E_{\mathrm{O_2}}^{\mathrm{total}} +
    \frac{1}{2} E_{\mathrm{O_2}}^{\mathrm{ZPE}} \Delta \mu_{\mathrm{O}} (T,P)
    $$

    $$ \Delta \mu_{\mathrm{O}} (T,P) = -\frac{1}{2}k_\mathrm{B}T 
    \left\{\begin{matrix}
    \overset{\mathrm{Translational}}{\overbrace{\ln \left [ \left ( 
    \frac{2 \pi m}{h^2}\right )^{3/2} \frac{ \left ( k_\mathrm{b} T \right ) 
    ^{5/2} }{P} \right ]}} + \overset{\mathrm{Rotational}}{\overbrace{ 
    \ln \left (\frac{k_\mathrm{B} T}{ \sigma^{\mathrm{sym}} B_0 }\right ) }}
    -\end{matrix}\right. $$

    $$\left.\begin{matrix}
    \overset{\mathrm{Vibrational}}{\overbrace{\ln \left [ 1 - \exp \left ( \frac{\hbar \omega_0}{k_\mathrm{B} T} \right ) \right ]}} + \overset{\mathrm{Electronic}}{\overbrace{\ln \left ( I^{\mathrm{spin}} \right )}}\end{matrix}\right\} $$

    where $E_{\mathrm{O_2}}^{\mathrm{ZPE}}$ is the zero point energy, $\hbar$ is the reduced Planck constant, $\sigma^{\mathrm{sym}}$ is the symmetry number for the O2 molecule, $B_0 = \hbar / 2 I$ is the rotational constant, $I$ is the the moment of inertia, $\omega_0$ is the vibrational mode of the O2 molecule considering an harmonic oscillator, and $I^{\mathrm{spin}}$ is the electronic spin degeneracy of the ground state.
    """
    pass
