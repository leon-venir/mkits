# -*- coding: utf-8 -*-
"""Quantum ESPRESSO input generation: pw.x and companion executables
(dos.x, bands.x, ph.x, neb.x, epw.x, pwcond.x, hp.x)."""

import os
import re
import shutil

import numpy as np

from mkits import database
from mkits import functions
from mkits import structure


# ================================================================== #
# namelist schema: one authoritative "valid keys" set per namelist.
# QE-specific, so it lives here rather than in mkits.database.
# ================================================================== #
_CONTROL_KEYS = {
    "calculation", "title", "verbosity", "restart_mode",
    "wf_collect", "nstep", "iprint", "tstress", "tprnfor", "dt",
    "outdir", "wfcdir", "prefix", "lkpoint_dir", "max_seconds",
    "etot_conv_thr", "forc_conv_thr", "disk_io", "pseudo_dir",
    "tefield", "dipfield", "lelfield", "nberrycyc", "lorbm",
    "lberry", "gdir", "nppstr", "gate", "lfcp", "trism"
}

_SYSTEM_KEYS = {
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
    "localization_thr", "lda_plus_u", "lda_plus_u_kind", "Hubbard_U",
    "Hubbard_V", "Hubbard_J", "Hubbard_J0", "Hubbard_occ", "Hubbard_alpha",
    "Hubbard_beta", "starting_ns_eigenvalue", "dmft",
    "dmft_prefix", "ensemble_energies", "edir", "emaxpos",
    "eopreg", "eamp", "angle1", "angle2", "lforcet",
    "constrained_magnetization", "fixed_magnetization", "lambda", "report",
    "lspinorb", "assume_isolated", "esm_bc", "esm_w", "esm_efield",
    "esm_nfit", "lgcscf", "gcscf_mu", "gcscf_conv_thr", "gcscf_beta",
    "vdw_corr", "london", "london_s6", "london_c6", "london_rvdw",
    "london_rcut", "dftd3_version", "dftd3_threebody", "ts_vdw_econv_thr",
    "ts_vdw_isolated", "xdm", "xdm_a1", "xdm_a2", "space_group", "uniqueb",
    "origin_choice", "rhombohedral", "zgate", "relaxz", "block", "block_1",
    "block_2", "block_height"
}

_ELECTRONS_KEYS = {
    "electron_maxstep", "scf_must_converge", "conv_thr",
    "adaptive_thr", "conv_thr_init", "conv_thr_multi",
    "mixing_mode", "mixing_beta", "mixing_ndim",
    "mixing_fixed_ns", "diagonalization", "diago_thr_init",
    "diago_cg_maxiter", "diago_ppcg_maxiter",
    "diago_david_ndim", "diago_rmm_ndim", "diago_rmm_conv",
    "diago_gs_nblock", "diago_full_acc", "efield",
    "efield_cart", "efield_phase", "startingpot",
    "startingwfc", "tqr", "real_space"
}

_IONS_KEYS = {
    "ion_positions", "ion_velocities", "ion_dynamics",
    "pot_extrapolation", "wfc_extrapolation", "remove_rigid_rot",
    "ion_temperature", "tempw", "tolp", "delta_t", "nraise",
    "refold_pos", "upscale", "bfgs_ndim", "trust_radius_max",
    "trust_radius_min", "trust_radius_ini", "w_1", "w_2",
    "fire_alpha_init", "fire_falpha", "fire_nmin", "fire_f_inc",
    "fire_f_dec", "fire_dtmax"
}

_CELL_KEYS = {
    "cell_dynamics", "press", "wmass", "cell_factor", "press_conv_thr",
    "cell_dofree"
}

# --- companion executables: one namelist each, so no per-calculation-type
# routing table (_CALC_NAMELISTS) is needed for them -- each generator
# below just opens a QEInput with its own single namelist name. ---

_DOS_KEYS = {
    "prefix", "outdir", "ngauss", "degauss", "Emin", "Emax", "DeltaE",
    "fildos", "bz_sum"
}

_BANDS_KEYS = {
    "prefix", "outdir", "filband", "spin_component", "lsigma",
    "lp", "filp", "lsym", "no_overlap", "plot_2d", "firstk", "lastk"
}

# Core, high-confidence ph.x (DFPT phonon) keywords only; ph.x supports
# several more specialized options (eg Allen-Heine-Cardona electron-phonon
# corrections, q-point plotting) not included here.
_INPUTPH_KEYS = {
    "amass", "outdir", "prefix", "fildyn", "fildrho", "fildvscf",
    "epsil", "trans", "ldisp", "nq1", "nq2", "nq3",
    "tr2_ph", "alpha_mix", "niter_ph", "recover",
    "start_q", "last_q", "start_irr", "last_irr",
    "search_sym", "reduce_io", "zeu",
    "electron_phonon", "el_ph_sigma", "el_ph_nsigma",
    "verbosity", "low_directory_check"
}

# MEDIUM CONFIDENCE: hp.x (Hubbard parameters via density-functional
# perturbation theory, the HP package -- merged into QE relatively
# recently) keyword set reconstructed from its documented &INPUTHP
# namelist, not verified against a live QE install. Only the core
# single-pass, full-q-mesh workflow is implemented by hp_geninput itself;
# the split-by-q-point/split-by-atom workflow (start_q/last_q,
# perturb_only_atom, a separate compute_hp=.true./sum_pertq=.true.
# collection run) is only reachable via `params` overrides, not built
# into hp_geninput -- verify against your installed QE version's HP/Doc
# before relying on this for production use.
_INPUTHP_KEYS = {
    "prefix", "outdir", "iverbosity", "max_seconds",
    "nq1", "nq2", "nq3", "skip_equivalence_q",
    "determine_num_pert_only", "determine_q_mesh_only", "find_atpert",
    "docc_thr", "skip_type", "equiv_type", "perturb_only_atom",
    "start_q", "last_q", "sum_pertq", "compute_hp",
    "conv_thr_chi", "thresh_init", "ethr_nscf", "niter_max",
    "alpha_mix", "nmix", "num_neigh", "lmin", "rmax"
}

# NEB's &PATH namelist; NEB does not use _NAMELIST_KEYS/QEInput at all (see
# qe_nebinput) because its file format nests a whole pw.x-like "engine"
# input plus one ATOMIC_POSITIONS block per image inside BEGIN/END
# markers, which does not fit the flat namelist+card document model.
_PATH_KEYS = {
    "restart_mode", "string_method", "nstep_path", "opt_scheme",
    "CI_scheme", "num_of_images", "k_max", "k_min", "path_thr",
    "ds", "first_last_opt", "use_masses", "use_freezing", "temp_req"
}

# epw.x has a very large, highly workflow-specific parameter set (Wannier90
# disentanglement, electron-phonon interpolation, transport, ...); only the
# keywords needed for the common "Wannierize + interpolate EPC" workflow
# are whitelisted here -- extend this set for anything beyond that.
_INPUTEPW_KEYS = {
    "prefix", "outdir", "dvscf_dir",
    "elph", "epbwrite", "epbread", "epwwrite", "epwread",
    "wannierize", "num_iter", "dis_win_min", "dis_win_max",
    "dis_froz_min", "dis_froz_max", "proj", "wdata", "bands_skipped",
    "nbndsub", "vme",
    "nk1", "nk2", "nk3", "nq1", "nq2", "nq3",
    "nkf1", "nkf2", "nkf3", "nqf1", "nqf2", "nqf3",
    "fsthick", "degaussw", "eps_acoustic",
    "etf_mem", "system_2d", "lpolar",
    "filkf", "filqf", "band_plot",
    # transport / mobility (Boltzmann transport equation) keywords -- LOWER
    # CONFIDENCE than the Wannierize+interpolate subset above, since EPW's
    # transport module output format has changed across releases; verify
    # against your installed EPW version.
    "scattering", "scattering_serta", "int_mob", "iterative_bte",
    "ncarrier", "mob_maxiter", "epmatkqread", "temps", "restart"
}

# LOW CONFIDENCE: pwcond.x is a niche, infrequently-used tool; this is a
# minimal skeleton of the keywords I am reasonably sure about, not a
# complete &INPUTCOND schema. It also does not set up the required
# separate lead/scattering-region scf calculations pwcond.x depends on.
# Verify carefully against the QE documentation before production use.
_INPUTCOND_KEYS = {
    "outdir", "prefixt", "prefixl", "prefixr", "prefixs", "tran_prefix",
    "save_file", "ikind", "iofspin", "energy0", "denergy", "nenergy",
    "ecut2d", "ewind", "epsproj", "orbj_in", "orbj_fin"
}

_NAMELIST_KEYS = {
    "control": _CONTROL_KEYS,
    "system": _SYSTEM_KEYS,
    "electrons": _ELECTRONS_KEYS,
    "ions": _IONS_KEYS,
    "cell": _CELL_KEYS,
    "dos": _DOS_KEYS,
    "bands": _BANDS_KEYS,
    "inputph": _INPUTPH_KEYS,
    "inputhp": _INPUTHP_KEYS,
    "inputepw": _INPUTEPW_KEYS,
    "inputcond": _INPUTCOND_KEYS,
}


# ================================================================== #
# value formatting
# ================================================================== #
def _coerce_param_value(value):
    """
    Coerce a CLI-supplied string value (eg from functions.parser_inputpara)
    to its likely native Python type, so _format_namelist_value formats it
    correctly. Values already given as their native Python type (bool,
    int, float) pass straight through unchanged.
    """
    if not isinstance(value, str):
        return value
    _lower = value.strip().lower()
    if _lower in (".true.", "true"):
        return True
    if _lower in (".false.", "false"):
        return False
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        pass
    return value


def _format_namelist_value(value):
    """
    Format a native Python value as a QE Fortran-namelist literal:
    bool -> .TRUE./.FALSE., numbers -> their plain string form, a
    space/comma-separated list of numbers (eg epw.x's "temps" array) ->
    left unquoted, other strings -> single-quoted (unless already quoted,
    or already a bare logical literal).
    """
    if isinstance(value, bool):
        return ".TRUE." if value else ".FALSE."
    if isinstance(value, (int, float)):
        return str(value)
    _text = str(value).strip()
    if _text.upper() in (".TRUE.", ".FALSE."):
        return _text.upper()
    if _text[:1] in ("'", '"'):
        return _text
    _tokens = _text.replace(",", " ").split()
    try:
        [float(_t) for _t in _tokens]
        _is_numeric_list = bool(_tokens)
    except ValueError:
        _is_numeric_list = False
    if _is_numeric_list:
        return _text
    return "'%s'" % _text


# ================================================================== #
# unified namelist + card document
# ================================================================== #
class QEInput(object):
    """
    In-memory representation of a pw.x input file: an ordered set of
    namelists plus an ordered list of already-rendered card blocks.
    `qe_geninput` builds the defaults, feature toggles (add_hubbard_u,
    add_magnetism, ...) mutate it, and `write` serializes everything.
    """

    def __init__(self, namelist_order):
        self.namelist_order = list(namelist_order)
        self.namelists = {_name: {} for _name in self.namelist_order}
        self.cards = []

    def set(self, key, value, namelist=None):
        """
        Set one namelist parameter. If `namelist` is not given, the key is
        routed to whichever active namelist declares it; an unknown key,
        or a key that belongs to a namelist not part of this document,
        raises a clear error instead of being silently dropped.
        """
        value = _coerce_param_value(value)
        _base_key = key.split("(")[0]

        if namelist is not None:
            self.namelists[namelist][key] = value
            return

        for _name in self.namelist_order:
            if _base_key in _NAMELIST_KEYS[_name]:
                self.namelists[_name][key] = value
                return

        for _name in _NAMELIST_KEYS:
            if _base_key in _NAMELIST_KEYS[_name]:
                raise functions.MkitsError(
                    "QE keyword '%s' belongs to &%s, which is not part of this %s calculation."
                    % (key, _name.upper(), "+".join(self.namelist_order))
                )
        raise functions.MkitsError("Unknown QE namelist keyword: %s" % key)

    def update(self, params):
        """Route every key of a flat dict (eg from functions.parser_inputpara)."""
        for _key, _value in params.items():
            self.set(_key, _value)

    def add_card(self, lines):
        """Append a fully-rendered card block (its own header line included)."""
        self.cards.append(lines)

    def to_lines(self):
        """Serialize all namelists and cards into pw.x input file lines."""
        _lines = []
        for _name in self.namelist_order:
            _lines.append("&%s\n" % _name.upper())
            for _key, _value in self.namelists[_name].items():
                _lines.append("   %s = %s\n" % (_key, _format_namelist_value(_value)))
            _lines.append("/\n\n")
        for _card in self.cards:
            _lines += _card
            _lines.append("\n")
        return _lines

    def write(self, fpath, fname):
        """Write the assembled pw.x input file."""
        with open(os.path.join(fpath, fname), "w", newline="\n") as f:
            f.writelines(self.to_lines())


# ================================================================== #
# card builders
# ================================================================== #
def _atomic_species_card(atomic_symbols, upf_info):
    """Build the ATOMIC_SPECIES card lines."""
    _lines = ["ATOMIC_SPECIES\n"]
    for _symbol in atomic_symbols:
        _mass = database.atom_data[
            database.symbol_map[_symbol]
        ][3] or 0.0
        _lines.append(
            "%-4s %10.5f %s\n" % (
                _symbol, _mass, upf_info[_symbol]["upf"]
            )
        )
    return _lines


def _cell_parameters_card(struct_obj):
    """Build the CELL_PARAMETERS card lines (always angstrom)."""
    _lines = ["CELL_PARAMETERS angstrom\n"]
    for _row in struct_obj.lattice9:
        _lines.append(
            "%18.10f %18.10f %18.10f\n" % (_row[0], _row[1], _row[2])
        )
    return _lines


def _atomic_positions_card(struct_obj, frac=True, dyn=False):
    """Build the ATOMIC_POSITIONS card lines, optionally with if_pos flags."""
    _tag = "crystal" if frac else "angstrom"
    _lines = ["ATOMIC_POSITIONS %s\n" % _tag]
    _coord = struct_obj.position[1:, 1:4] if frac else struct_obj.position[1:, 4:7]
    for _i in range(struct_obj.total_atom):
        _z = int(struct_obj.position[_i + 1, 0])
        _symbol = database.atom_data[_z][1]
        _x, _y, _w = _coord[_i]
        if dyn:
            _fx, _fy, _fz = struct_obj.position[_i + 1, 7:10]
            _lines.append("%-4s %18.10f %18.10f %18.10f %d %d %d\n" % (
                _symbol, _x, _y, _w, int(_fx > 0.5), int(_fy > 0.5), int(_fz > 0.5)
            ))
        else:
            _lines.append("%-4s %18.10f %18.10f %18.10f\n" % (_symbol, _x, _y, _w))
    return _lines


def _kpoints_automatic_card(mesh):
    """Build a K_POINTS automatic card (Gamma-centered Monkhorst-Pack, unshifted)."""
    return ["K_POINTS automatic\n", "%d %d %d 0 0 0\n" % mesh]


def _kpoints_crystal_card(kpoints):
    """
    Build a K_POINTS crystal card: an explicit list of fractional
    (reciprocal-lattice) k-points, equal weight -- used for a
    non-self-consistent run over an arbitrary k-point set (eg
    structure.local_kbox's local effective-mass-tensor sampling) rather than
    a Monkhorst-Pack mesh or a high-symmetry path. Unlike VASP's
    equivalent explicit-list KPOINTS mode, QE's "crystal" unit is already
    fractional -- no unit conversion is needed on the way in (only on the
    way back out when reading pw.x's own tpiba-unit stdout, see
    read_pwout_bands).
    """
    _lines = ["K_POINTS crystal\n", "%d\n" % len(kpoints)]
    for _k in kpoints:
        _lines.append("%16.10f %16.10f %16.10f 1.0\n" % tuple(_k))
    return _lines


# ================================================================== #
# pseudopotentials (SSSP v2.0 only, for now)
# ================================================================== #
def _upf_valence_electrons(lines):
    """
    Read the number of valence electrons declared in a UPF file's header.
    Handles both the plain-text UPF v1 header ("<value>   Z valence") and
    the XML-based UPF v2 header ('z_valence="<value>"').
    """
    _z_val = 0.0
    for _line in lines:
        if "Z valence" in _line:
            _clean = _line.replace("=", " ").replace(",", " ").replace('"', " ")
            _z_val = float(_clean.split()[0])
        elif "z_valence" in _line:
            _clean = _line.replace("=", " ").replace(",", " ").replace('"', " ")
            _z_val = float(_clean.split()[-1])
    return _z_val


def _upf_is_uspp(lines):
    """
    Detect whether a UPF file is ultrasoft, needed to pick a safe default
    ecutrho/ecutwfc ratio. Checks the UPF v2 XML `is_ultrasoft="T"/"F"`
    attribute value (not just its presence, since the attribute *name*
    contains the substring "ultrasoft" regardless of its value); falls
    back to the "Ultrasoft"/"USPP" descriptive text used in the older
    plain-text UPF v1 header.
    """
    for _line in lines:
        _match = re.search(r'is_ultrasoft\s*=\s*"?\s*([TFtf])', _line)
        if _match:
            return _match.group(1).upper() == "T"
        if "Ultrasoft" in _line or "USPP" in _line:
            return True
    return False


def sssp_v2_resolve_upf(atomic_symbols, upfpath, dest_dir):
    """
    Resolve one UPF pseudopotential per element from a local SSSP v2.0
    library, copy each into `dest_dir`, and read its valence electron
    count and ultrasoft/PAW-vs-norm-conserving nature from the file
    content (not from the filename).

    Only the SSSP v2.0 "<Element>.<...>.upf" naming convention is
    supported (eg "F.us.pbesol.z_7.ld1.psl.v0.1.upf"): the element symbol
    is matched as the exact token before the first "."; other
    pseudopotential libraries are not handled yet.

    :param atomic_symbols: element symbols present in the structure
    :param upfpath: directory containing the SSSP v2.0 UPF files
    :param dest_dir: destination directory (eg "<wkdir>/pseudo")

    Return
    ------
    dict: {symbol: {"upf": filename, "valence": float, "is_uspp": bool}}
    """
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)

    _by_symbol = {}
    for _fname in os.listdir(upfpath):
        if not _fname.lower().endswith(".upf"):
            continue
        _symbol = _fname.split(".")[0]
        _by_symbol.setdefault(_symbol, _fname)

    _resolved = {}
    for _symbol in atomic_symbols:
        if _symbol not in _by_symbol:
            raise functions.MkitsError("No SSSP v2.0 UPF file found for element %s in %s." % (_symbol, upfpath))
        _fname = _by_symbol[_symbol]
        _src = os.path.join(upfpath, _fname)
        with open(_src, "r", errors="ignore") as f:
            _upflines = f.readlines()

        _is_uspp = _upf_is_uspp(_upflines)
        _valence = _upf_valence_electrons(_upflines)

        shutil.copy(_src, os.path.join(dest_dir, _fname))
        _resolved[_symbol] = {"upf": _fname, "valence": _valence, "is_uspp": _is_uspp}

    return _resolved


# ================================================================== #
# feature toggles
# ================================================================== #
def add_hubbard_u(qeinput, ggau, atomic_type):
    """
    Enable a simplified DFT+U treatment using the legacy
    lda_plus_u=.true. + Hubbard_U(i) syntax (supported by every QE
    version, unlike the newer HUBBARD card).

    :param qeinput: QEInput
    :param ggau: "Ti:4.2,Cu:6.2" -- U value (eV) per element symbol
    :param atomic_type: ordered list of atomic numbers matching ATOMIC_SPECIES
    """
    _u_values = functions.parser_inputpara(ggau)
    qeinput.set("lda_plus_u", True)
    for _i, _z in enumerate(atomic_type):
        _symbol = database.atom_data[_z][1]
        _u = float(_u_values[_symbol]) if _symbol in _u_values else 0.0
        qeinput.set("Hubbard_U(%d)" % (_i + 1), _u)


def add_magnetism(qeinput, atomic_type):
    """
    Enable spin polarization (nspin=2) and seed a starting_magnetization
    per species from mkits.database's default magnetic moments. This is a
    coarse initial guess to help the SCF cycle break spin symmetry, not a
    physically converged magnitude, and it is ferromagnetic-like (uniform
    sign per species); antiferromagnetic seeding needs per-atom overrides.

    :param qeinput: QEInput
    :param atomic_type: ordered list of atomic numbers matching ATOMIC_SPECIES
    """
    qeinput.set("nspin", 2)
    for _i, _z in enumerate(atomic_type):
        _seed = 0.4 if database.atom_data[_z][4] != 0 else 0.0
        qeinput.set("starting_magnetization(%d)" % (_i + 1), _seed)


# ================================================================== #
# qe_geninput
# ================================================================== #
_CALC_LABEL = {
    "opt": "relax",
    "vcopt": "vc-relax",
    "scf": "scf",
    "nscf": "nscf",
    "band": "bands",
}

_CALC_NAMELISTS = {
    "scf": ("control", "system", "electrons"),
    "nscf": ("control", "system", "electrons"),
    "band": ("control", "system", "electrons"),
    "opt": ("control", "system", "electrons", "ions"),
    "vcopt": ("control", "system", "electrons", "ions", "cell"),
}


def qe_geninput(
        calculation="scf",
        struct_inp="POSCAR",
        wpath="./",
        wname="none",
        upfpath="./",
        metal=True,
        mag=False,
        execode="mpirun -np $SLURM_NTASKS pw.x",
        params="gga=pbe",
        label=None,
        kpoints_override=None,
        **kwargs
):
    """
    Generate a single Quantum ESPRESSO pw.x calculation (input file + run
    script), written as "pw<label>.in" -- eg "pwscf.in", "pwband.in" --
    mirroring mkits.vasp's "INCAR_<dft>" convention closely enough that
    several calculation types can share one directory. Chaining several
    calculations (opt -> scf -> band -> dos) is left to mkits.workflow.

    :param calculation: one of opt, vcopt, scf, nscf, band -- the same
        vocabulary as mkits.vasp.vasp_gen_input's `dft` where the concept
        exists in both codes ("opt" writes calculation='relax', "vcopt"
        writes calculation='vc-relax', "band" writes calculation='bands')
    :param struct_inp: path to the input structure (any mkits.structure format)
    :param wpath: root working directory
    :param wname: calculation subdirectory name; defaults to `label`
    :param upfpath: directory containing SSSP v2.0 UPF files
    :param metal: metallic (Gaussian smearing) vs insulating (fixed occupations)
    :param mag: enable spin polarization with a heuristic starting guess
    :param execode: the mpirun/srun invocation, without -in/-out redirection
    :param params: "key=value,key=value" overrides; keys matching a QE
        namelist keyword are routed there (raising a clear error if none
        matches); "kspacing", "oddeven", "kfix", "kpoints_per_segment" are
        mkits-specific k-mesh controls, not QE keywords
    :param label: overrides the filename/directory role name (defaults to
        `calculation`) without changing the real QE calculation written
        inside the file -- eg mkits.workflow generates the DOS-densifying
        step as calculation="nscf", label="dos", so it is named
        "pwdos.in" while still containing calculation='nscf'
    :param kpoints_override: optional (n, 3) array of explicit fractional
        k-points; when given, this bypasses both structure.kmesh() and (for
        calculation="band") seekpath, writing this exact k-point set as a
        "K_POINTS crystal" card instead -- used by
        mkits.workflow.qe_workflow_effective_mass to sample a
        structure.local_kbox() neighborhood around a band extremum

    kwargs
    ------
    ggau: str
        "Ti:4.2,Cu:6.2" -- enable DFT+U for the given elements
    dynrange: str
        "xmin=..,xmax=..,fix=O" -- selective-dynamics-style atom constraints
    """
    if calculation not in _CALC_NAMELISTS:
        raise functions.MkitsError("Unknown calculation type: %s" % calculation)
    _label = label if label is not None else calculation

    _params = functions.parser_inputpara(params)
    _kspacing = float(_params.pop("kspacing", 0.25 if calculation == "scf" else 0.30))
    _oddeven = _params.pop("oddeven", "none")
    _kfix = [int(v) for v in _params.pop("kfix").split()] if "kfix" in _params else None
    _kpoints_per_segment = int(_params.pop("kpoints_per_segment", 20))

    _struct = structure.struct(struct_inp)

    _atomic_indices = [
        int(_) for _ in _struct.position[1:, 0]
    ]
    _atomic_type = sorted(set(_atomic_indices))
    _atomic_num = [
        _atomic_indices.count(_) for _ in _atomic_type
    ]
    _atomic_symbols = [
        database.atom_data[_][1] for _ in _atomic_type
    ]

    if wpath == "./":
        wpath = os.path.abspath("./")
    _wkdir = os.path.join(wpath, wname if wname != "none" else _label)
    os.makedirs(_wkdir, exist_ok=True)

    _upf_info = sssp_v2_resolve_upf(_atomic_symbols, upfpath, os.path.join(_wkdir, "pseudo"))
    _val_electron = sum(
        _upf_info[_atomic_symbols[_i]]["valence"] * _atomic_num[_i] for _i in range(len(_atomic_symbols))
    )
    _is_uspp = any(_upf_info[_s]["is_uspp"] for _s in _atomic_symbols)

    _dyn = False
    if "dynrange" in kwargs:
        _dyn = True
        functions.apply_dynrange(_struct, kwargs["dynrange"])

    _qe = QEInput(_CALC_NAMELISTS[calculation])

    # --- CONTROL ---
    _qe.set("calculation", _CALC_LABEL[calculation])
    _qe.set("restart_mode", "from_scratch")
    _qe.set("pseudo_dir", "./pseudo")
    _qe.set("outdir", "./outdir")
    _qe.set("verbosity", "low")
    if calculation in ("opt", "vcopt"):
        _qe.set("forc_conv_thr", 1.0e-4)

    # --- SYSTEM ---
    _qe.set("ibrav", 0)
    _qe.set("nat", _struct.total_atom)
    _qe.set("ntyp", len(_atomic_type))
    _qe.set("input_dft", _params.pop("gga", "pbe"))
    _ecutwfc = float(_params.pop("ecutwfc", 60.0))
    _qe.set("ecutwfc", _ecutwfc)
    _qe.set(
        "ecutrho", float(
            _params.pop("ecutrho", _ecutwfc * (
                10.0 if _is_uspp else 4.0)
            )
        )
    )

    if metal:
        _qe.set("occupations", "smearing")
        _qe.set("smearing", "gaussian")
        _qe.set("degauss", float(_params.pop("degauss", 0.01)))
        _qe.set("nbnd", int(_val_electron / 2.0 * float(_params.pop("nbndfactor", 1.3))))
    else:
        _qe.set("occupations", "fixed")
        if calculation in ("nscf", "band"):
            _qe.set("nbnd", int(_val_electron / 2.0 * float(_params.pop("nbndfactor", 1.5))))

    if mag:
        add_magnetism(_qe, _atomic_type)
    if "ggau" in kwargs:
        add_hubbard_u(_qe, kwargs["ggau"], _atomic_type)

    # any remaining user-supplied keys are routed directly to their namelist
    _qe.update(_params)

    # --- cards ---
    _qe.add_card(_atomic_species_card(_atomic_symbols, _upf_info))
    _qe.add_card(_cell_parameters_card(_struct))
    _qe.add_card(_atomic_positions_card(_struct, frac=True, dyn=_dyn))

    if kpoints_override is not None:
        _qe.add_card(_kpoints_crystal_card(kpoints_override))
    elif calculation == "band":
        _sp = structure.seekpath(struct_inp)
        _qe.add_card(_sp.write_kpath(code="qe", kpoints_per_segment=_kpoints_per_segment, write2file=False))
    else:
        _mesh = _struct.kmesh(kspacing=_kspacing, oddeven=_oddeven, kfix=_kfix)
        _qe.add_card(_kpoints_automatic_card(_mesh))

    _qe.write(_wkdir, "pw%s.in" % _label)

    _cmd = "# %s calculation\n" % _label
    _cmd += "%s -in pw%s.in > pw%s.out 2>&1\n" % (execode, _label, _label)
    functions.write_runsh(os.path.join(_wkdir, "run_pw%s.sh" % _label), _cmd)


# ================================================================== #
# companion executables: dos.x, bands.x, ph.x, neb.x, epw.x, pwcond.x
# ================================================================== #
def qepost_geninput(
        tool="dos",
        prefix="pwscf",
        outdir="../scf/outdir",
        wpath="./",
        wname="none",
        execode="mpirun -np $SLURM_NTASKS",
        params=""
):
    """
    Generate the input file (+ run script) for dos.x or bands.x, the QE
    post-processing tools that read a finished (n)scf calculation's save
    directory and produce plottable band/DOS data. Unlike VASP, where band
    and DOS data already sit in vasprun.xml after a single run, QE splits
    band-sorting (bands.x) and DOS smearing (dos.x) into separate
    post-processing steps because pw.x's own per-k-point band ordering is
    not continuous across k-points, and DOS needs its own broadening.

    Written as "post_<tool>.in" (eg "post_band.in", "post_dos.in") so it
    never collides with the "pw<label>.in" pw.x input it reads the output
    of (eg qe_geninput(calculation="band")'s "pwband.in") -- both can then
    share one directory.

    :param tool: "dos" or "band" (bands.x is still the executable name)
    :param prefix: must match the &CONTROL "prefix" of the finished (n)scf run
    :param outdir: the &CONTROL "outdir" of that (n)scf run
    :param params: "key=value,key=value" overrides for &DOS / &BANDS
    """
    if tool not in ("dos", "band"):
        raise functions.MkitsError("Unknown post-processing tool: %s (expected 'dos' or 'band')" % tool)

    _params = functions.parser_inputpara(params) if params else {}

    if wpath == "./":
        wpath = os.path.abspath("./")
    _wkdir = os.path.join(wpath, wname if wname != "none" else tool)
    os.makedirs(_wkdir, exist_ok=True)

    _namelist = "dos" if tool == "dos" else "bands"
    _qe = QEInput((_namelist,))
    _qe.set("prefix", prefix)
    _qe.set("outdir", outdir)

    if tool == "dos":
        _qe.set("fildos", _params.pop("fildos", "%s.dos" % prefix))
        _qe.set("degauss", float(_params.pop("degauss", 0.01)))
        _qe.set("DeltaE", float(_params.pop("DeltaE", 0.01)))
        _exe = "dos.x"
    else:
        _qe.set("filband", _params.pop("filband", "filband"))
        _qe.set("lsym", True)
        _exe = "bands.x"

    _qe.update(_params)
    _qe.write(_wkdir, "post_%s.in" % tool)

    _cmd = "# %s post-processing\n" % tool
    _cmd += "%s %s -in post_%s.in > post_%s.out 2>&1\n" % (execode, _exe, tool, tool)
    functions.write_runsh(os.path.join(_wkdir, "run_post_%s.sh" % tool), _cmd)


def qe_phinput(
        struct_inp="POSCAR",
        prefix="pwscf",
        outdir="../scf/outdir",
        wpath="./",
        wname="none",
        execode="mpirun -np $SLURM_NTASKS",
        qmesh=(1, 1, 1),
        params=""
):
    """
    Generate a ph.x (DFPT phonon) input file for a full q-point grid
    (ldisp=.true.). `struct_inp` must be the same structure used by the
    preceding scf run, so amass(i) lines up with that run's ATOMIC_SPECIES
    ordering (sorted by atomic number -- the same convention qe_geninput
    uses for ATOMIC_SPECIES).

    :param outdir: the &CONTROL "outdir" of the finished scf run
    :param qmesh: (nq1, nq2, nq3) q-point grid
    """
    _params = functions.parser_inputpara(params) if params else {}

    _struct = structure.struct(struct_inp)
    _atomic_indices = [int(_) for _ in _struct.position[1:, 0]]
    _atomic_type = sorted(set(_atomic_indices))

    if wpath == "./":
        wpath = os.path.abspath("./")
    _wkdir = os.path.join(wpath, wname if wname != "none" else "ph")
    os.makedirs(_wkdir, exist_ok=True)

    _qe = QEInput(("inputph",))
    _qe.set("prefix", prefix)
    _qe.set("outdir", outdir)
    _qe.set("fildyn", _params.pop("fildyn", "%s.dyn" % prefix))
    _qe.set("ldisp", True)
    _qe.set("nq1", qmesh[0])
    _qe.set("nq2", qmesh[1])
    _qe.set("nq3", qmesh[2])
    _qe.set("tr2_ph", float(_params.pop("tr2_ph", 1.0e-14)))
    for _i, _z in enumerate(_atomic_type):
        _qe.set("amass(%d)" % (_i + 1), database.atom_data[_z][3] or 0.0)
    _qe.update(_params)
    _qe.write(_wkdir, "ph.in")

    _cmd = "# ph.x phonon calculation\n"
    _cmd += "%s ph.x -in ph.in > ph.out 2>&1\n" % execode
    functions.write_runsh(os.path.join(_wkdir, "run_ph.sh"), _cmd)


def hp_geninput(
        prefix="pwscf",
        outdir="../scf/outdir",
        wpath="./",
        wname="none",
        execode="mpirun -np $SLURM_NTASKS",
        qmesh=(1, 1, 1),
        params=""
):
    """
    Generate an hp.x (Hubbard parameters via density-functional
    perturbation theory, the HP package) input file for a full q-point
    grid -- computes self-consistent Hubbard U (and, for an extended
    Hubbard scf, V) parameters in a single pass over `qmesh`. Splitting
    the calculation by q-point/atom for large systems (start_q/last_q,
    perturb_only_atom, plus a separate compute_hp=.true./sum_pertq=.true.
    collection run) is not built into this function -- reachable only via
    `params`, see _INPUTHP_KEYS's confidence caveat.

    The preceding scf run must already have DFT+U enabled (see
    qe_geninput's `ggau` kwarg / mkits.qe.add_hubbard_u): hp.x reads the
    Hubbard_U(i) values to perturb from that run's own save directory,
    not from anything generated here -- no `struct_inp` is needed.

    :param prefix: must match the &CONTROL "prefix" of the finished
        (DFT+U) scf run
    :param outdir: the &CONTROL "outdir" of that scf run
    :param qmesh: (nq1, nq2, nq3) q-point grid for the linear-response calculation
    :param params: "key=value,key=value" overrides for &INPUTHP
    """
    _params = functions.parser_inputpara(params) if params else {}

    if wpath == "./":
        wpath = os.path.abspath("./")
    _wkdir = os.path.join(wpath, wname if wname != "none" else "hp")
    os.makedirs(_wkdir, exist_ok=True)

    _qe = QEInput(("inputhp",))
    _qe.set("prefix", prefix)
    _qe.set("outdir", outdir)
    _qe.set("nq1", qmesh[0])
    _qe.set("nq2", qmesh[1])
    _qe.set("nq3", qmesh[2])
    _qe.update(_params)
    _qe.write(_wkdir, "hp.in")

    _cmd = "# hp.x Hubbard parameters calculation\n"
    _cmd += "%s hp.x -in hp.in > hp.out 2>&1\n" % execode
    functions.write_runsh(os.path.join(_wkdir, "run_hp.sh"), _cmd)


def _render_namelist(name, params):
    """
    Render a single "&NAME ... /" namelist block. Used by qe_nebinput,
    which does not fit the QEInput namelist+card document model (see
    _PATH_KEYS).
    """
    _lines = ["&%s\n" % name.upper()]
    for _key, _value in params.items():
        _lines.append("   %s = %s\n" % (_key, _format_namelist_value(_value)))
    _lines.append("/\n\n")
    return _lines


def qe_nebinput(
        first_struct,
        last_struct,
        num_of_images=7,
        prefix="pwscf",
        wpath="./",
        wname="neb",
        upfpath="./",
        metal=True,
        execode="mpirun -np $SLURM_NTASKS",
        params=""
):
    """
    Generate a climbing-image NEB input file. `first_struct`/`last_struct`
    become the first/last images; intermediate images are seeded by linear
    interpolation of fractional coordinates (a common but crude initial
    guess -- replace with a better interpolation, eg IDPP, for production
    runs where the linear path crosses atoms).

    CAVEAT: QE's NEB input uses a sectioned
    BEGIN/BEGIN_PATH_INPUT/BEGIN_ENGINE_INPUT/BEGIN_POSITIONS format that
    has had some syntax churn across QE releases; verify it against the
    documentation of your installed QE version before a production run.
    Run with neb.x (or "pw.x -input neb.in" on QE builds where NEB was
    merged into pw.x, depending on your version).

    :param first_struct: initial-state structure file
    :param last_struct: final-state structure file (same atom count/order)
    """
    _params = functions.parser_inputpara(params) if params else {}

    _first = structure.struct(first_struct)
    _last = structure.struct(last_struct)
    if _first.total_atom != _last.total_atom:
        raise functions.MkitsError("NEB endpoints must have the same number of atoms.")

    _atomic_indices = [int(_) for _ in _first.position[1:, 0]]
    _atomic_type = sorted(set(_atomic_indices))
    _atomic_symbols = [database.atom_data[_][1] for _ in _atomic_type]

    if wpath == "./":
        wpath = os.path.abspath("./")
    _wkdir = os.path.join(wpath, wname)
    os.makedirs(_wkdir, exist_ok=True)

    _upf_info = sssp_v2_resolve_upf(_atomic_symbols, upfpath, os.path.join(_wkdir, "pseudo"))
    _val_electron = sum(_upf_info[_s]["valence"] for _s in _atomic_symbols)
    _is_uspp = any(_upf_info[_s]["is_uspp"] for _s in _atomic_symbols)

    # --- intermediate images: linear interpolation in fractional coordinates ---
    _n_mid = max(0, num_of_images - 2)
    _images = [_first]
    for _i in range(1, _n_mid + 1):
        _t = _i / float(num_of_images - 1)
        _mid = structure.struct("none")
        _mid.lattice9 = _first.lattice9 * (1 - _t) + _last.lattice9 * _t
        _mid.lattice6 = functions.lattice_conversion(_mid.lattice9)
        _mid.total_atom = _first.total_atom
        _frac = _first.position[1:, 1:4] * (1 - _t) + _last.position[1:, 1:4] * _t
        _cart = functions.frac2cart(_mid.lattice9, _frac)
        _rest = _first.position[1:, 7:11]
        _mid.position = np.vstack((
            np.zeros((1, 11)),
            np.hstack((_first.position[1:, 0:1], _frac, _cart, _rest))
        ))
        _images.append(_mid)
    _images.append(_last)

    # --- &PATH ---
    _path_params = {
        "restart_mode": "from_scratch",
        "string_method": "neb",
        "nstep_path": int(_params.pop("nstep_path", 100)),
        "opt_scheme": _params.pop("opt_scheme", "broyden"),
        "num_of_images": num_of_images,
        "CI_scheme": _params.pop("CI_scheme", "auto"),
        "k_max": float(_params.pop("k_max", 0.3)),
        "k_min": float(_params.pop("k_min", 0.2)),
        "path_thr": float(_params.pop("path_thr", 0.1)),
        "ds": float(_params.pop("ds", 2.0)),
        "first_last_opt": False,
    }
    for _key in list(_params):
        if _key in _PATH_KEYS:
            _path_params[_key] = _coerce_param_value(_params.pop(_key))

    # --- engine CONTROL/SYSTEM/ELECTRONS: a stripped-down scf-like setup ---
    _ecutwfc = float(_params.pop("ecutwfc", 60.0))
    _control = {
        "calculation": "scf", "restart_mode": "from_scratch",
        "pseudo_dir": "./pseudo", "outdir": "./outdir", "prefix": prefix,
    }
    _system = {
        "ibrav": 0, "nat": _first.total_atom, "ntyp": len(_atomic_type),
        "input_dft": _params.pop("gga", "pbe"),
        "ecutwfc": _ecutwfc, "ecutrho": _ecutwfc * (10.0 if _is_uspp else 4.0),
    }
    if metal:
        _system.update({
            "occupations": "smearing", "smearing": "gaussian",
            "degauss": float(_params.pop("degauss", 0.01)),
            "nbnd": int(_val_electron / 2.0 * float(_params.pop("nbndfactor", 1.3))),
        })
    else:
        _system.update({"occupations": "fixed"})
    _electrons = {"conv_thr": float(_params.pop("conv_thr", 1.0e-8))}

    _lines = ["BEGIN\n", "BEGIN_PATH_INPUT\n"]
    _lines += _render_namelist("PATH", _path_params)
    _lines += ["END_PATH_INPUT\n", "BEGIN_ENGINE_INPUT\n"]
    _lines += _render_namelist("CONTROL", _control)
    _lines += _render_namelist("SYSTEM", _system)
    _lines += _render_namelist("ELECTRONS", _electrons)
    _lines += _atomic_species_card(_atomic_symbols, _upf_info) + ["\n"]
    _mesh = _first.kmesh(kspacing=float(_params.pop("kspacing", 0.30)))
    _lines += _kpoints_automatic_card(_mesh) + ["\n"]
    _lines += ["BEGIN_POSITIONS\n"]
    _labels = ["FIRST_IMAGE"] + ["INTERMEDIATE_IMAGE"] * _n_mid + ["LAST_IMAGE"]
    for _label, _img in zip(_labels, _images):
        _lines.append(_label + "\n")
        _lines += _atomic_positions_card(_img, frac=True, dyn=False)
    _lines += ["END_POSITIONS\n\n"]
    _lines += _cell_parameters_card(_first)
    _lines += ["END_ENGINE_INPUT\n", "END\n"]

    with open(os.path.join(_wkdir, "neb.in"), "w", newline="\n") as f:
        f.writelines(_lines)

    _cmd = "# NEB calculation\n"
    _cmd += "%s neb.x -in neb.in > neb.out 2>&1\n" % execode
    functions.write_runsh(os.path.join(_wkdir, "run_neb.sh"), _cmd)


def qe_epwinput(
        prefix="pwscf",
        outdir="../scf/outdir",
        dvscf_dir="../ph/save",
        wpath="./",
        wname="epw",
        execode="mpirun -np $SLURM_NTASKS",
        nbndsub=None,
        proj="",
        wdata=None,
        transport=False,
        temps=None,
        ncarrier=None,
        params=""
):
    """
    Generate an epw.x input. Two modes, meant to run as a pair in separate
    directories (transport reads the Wannier data the non-transport run
    produces, via a matching `outdir`/`prefix`):

    - transport=False (default): the Wannierize + interpolate-EPC step
      (epwwrite=.true.), covering the common workflow of building the
      Wannier-interpolated electron-phonon matrix elements.
    - transport=True: a mobility/Boltzmann-transport-equation step
      (epwread=.true., wannierize=.false., scattering=.true.,
      int_mob=.true.) that computes phonon-limited carrier mobility --
      LOWER CONFIDENCE than the Wannierize step (see _INPUTEPW_KEYS's
      transport-keyword comment: EPW's transport module output format has
      changed across releases, verify against your installed version).
      `mkits.mobility.epw_mobility_extractor` reads its "epw.out".

    epw.x has a very large, highly workflow-specific parameter set; only
    the commonly used keys are whitelisted (see _INPUTEPW_KEYS) -- extend
    that set for parameters not covered here.

    :param dvscf_dir: the "save" directory produced by a preceding ph.x run
    :param nbndsub: number of Wannier functions/bands to disentangle (must
        match between the transport=False and transport=True runs)
    :param proj: Wannier90-style projection string, eg "Fe:d,As:p"
        (transport=False only)
    :param wdata: optional list of raw Wannier90-block lines, written as
        wdata(1), wdata(2), ... (transport=False only, eg ["bands_plot = .true."])
    :param temps: list of temperatures (K) for the transport run, eg
        [100, 200, 300] -- written unquoted as epw.x's "temps" array
    :param ncarrier: signed carrier concentration for the transport run,
        eg "1E13" (electrons) or "-1E13" (holes), cm^-3 for 3D / cm^-2 for
        2D -- written as a plain decimal literal (the same value, just not
        the original scientific-notation text: it goes through the same
        string->number coercion as any other QEInput.set value)
    """
    _params = functions.parser_inputpara(params) if params else {}

    if wpath == "./":
        wpath = os.path.abspath("./")
    _wkdir = os.path.join(wpath, wname)
    os.makedirs(_wkdir, exist_ok=True)

    _qe = QEInput(("inputepw",))
    _qe.set("prefix", prefix)
    _qe.set("outdir", outdir)
    _qe.set("dvscf_dir", dvscf_dir)
    _qe.set("elph", True)
    if nbndsub is not None:
        _qe.set("nbndsub", nbndsub)

    if transport:
        _qe.set("epwread", _coerce_param_value(_params.pop("epwread", True)))
        _qe.set("epwwrite", _coerce_param_value(_params.pop("epwwrite", False)))
        _qe.set("wannierize", _coerce_param_value(_params.pop("wannierize", False)))
        _qe.set("scattering", _coerce_param_value(_params.pop("scattering", True)))
        _qe.set("int_mob", _coerce_param_value(_params.pop("int_mob", True)))
        if temps is not None:
            _qe.set("temps", " ".join(str(_t) for _t in temps))
        if ncarrier is not None:
            _qe.set("ncarrier", ncarrier)
    else:
        _qe.set("epbwrite", _coerce_param_value(_params.pop("epbwrite", True)))
        _qe.set("epbread", _coerce_param_value(_params.pop("epbread", False)))
        _qe.set("wannierize", _coerce_param_value(_params.pop("wannierize", True)))
        _qe.set("num_iter", int(_params.pop("num_iter", 300)))
        if proj:
            _qe.set("proj", proj)
        if wdata:
            for _i, _line in enumerate(wdata):
                _qe.set("wdata(%d)" % (_i + 1), _line)

    _qe.update(_params)
    _qe.write(_wkdir, "epw.in")

    _cmd = "# epw.x electron-phonon calculation\n"
    _cmd += "%s epw.x -in epw.in > epw.out 2>&1\n" % execode
    functions.write_runsh(os.path.join(_wkdir, "run_epw.sh"), _cmd)


def qe_pwcondinput(
        tran_prefix="cond",
        outdir="./outdir",
        ikind=0,
        wpath="./",
        wname="pwcond",
        execode="mpirun -np $SLURM_NTASKS",
        params=""
):
    """
    Generate a minimal pwcond.x (ballistic conductance) &INPUTCOND input.

    LOW CONFIDENCE / INCOMPLETE: pwcond.x needs a multi-part setup
    (separate scf runs for the left lead, right lead, and scattering
    region, with matching in-plane cells) that this function does not
    build -- it only writes the &INPUTCOND namelist itself, using a small
    keyword whitelist (_INPUTCOND_KEYS) I am less confident is complete
    or fully accurate than the rest of mkits.qe. Cross-check every
    key/value against the pwcond.x documentation of your installed QE
    version before relying on this for production runs.

    :param ikind: 0 = complex band structure only, 1 = one-lead
        transmission, 2 = two-lead (scattering) transmission
    """
    _params = functions.parser_inputpara(params) if params else {}

    if wpath == "./":
        wpath = os.path.abspath("./")
    _wkdir = os.path.join(wpath, wname)
    os.makedirs(_wkdir, exist_ok=True)

    _qe = QEInput(("inputcond",))
    _qe.set("outdir", outdir)
    _qe.set("tran_prefix", tran_prefix)
    _qe.set("ikind", ikind)
    _qe.set("energy0", float(_params.pop("energy0", 0.0)))
    _qe.set("denergy", float(_params.pop("denergy", 0.01)))
    _qe.set("nenergy", int(_params.pop("nenergy", 100)))
    _qe.update(_params)
    _qe.write(_wkdir, "pwcond.in")

    _cmd = "# pwcond.x ballistic conductance calculation\n"
    _cmd += "%s pwcond.x -in pwcond.in > pwcond.out 2>&1\n" % execode
    functions.write_runsh(os.path.join(_wkdir, "run_pwcond.sh"), _cmd)


# ================================================================== #
# post-processing: dos.x / bands.x output (companions to qepost_geninput)
# ================================================================== #
_EFERMI_RE = re.compile(r"EFermi\s*=\s*(-?\d+\.\d+)")


def _read_alat_bohr(pwout):
    """Read the 'lattice parameter (alat) = ... a.u.' line from a pw.x stdout log."""
    with open(pwout, "r") as f:
        for _line in f:
            if "lattice parameter (alat)" in _line:
                return float(_line.split("=")[1].split()[0])
    raise functions.MkitsError("Could not find 'lattice parameter (alat)' in %s." % pwout)


def _read_fildos_efermi(fildos):
    """Read the 'EFermi = ... eV' value from a dos.x fildos header line."""
    with open(fildos, "r") as f:
        _header = f.readline()
    _match = _EFERMI_RE.search(_header)
    if _match is None:
        raise functions.MkitsError("No 'EFermi = ... eV' header found in %s." % fildos)
    return float(_match.group(1))


def dos_extractor(fildos="pwscf.dos", shift_fermi=True, write_data=True, write_to="dos.dat"):
    """
    Extract total DOS from a dos.x fildos output file -- the QE analogue
    of mkits.vasp.dos_extractor's total-DOS part. Unlike VASP's DOSCAR,
    dos.x only ever produces the total DOS (spin-resolved when nspin=2);
    per-atom/per-orbital partial DOS is a separate projwfc.x output, not
    covered here.

    :param shift_fermi: subtract the Fermi energy printed in fildos's own
        header ("EFermi = ... eV")
    :param write_to: written as whitespace-separated columns with a
        one-line header comment

    Return
    ------
    (data, headers): data is a (nedos, ncol) array, headers names each column.
    """
    _efermi = _read_fildos_efermi(fildos)
    _data = np.loadtxt(fildos, comments="#")
    _energy = _data[:, 0] - _efermi if shift_fermi else _data[:, 0]
    # columns after energy are [dos_up, (dos_dn), intdos_up, (intdos_dn)] -- keep dos only
    _nspin = 2 if _data.shape[1] - 1 >= 3 else 1
    _spin_names = ["spin1", "spin2"][:_nspin]

    _columns = [_energy]
    _headers = ["Energy"]
    for _s, _sname in enumerate(_spin_names):
        _columns.append(_data[:, 1 + _s])
        _headers.append("%s_Total" % _sname)

    _final = np.column_stack(_columns)

    if write_data:
        with open(write_to, "w") as f:
            f.write("# " + " ".join("%15s" % _h for _h in _headers) + "\n")
            np.savetxt(f, _final, fmt="%16.8e")

    return _final, _headers


def _read_filband(filband):
    """
    Bulk-read bands.x's filband output: an ordered list of k-points (2*pi/
    alat "tpiba" Cartesian units) and their band energies (eV). Tokens are
    consumed positionally according to the header's declared nbnd/nks,
    rather than split per line, since bands.x wraps each k-point's
    band-energy row across a fixed-width line count that varies with nbnd.

    LIMITATION: only a single spin channel is supported. A spin-polarized
    run writes one filband per spin channel (named depending on how the
    bands.x call was parameterized) -- call this once per file.
    """
    with open(filband, "r") as f:
        _lines = f.readlines()

    _match = re.search(r"nbnd\s*=\s*(\d+)\s*,\s*nks\s*=\s*(\d+)", _lines[0])
    if _match is None:
        raise functions.MkitsError("Could not parse 'nbnd=.., nks=..' header in %s." % filband)
    _nbnd, _nks = int(_match.group(1)), int(_match.group(2))

    _tokens = [float(_t) for _line in _lines[1:] for _t in _line.split()]
    if len(_tokens) != _nks * (3 + _nbnd):
        raise functions.MkitsError(
            "%s does not contain the expected nks*(3+nbnd)=%d data values (found %d)."
            % (filband, _nks * (3 + _nbnd), len(_tokens))
        )

    _kpoints_tpiba = np.zeros((_nks, 3))
    _energies = np.zeros((_nks, _nbnd))
    _i = 0
    for _k in range(_nks):
        _kpoints_tpiba[_k] = _tokens[_i:_i + 3]
        _i += 3
        _energies[_k] = _tokens[_i:_i + _nbnd]
        _i += _nbnd

    return _kpoints_tpiba, _energies


def band_extractor(
        filband="filband",
        pwout="pwband.out",
        efermi=None,
        fildos=None,
        shift_fermi=True,
        write_data=True,
        write_to="band.dat"
):
    """
    Extract band eigenvalues (a k-path distance + energy table) from a
    bands.x filband file -- the QE analogue of mkits.vasp.band_extractor.
    Unlike VASP's EIGENVAL (already fractional-coordinate), filband's
    k-points are in 2*pi/alat Cartesian ("tpiba") units, so `pwout` (the
    pw.x stdout of the SAME band-structure run, eg qe_geninput's own
    "pwband.out") is also needed, to read alat and convert them to the
    Cartesian 1/angstrom (2*pi convention) axis mkits.vasp.band_extractor
    also uses.

    filband carries no Fermi energy either. If `efermi` is not given and
    `shift_fermi` is True, it is read from `fildos`'s own header (see
    dos_extractor); pass `efermi` explicitly, or `shift_fermi=False`, when
    no dos.x output is at hand.

    LIMITATION: single spin channel only, see _read_filband.

    :param write_to: written as "Kpath Spin1_Band1 ..." columns, with a
        comment line listing the cumulative k-path distance of each
        detected high-symmetry point (a kink in the k-path direction)

    Return
    ------
    (kpath, energies, high_sym_points)
        kpath: (nk,) cumulative k-path distance (1/angstrom, 2*pi convention)
        energies: (1, nk, nb) band energies (Fermi-shifted if requested)
        high_sym_points: list of kpath distances where the path direction changes
    """
    _kpoints_tpiba, _energies = _read_filband(filband)

    if shift_fermi and efermi is None:
        if fildos is not None and os.path.exists(fildos):
            efermi = _read_fildos_efermi(fildos)
        else:
            raise functions.MkitsError(
                "shift_fermi=True but no efermi was given and no fildos was provided; "
                "pass efermi explicitly, fildos, or set shift_fermi=False."
            )
    if shift_fermi:
        _energies = _energies - efermi

    _alat_bohr = _read_alat_bohr(pwout)
    _tpiba2ang = (2.0 * np.pi / _alat_bohr) / database.uc_bohr2ang
    _k_cart = _kpoints_tpiba * _tpiba2ang

    _nk = len(_k_cart)
    _kpath = np.zeros(_nk)
    _high_sym = [0.0]
    for _i in range(1, _nk):
        _kpath[_i] = _kpath[_i - 1] + np.linalg.norm(_k_cart[_i] - _k_cart[_i - 1])
        if _i < _nk - 1:
            _v1 = _k_cart[_i] - _k_cart[_i - 1]
            _v2 = _k_cart[_i + 1] - _k_cart[_i]
            _n1, _n2 = np.linalg.norm(_v1), np.linalg.norm(_v2)
            if _n1 > 1e-6 and _n2 > 1e-6 and np.dot(_v1, _v2) / (_n1 * _n2) < 0.999:
                _high_sym.append(_kpath[_i])
    _high_sym.append(_kpath[-1])

    _energies = _energies[np.newaxis, :, :]

    if write_data:
        _nb = _energies.shape[2]
        _headers = ["Kpath"] + ["Spin1_Band%d" % (_b + 1) for _b in range(_nb)]
        _cols = [_kpath] + [_energies[0, :, _b] for _b in range(_nb)]
        with open(write_to, "w") as f:
            f.write("# high-symmetry points: " + " ".join("%.6f" % _h for _h in _high_sym) + "\n")
            f.write("# " + " ".join("%15s" % _h for _h in _headers) + "\n")
            np.savetxt(f, np.column_stack(_cols), fmt="%16.8e")

    return _kpath, _energies, _high_sym


# ================================================================== #
# post-processing: total energy / pw.x stdout band eigenvalues
# ================================================================== #
def total_energy_extractor(pwout):
    """
    Read the converged total energy (eV) from pw.x's own converged-energy
    line -- QE marks a converged (as opposed to intermediate) SCF energy
    with a leading "!", printing one such line per ionic step, so the
    last occurrence in the file is taken (mirroring
    mkits.vasp.total_energy_extractor's OUTCAR TOTEN convention). QE
    reports energies in Ry; converted to eV via database.uc_ry2ev. Used
    by mkits.mobility.elastic_modulus's strain-energy fits, matching
    mkits.vasp.total_energy_extractor's role for VASP.
    """
    _energy_ry = None
    with open(pwout, "r") as f:
        for _line in f:
            _stripped = _line.strip()
            if _stripped.startswith("!") and "total energy" in _stripped:
                _energy_ry = float(_stripped.split("=")[1].split()[0])
    if _energy_ry is None:
        raise functions.MkitsError("No converged ('!'-marked) total energy line found in %s." % pwout)
    return _energy_ry * database.uc_ry2ev


# ================================================================== #
# post-processing: pw.x stdout band eigenvalues (the QE analogue of
# mkits.vasp's EIGENVAL-based post-processing, for mkits.mobility's
# effective-mass pipeline)
# ================================================================== #
_KLINE_RE = re.compile(r"k\s*=\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")


def _take_until_blank(lines):
    """
    Collect whitespace-separated float tokens from the first non-blank
    line up to (not including) the next blank line -- skips any leading
    blank lines (eg the one QE always prints right after a "k = ..."
    header) instead of stopping on them.
    """
    _values = []
    _started = False
    for _line in lines:
        if not _line.strip():
            if _started:
                break
            continue
        _started = True
        _values += [float(_v) for _v in _line.split()]
    return _values


def read_pwout_bands(pwout, rec_basis, only_after_last_scf=True):
    """
    Parse per-k-point band eigenvalues/occupations from a pw.x stdout log
    ("k = kx ky kz (... PWs) bands (ev):" blocks) -- the QE analogue of
    mkits.vasp.read_eigenval.

    LOWER CONFIDENCE than mkits.vasp's EIGENVAL reader: this parses pw.x's
    human-readable stdout, not a fixed-format data file, and only
    non-spin-polarized output is handled (nspin=1 in the return, matching
    read_eigenval's shape convention but never 2). The per-k eigenvalue
    printout for a self-consistent calculation only appears at all when
    &CONTROL verbosity='high' is set (mkits.workflow.qe_workflow_effective_mass
    forces this for its coarse scan) -- verify against your installed QE
    version if this doesn't parse cleanly.

    QE prints "k = ..." in Cartesian tpiba units (2*pi/alat, alat read
    from this same file's "lattice parameter (alat)" line, in bohr) --
    this function converts them back to fractional (crystal) coordinates
    using `rec_basis`, so its return matches read_eigenval's
    fractional-coordinate convention.

    :param rec_basis: (3, 3) reciprocal lattice basis WITH the 2*pi factor
        (mkits.structure.struct.get_rec_basis(with_2pi=True)) -- the SAME
        structure (and therefore the same lattice) that was fed into the
        pw.x calculation this file is the output of
    :param only_after_last_scf: if the file contains one or more "End of
        self-consistent calculation" markers (a real scf run with
        verbosity='high' reprints every k-point's bands at every SCF
        iteration), only the k-point blocks after the LAST such marker
        (the converged, final pass) are read. A non-self-consistent run
        (eg calculation='bands') has no such marker and prints each
        k-point's bands exactly once, so this flag has no effect there.

    Return
    ------
    (kpoints_frac, weights, energies, occupations) -- same shapes/
    conventions as mkits.vasp.read_eigenval (weights are always 1.0,
    since pw.x's stdout band printout doesn't restate k-point weights;
    nspin is always 1).
    """
    with open(pwout, "r") as f:
        _lines = f.readlines()

    _alat_bohr = _read_alat_bohr(pwout)

    _scan_lines = _lines
    if only_after_last_scf:
        _scf_idx = [_i for _i, _l in enumerate(_lines) if "End of self-consistent calculation" in _l]
        if _scf_idx:
            _scan_lines = _lines[_scf_idx[-1]:]

    _k_idx = [_i for _i, _l in enumerate(_scan_lines) if _KLINE_RE.search(_l)]
    if not _k_idx:
        raise functions.MkitsError("No 'k = ...' band blocks found in %s." % pwout)

    _tpiba2ang = (2.0 * np.pi / _alat_bohr) / database.uc_bohr2ang
    _rec_inv = np.linalg.inv(np.asarray(rec_basis, dtype=float))

    _kpoints, _energies_list, _occ_list = [], [], []
    for _n, _start in enumerate(_k_idx):
        _end = _k_idx[_n + 1] if _n + 1 < len(_k_idx) else len(_scan_lines)
        _block = _scan_lines[_start + 1:_end]

        _k_tpiba = np.array([float(_v) for _v in _KLINE_RE.search(_scan_lines[_start]).groups()])
        _k_frac = (_k_tpiba * _tpiba2ang) @ _rec_inv

        _occ_i = next((_i for _i, _l in enumerate(_block) if "occupation numbers" in _l), None)
        _bands = _take_until_blank(_block[:_occ_i] if _occ_i is not None else _block)
        _occ = _take_until_blank(_block[_occ_i + 1:]) if _occ_i is not None else [1.0] * len(_bands)

        _kpoints.append(_k_frac)
        _energies_list.append(_bands)
        _occ_list.append(_occ)

    _nk = len(_kpoints)
    _nb = min(len(_e) for _e in _energies_list)
    _kpoints_arr = np.array(_kpoints)
    _weights = np.ones(_nk)
    _energies = np.zeros((1, _nk, _nb))
    _occupations = np.zeros((1, _nk, _nb))
    for _k in range(_nk):
        _energies[0, _k] = _energies_list[_k][:_nb]
        _occupations[0, _k] = _occ_list[_k][:_nb]

    return _kpoints_arr, _weights, _energies, _occupations


def _qe_min_image_dist(k1, k2):
    """Fractional k-point distance under the minimum-image convention (periodic wraparound at 0/1)."""
    _d = np.asarray(k1) - np.asarray(k2)
    _d = _d - np.round(_d)
    return np.linalg.norm(_d)


def find_band_extrema(
        pwout, 
        rec_basis, 
        tol=0.1, 
        min_kdist=0.05
):
    """
    Scan a full (symmetry-irreducible) coarse-mesh pw.x scf output for the
    global CBM/VBM and any near-degenerate secondary valleys within `tol`
    eV -- the QE analogue of mkits.vasp.find_band_extrema, used by
    mkits.workflow.qe_workflow_effective_mass. See read_pwout_bands for
    the confidence caveat and the required &CONTROL verbosity='high'.

    :param rec_basis: (3, 3) reciprocal lattice basis WITH the 2*pi factor
        (mkits.structure.struct.get_rec_basis(with_2pi=True)), passed through
        to read_pwout_bands
    :param tol: energy window (eV) below the true VBM / above the true
        CBM within which a k-point is reported as a candidate valley
    :param min_kdist: minimum fractional-coordinate separation (minimum-
        image convention) between reported candidates

    Return
    ------
    dict: {"vbm": [...], "cbm": [...]} -- same shape as
        mkits.vasp.find_band_extrema's return (each entry also carries a
        "spin" key, always 0, for API parity).
    """
    _kpoints, _weights, _energies, _occupations = read_pwout_bands(pwout, rec_basis, only_after_last_scf=True)
    _nk, _nb = _energies.shape[1], _energies.shape[2]

    _vbm_all, _cbm_all = [], []
    for _k in range(_nk):
        _occ_mask = _occupations[0, _k] > 0.5
        if np.any(_occ_mask):
            _b = int(np.argmax(np.where(_occ_mask, _energies[0, _k], -np.inf)))
            _vbm_all.append({"kpoint": _kpoints[_k], "energy": float(_energies[0, _k, _b]), "band_index": _b, "spin": 0})
        if np.any(~_occ_mask):
            _b = int(np.argmin(np.where(~_occ_mask, _energies[0, _k], np.inf)))
            _cbm_all.append({"kpoint": _kpoints[_k], "energy": float(_energies[0, _k, _b]), "band_index": _b, "spin": 0})

    if not _vbm_all or not _cbm_all:
        raise functions.MkitsError(
            "Could not find both occupied and unoccupied bands in %s "
            "(check that this is a semiconductor/insulator scf, not a metal)." % pwout
        )

    def _select(_candidates, _is_max):
        _candidates = sorted(_candidates, key=lambda _c: -_c["energy"] if _is_max else _c["energy"])
        _extremum_e = _candidates[0]["energy"]
        _window = [_c for _c in _candidates if abs(_c["energy"] - _extremum_e) < tol]
        _accepted = []
        for _c in _window:
            if all(_qe_min_image_dist(_c["kpoint"], _a["kpoint"]) >= min_kdist for _a in _accepted):
                _accepted.append(_c)
        return _accepted

    return {"vbm": _select(_vbm_all, True), "cbm": _select(_cbm_all, False)}
