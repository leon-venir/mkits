# -*- coding: utf-8 -*

import re
import xml.etree.ElementTree as ET

import numpy as np
import spglib as spg

from mkits import database
from mkits import functions


def _spglib_field(dataset, key):
    """
    Read a field from a spglib symmetry dataset regardless of whether the
    installed spglib version returns an object with attribute access
    (>=2.5) or a plain dict (<2.5, still supported via a deprecated
    __getitem__ shim on newer versions).
    """
    if dataset is None:
        return None
    try:
        return getattr(dataset, key)
    except AttributeError:
        return dataset[key]


class struct(object):
    """
    Structure object.

    Supported formats
    ------------------
    poscar, cif, xsd, xyz, qein, qeout, lammps, wien2k
    (legacy VASP 4.x POSCARs without a species line are not supported)

    Attributes
    ----------
    inp
        The name of the input file, or a dictionary describing the structure.
    calculator: str
        Detected input format.
    title: str
        Structure title/comment line.
    total_atom: int
        Total number of atoms.
    lattice9: 3x3 numpy array
        Cartesian lattice vectors [[ax,ay,az], [bx,by,bz], [cx,cy,cz]].
    lattice6: 6x1 numpy array
        Direct lattice parameters [a, b, c, alpha, beta, gamma].
    position: (n+1)x11 numpy array
        [atomic_index, frac_x, _y, _z, cart_x, _y, _z, dyn_x, _y, _z, magmom]
        Row 0 is an internal buffer row; real atoms start at row 1.

    Functions
    ---------
    sort_atoms, add_dyn, set_magmom, scale_lattice9, write_struct,
    get_rec_basis, add_atom, replace_atom, supercell, transit_axis,
    get_spacegroup, get_symmetry, get_primitive_cell
    """

    def __init__(self, inp):
        """
        :param self: Initialize the class
        :param inp: The input of structures: a file path, a structure dict,
                    or the literal string "none" for an empty structure.
        """
        self.inp = inp
        self.calculator = "none"
        self.title = "none"
        self.lattice9 = np.array([[]])
        self.lattice6 = np.array([])
        self.position = np.zeros([1, 11])
        self.total_atom = 0

        if isinstance(inp, dict):
            self.calculator = inp["calculator"]
        elif isinstance(inp, str) and inp != "none" and "\n" in inp:
            raise functions.MkitsError("Parsing a structure directly from a string is not supported yet.")
        elif inp == "none":
            pass
        else:
            try:
                with open(inp, "r") as f:
                    _struct_lines = f.readlines()
            except OSError as e:
                raise functions.MkitsError("Cannot open the structure file %s: %s" % (inp, e))
            self.__parse_calculator(_struct_lines)
            self.__parse_structure(_struct_lines)
            self.sort_atoms()

    # ================================================================== #
    # format detection
    # ================================================================== #
    def __parse_calculator(self, lines):
        """
        Detect the structure format from the raw file content.

        Return
        ------
        A string: poscar, cif, xsd, xyz, qein, qeout, lammps, wien2k.
        """
        if "open-source Quantum ESPRESSO" in lines[3]:
            self.calculator = "qeout"
        elif "&SYSTEM" in "".join(lines) or "&system" in "".join(lines):
            self.calculator = "qein"
        elif "_cell_length_a" in "".join(lines):
            self.calculator = "cif"
        elif len(lines) > 1 and "LATTICE" in lines[1] and "NONEQUIV.ATOMS" in lines[1]:
            self.calculator = "wien2k"
        elif self.__is_lammps_data(lines):
            self.calculator = "lammps"
        elif len(lines[2].split()) == 3:
            self.calculator = "poscar"
        elif len(lines) >= 3 and lines[0].strip().isdigit() and len(lines[2].split()) >= 4:
            self.calculator = "xyz"
        elif "<?xml" in "".join(lines[:20]) or "AtomisticTreeRoot" in "".join(lines[:20]) \
                or "IdentityMapping" in "".join(lines):
            self.calculator = "xsd"
        else:
            raise functions.MkitsError(
                "Unsupported structure format, supported: poscar, cif, xsd, xyz, qein, qeout, lammps, wien2k."
            )

    @staticmethod
    def __is_lammps_data(lines):
        """Heuristically detect a LAMMPS data file by its header keywords."""
        _head = " ".join(lines[:30])
        return bool(re.search(r"\d+\s+atoms", _head)) and \
            bool(re.search(r"\d+\s+atom types", _head)) and \
            "xlo xhi" in _head

    # ================================================================== #
    # structure parsing
    # ================================================================== #
    def __parse_structure(self, lines):
        """Dispatch to the format-specific parser and populate self.position."""
        if self.calculator == "poscar":
            self.__parse_poscar(lines)
        elif self.calculator == "cif":
            self.__parse_cif(lines)
        elif self.calculator == "xsd":
            self.__parse_xsd(lines)
        elif self.calculator == "xyz":
            self.__parse_xyz(lines)
        elif self.calculator == "qein":
            self.__parse_qein(lines)
        elif self.calculator == "qeout":
            self.__parse_qeout(lines)
        elif self.calculator == "lammps":
            self.__parse_lammps(lines)
        elif self.calculator == "wien2k":
            self.__parse_wien2k(lines)

    def __default_magmom(self, atomic_indices):
        """Look up the default initial magnetic moment for a list of atomic numbers."""
        return np.array([[database.atom_data[int(_z)][4]] for _z in atomic_indices])

    # ------------------------------------------------------------------ #
    # POSCAR
    # ------------------------------------------------------------------ #
    def __parse_poscar(self, lines):
        """Parse a VASP POSCAR file (species-line format only)."""
        self.title = lines[0][:-1]
        _ratio = float(lines[1][:-1])
        _lattice_x = [float(v) for v in lines[2].split()]
        _lattice_y = [float(v) for v in lines[3].split()]
        _lattice_z = [float(v) for v in lines[4].split()]
        self.lattice9 = np.array([_lattice_x, _lattice_y, _lattice_z]) * _ratio
        self.lattice6 = functions.lattice_conversion(self.lattice9)

        # read elements from the 6th/7th line; legacy POSCAR 
        # no species line, VASP 4.x style are not supported
        try:
            _atom_type = lines[5].split()
            _atom_num = np.array([int(v) for v in lines[6].split()])
            _atoms = functions.listcross(_atom_type, _atom_num)
            self.total_atom = int(np.sum(_atom_num))
        except (ValueError, IndexError):
            raise functions.MkitsError(
                "POSCAR without a species line (VASP 4.x style) is not supported."
            )

        _atom_index = np.array([
            database.symbol_map[i] for i in _atoms
        ]).reshape(-1, 1)

        if functions.rmspace(lines[7]).upper()[0] == "S" and len(lines[9].split()) >= 6:
            for _i in range(self.total_atom):
                lines[9 + _i] = lines[9 + _i].upper().replace("F", "0").replace("T", "1")
            _position = np.array([[
                float(v) for v in line.split()[:3]
            ] for line in lines[9: 9 + self.total_atom]])
            _dyn = np.array([[
                float(v) for v in line.split()[3:6]
            ] for line in lines[9: 9 + self.total_atom]])
            _coord_tag = lines[8]
        elif functions.rmspace(lines[7]).upper()[0] in ("D", "C"):
            _position = np.array([[
                float(v) for v in line.split()[:3]
            ] for line in lines[8: 8 + self.total_atom]])
            _dyn = np.ones([self.total_atom, 3])
            _coord_tag = lines[7]
        else:
            raise functions.MkitsError("Cannot determine the selective-dynamics flag of the POSCAR file.")

        _magmom = self.__default_magmom(_atom_index.flatten())

        if functions.rmspace(_coord_tag).upper()[0] == "D":
            _cart = functions.frac2cart(self.lattice9, _position)
            _block = np.hstack((
                _atom_index, _position, _cart, _dyn, _magmom
            ))
        elif functions.rmspace(_coord_tag).upper()[0] == "C":
            _frac = functions.cart2frac(self.lattice6, _position)
            _block = np.hstack((
                _atom_index, _frac, _position, _dyn, _magmom
            ))
        else:
            raise functions.MkitsError("Cannot determine whether POSCAR coordinates are fractional or Cartesian.")

        self.position = np.vstack((self.position, _block))

    # ------------------------------------------------------------------ #
    # CIF
    # ------------------------------------------------------------------ #
    @staticmethod
    def __find_cif_loop(lines, required_headers):
        """
        Find a CIF loop_ block containing at least one of the required headers.

        Return
        ------
        (headers, first_data_line_index), or (None, None) if not found.
        """
        for i, line in enumerate(lines):
            if "loop_" in line:
                _headers = []
                j = i + 1
                while j < len(lines):
                    _stripped = lines[j].strip()
                    if _stripped.startswith("_"):
                        _headers.append(_stripped)
                        j += 1
                    elif _stripped.startswith("#") or _stripped == "":
                        j += 1
                    else:
                        break
                if any(h in _headers for h in required_headers):
                    return _headers, j
        return None, None

    @staticmethod
    def __read_cif_loop_rows(lines, start):
        """Yield the raw data rows of a CIF loop_ block starting at `start`."""
        for line in lines[start:]:
            _stripped = line.strip()
            if not _stripped:
                continue
            if _stripped.startswith("#") or _stripped.startswith("loop_") or _stripped.startswith("_"):
                break
            yield _stripped

    @staticmethod
    def __symop_str_to_matrix(op_str):
        """
        Convert a CIF-style symmetry operation string, eg "-x+1/2,y,z", into
        a (3x3 rotation, 3 translation) pair, without using eval().
        """
        _op_str = op_str.lower().replace("'", "").replace('"', "").strip()
        _parts = _op_str.split(",")
        _rot = np.zeros((3, 3))
        _trans = np.zeros(3)
        if len(_parts) != 3:
            return np.eye(3), np.zeros(3)

        _axis_index = {"x": 0, "y": 1, "z": 2}
        for _row, _expr in enumerate(_parts):
            _expr = _expr.replace(" ", "")
            for _term in re.findall(r"[+-]?[^+-]+", _expr):
                _sign = -1.0 if _term.startswith("-") else 1.0
                _term = _term.lstrip("+-")
                if not _term:
                    continue
                if _term[0] in _axis_index:
                    _rot[_row, _axis_index[_term[0]]] = _sign
                elif "/" in _term:
                    _num, _den = _term.split("/")
                    _trans[_row] += _sign * float(_num) / float(_den)
                else:
                    _trans[_row] += _sign * float(_term)
        return _rot, _trans

    @staticmethod
    def __symops_from_spacegroup_number(sg_number):
        """Look up symmetry operations for a space-group number via spglib's database."""
        for _hall_number in range(1, 531):
            _sg_type = spg.get_spacegroup_type(_hall_number)
            if _sg_type is not None and _spglib_field(_sg_type, "number") == sg_number:
                _dataset = spg.get_symmetry_from_database(_hall_number)
                return list(_spglib_field(_dataset, "rotations")), list(_spglib_field(_dataset, "translations"))
        return [np.eye(3)], [np.zeros(3)]

    @staticmethod
    def __expand_by_symmetry(frac_coords, symbols, rotations, translations, tol=1e-4):
        """
        Apply a set of symmetry operations to representative (asymmetric-unit)
        fractional coordinates and remove duplicate images, building the full
        list of atoms in the unit cell. Shared by the CIF and WIEN2k parsers.

        Parameters
        ----------
        frac_coords: (n, 3) array-like
        symbols: list of str, length n
        rotations: list of (3, 3) array-like
        translations: list of (3,) array-like

        Return
        ------
        (expanded_symbols, expanded_frac_coords)
        """
        _expanded_symbols = []
        _expanded_coords = []
        for _sym, _frac in zip(symbols, frac_coords):
            _site_coords = []
            for _rot, _trans in zip(rotations, translations):
                _new_pos = (np.array(_rot, dtype=float) @ np.array(_frac, dtype=float)
                            + np.array(_trans, dtype=float)) % 1.0
                _is_duplicate = False
                for _existing in _site_coords:
                    _diff = np.abs(_new_pos - _existing)
                    _diff = np.minimum(_diff, 1.0 - _diff)
                    if np.sum(_diff ** 2) < tol:
                        _is_duplicate = True
                        break
                if not _is_duplicate:
                    _site_coords.append(_new_pos)
            for _c in _site_coords:
                _expanded_symbols.append(_sym)
                _expanded_coords.append(_c)
        return _expanded_symbols, np.array(_expanded_coords)

    def __parse_cif(self, lines):
        """Parse a CIF file, expanding the asymmetric unit into the full atom list."""
        _params = {}
        _sg_number = None
        for line in lines:
            _line = line.strip()
            if _line.startswith("_cell_"):
                _parts = _line.split()
                if len(_parts) >= 2:
                    _params[_parts[0]] = float(re.sub(r"\(.+\)", "", _parts[1]))
            elif _line.lower().startswith("_symmetry_int_tables_number") or \
                    _line.lower().startswith("_space_group_it_number"):
                _parts = _line.split()
                if len(_parts) >= 2:
                    _sg_number = int(re.sub(r"\D", "", _parts[1]))

        _a = _params.get("_cell_length_a")
        _b = _params.get("_cell_length_b")
        _c = _params.get("_cell_length_c")
        _alpha = _params.get("_cell_angle_alpha", 90.0)
        _beta = _params.get("_cell_angle_beta", 90.0)
        _gamma = _params.get("_cell_angle_gamma", 90.0)
        if None in (_a, _b, _c):
            raise functions.MkitsError("CIF file is missing one of the cell length tags.")

        self.lattice6 = np.array([_a, _b, _c, _alpha, _beta, _gamma])
        self.lattice9 = functions.lattice_conversion(self.lattice6)

        # symmetry operations: explicit loop, else space-group number via
        # spglib, else P1
        _sym_headers, _sym_start = self.__find_cif_loop(
            lines, ["_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz"]
        )
        if _sym_start is not None:
            _idx_op = _sym_headers.index(
                "_space_group_symop_operation_xyz" if "_space_group_symop_operation_xyz" in _sym_headers
                else "_symmetry_equiv_pos_as_xyz"
            )
            _rotations, _translations = [], []
            for _row in self.__read_cif_loop_rows(lines, _sym_start):
                if "'" in _row:
                    _op_str = _row.split("'")[1]
                else:
                    _parts = _row.split()
                    if len(_parts) <= _idx_op:
                        continue
                    _op_str = _parts[_idx_op]
                _rot, _trans = self.__symop_str_to_matrix(_op_str)
                _rotations.append(_rot)
                _translations.append(_trans)
        elif _sg_number and _sg_number > 1:
            _rotations, _translations = self.__symops_from_spacegroup_number(_sg_number)
        else:
            _rotations, _translations = [np.eye(3)], [np.zeros(3)]

        # asymmetric-unit atomic positions
        _loop_headers, _loop_start = self.__find_cif_loop(lines, ["_atom_site_fract_x"])
        if _loop_start is None:
            raise functions.MkitsError("CIF file does not contain an atomic coordinates loop.")

        _idx_x = _loop_headers.index("_atom_site_fract_x")
        _idx_y = _loop_headers.index("_atom_site_fract_y")
        _idx_z = _loop_headers.index("_atom_site_fract_z")
        if "_atom_site_type_symbol" in _loop_headers:
            _idx_sym = _loop_headers.index("_atom_site_type_symbol")
        elif "_atom_site_label" in _loop_headers:
            _idx_sym = _loop_headers.index("_atom_site_label")
        else:
            raise functions.MkitsError("CIF file does not contain an atom-symbol column.")

        _rep_symbols, _rep_coords = [], []
        for _row in self.__read_cif_loop_rows(lines, _loop_start):
            _parts = _row.split()
            if len(_parts) < len(_loop_headers):
                continue
            _sym = re.sub(r"[^a-zA-Z]", "", _parts[_idx_sym])
            try:
                _fx = float(re.sub(r"\(.+\)", "", _parts[_idx_x]))
                _fy = float(re.sub(r"\(.+\)", "", _parts[_idx_y]))
                _fz = float(re.sub(r"\(.+\)", "", _parts[_idx_z]))
            except ValueError:
                continue
            _rep_symbols.append(_sym)
            _rep_coords.append([_fx, _fy, _fz])

        _symbols, _frac_coords = self.__expand_by_symmetry(
            _rep_coords, 
            _rep_symbols, 
            _rotations, 
            _translations
        )

        self.total_atom = len(_symbols)
        _atom_index = np.array([
            database.symbol_map[s] for s in _symbols
        ]).reshape(-1, 1)
        _cart_coords = functions.frac2cart(self.lattice9, _frac_coords)
        _dyn = np.ones((self.total_atom, 3))
        _magmom = self.__default_magmom(_atom_index.flatten())

        _block = np.hstack((_atom_index, _frac_coords, _cart_coords, _dyn, _magmom))
        self.position = np.vstack((self.position, _block))

    # ------------------------------------------------------------------ #
    # XSD (Materials Studio)
    # ------------------------------------------------------------------ #
    def __parse_xsd(self, lines):
        """Parse a Materials Studio .xsd file (crystal or isolated molecule)."""
        try:
            _root = ET.fromstring("".join(lines))

            _space_group = None
            for _elem in _root.iter():
                if "SpaceGroup" in _elem.tag:
                    _space_group = _elem
                    break

            if _space_group is not None:
                try:
                    _vec_a = [
                        float(v) for v in _space_group.attrib["AVector"].split(",")
                    ]
                    _vec_b = [
                        float(v) for v in _space_group.attrib["BVector"].split(",")
                    ]
                    _vec_c = [
                        float(v) for v in _space_group.attrib["CVector"].split(",")
                    ]
                    self.lattice9 = np.array([_vec_a, _vec_b, _vec_c])
                    self.lattice6 = functions.lattice_conversion(self.lattice9)
                except KeyError:
                    raise functions.MkitsError("XSD file has a SpaceGroup tag but is missing its lattice vectors.")
            else:
                # isolated molecule: no periodicity, use a large dummy box
                self.lattice9 = np.eye(3) * 100.0
                self.lattice6 = np.array([100.0, 100.0, 100.0, 90.0, 90.0, 90.0])

            _atom_indices = []
            _cart_coords = []
            _found_atoms = False
            for _elem in _root.iter():
                if "Atom3d" in _elem.tag:
                    _found_atoms = True
                    _attrib = _elem.attrib

                    _sym = _attrib.get("Components") or _attrib.get("UserSymbol")
                    if not _sym:
                        _sym = re.sub(r"[^a-zA-Z]", "", _attrib.get("Name", ""))
                    if _sym and "," in _sym:
                        _sym = _sym.split(",")[0]

                    _atom_indices.append(database.symbol_map.get(_sym, 0))
                    if "XYZ" in _attrib:
                        _cart_coords.append([
                            float(v) for v in _attrib["XYZ"].split(",")
                        ])
                    else:
                        _cart_coords.append([0.0, 0.0, 0.0])

            if not _found_atoms:
                raise functions.MkitsError("XSD file does not contain any Atom3d tags.")

            self.total_atom = len(_atom_indices)
            _atom_indices = np.array(_atom_indices).reshape(self.total_atom, 1)
            _cart_coords = np.array(_cart_coords)

            if _space_group is not None:
                _frac_coords = functions.cart2frac(self.lattice9, _cart_coords)
            else:
                _frac_coords = np.zeros((self.total_atom, 3))

            _dyn = np.ones((self.total_atom, 3))
            _magmom = self.__default_magmom(_atom_indices.flatten())

            _block = np.hstack((_atom_indices, _frac_coords, _cart_coords, _dyn, _magmom))
            self.position = np.vstack((self.position, _block))

        except functions.MkitsError:
            raise
        except Exception as e:
            raise functions.MkitsError("Failed to parse XSD file: %s" % e)

    # ------------------------------------------------------------------ #
    # XYZ
    # ------------------------------------------------------------------ #
    def __parse_xyz(self, lines):
        """Parse a plain XYZ file (no periodicity, Cartesian coordinates only)."""
        try:
            self.total_atom = int(lines[0].strip())
            self.title = lines[1].strip()

            # XYZ carries no lattice information; use a large dummy box
            self.lattice9 = np.eye(3) * 100.0
            self.lattice6 = np.array([100.0, 100.0, 100.0, 90.0, 90.0, 90.0])

            _atom_indices = []
            _cart_coords = []
            for line in lines[2: 2 + self.total_atom]:
                _parts = line.split()
                if len(_parts) < 4:
                    continue
                _sym_raw = _parts[0]
                if _sym_raw.isdigit():
                    _idx = int(_sym_raw)
                else:
                    _sym = re.sub(r"[^a-zA-Z]", "", _sym_raw)
                    _idx = database.symbol_map.get(_sym, 0)
                _atom_indices.append(_idx)
                _cart_coords.append([float(_parts[1]), float(_parts[2]), float(_parts[3])])

            _n = len(_atom_indices)
            _atom_indices = np.array(_atom_indices).reshape(_n, 1)
            _cart_coords = np.array(_cart_coords)
            _frac_coords = np.zeros((_n, 3))
            _dyn = np.ones((_n, 3))
            _magmom = self.__default_magmom(_atom_indices.flatten())

            _block = np.hstack((
                _atom_indices, 
                _frac_coords, 
                _cart_coords, 
                _dyn, 
                _magmom
            ))
            self.position = np.vstack((self.position, _block))

        except functions.MkitsError:
            raise
        except Exception as e:
            raise functions.MkitsError("Failed to parse XYZ file: %s" % e)

    # ------------------------------------------------------------------ #
    # Quantum ESPRESSO input/output
    # ------------------------------------------------------------------ #
    def __parse_qein_block(self, lines):
        """Split a QE input file into its &BLOCK ... / namelists and cards."""
        _qeinblock = {}
        _blockname = []
        _start_idx = []
        _end_idx = []
        for i, line in enumerate(lines):
            _stripped = functions.rmspace(line)
            if _stripped and _stripped[0] == "&":
                _blockname.append(_stripped[1:].rstrip("\n"))
                _start_idx.append(i)
        for i in _start_idx:
            for j in range(i, len(lines)):
                if functions.rmspace(lines[j]).startswith("/"):
                    _end_idx.append(j)
                    break
        for _name, _s, _e in zip(_blockname, _start_idx, _end_idx):
            _qeinblock[_name] = lines[_s + 1: _e]

        if "SYSTEM" not in _qeinblock:
            raise functions.MkitsError("QE input file does not contain an &SYSTEM namelist.")

        _system_block = functions.parse_inputfile(
            lines=_qeinblock["SYSTEM"], comment_sym="!", assign_sym="=", seprate_sym=","
        )
        try:
            _tot_atom = int(_system_block["nat"])
        except KeyError:
            raise functions.MkitsError("The QE input file does not define 'nat' in &SYSTEM.")

        for i, line in enumerate(lines):
            _upper = line.upper()
            if "ATOMIC_POSITIONS" in _upper:
                _qeinblock["ATOMIC_POSITIONS"] = lines[i: i + _tot_atom + 1]
            elif "K_POINTS" in _upper:
                _qeinblock["K_POINTS"] = lines[i: i + 2]
            elif "CELL_PARAMETERS" in _upper:
                _qeinblock["CELL_PARAMETERS"] = lines[i: i + 4]

        return _qeinblock

    def __parse_qein(self, lines):
        """Parse a Quantum ESPRESSO pw.x input file."""
        _block = self.__parse_qein_block(lines)

        _lattice_lines = _block["CELL_PARAMETERS"]
        self.lattice9 = np.array([[float(v) for v in line.split()] for line in _lattice_lines[1:]])
        if "bohr" in _lattice_lines[0]:
            self.lattice9 *= database.uc_bohr2ang
        self.lattice6 = functions.lattice_conversion(self.lattice9)

        _position_lines = _block["ATOMIC_POSITIONS"]
        self.total_atom = len(_position_lines) - 1

        _symbols = [line.split()[0] for line in _position_lines[1:]]
        _atom_index = np.array([database.symbol_map[s] for s in _symbols]).reshape(-1, 1)

        _ncol = len(_position_lines[1].split())
        if _ncol == 4:
            _coords = np.array([
                [
                    float(v) for v in line.split()[1:4]
                ] for line in _position_lines[1:]
            ])
            _dyn = np.ones((self.total_atom, 3))
        elif _ncol == 7:
            _coords = np.array([
                [
                    float(v) for v in line.split()[1:4]
                ] for line in _position_lines[1:]
            ])
            _dyn = np.array([
                [
                    float(v) for v in line.split()[4:7]
                ] for line in _position_lines[1:]
            ])
        else:
            raise functions.MkitsError("QE ATOMIC_POSITIONS block must have 4 or 7 columns per line.")

        _tag = _position_lines[0]
        if "crystal" in _tag:
            _frac = _coords
            _cart = functions.frac2cart(self.lattice9, _frac)
        elif "bohr" in _tag:
            _cart = _coords * database.uc_bohr2ang
            _frac = functions.cart2frac(self.lattice6, _cart)
        else:
            # angstrom is QE's default unit when no keyword is given
            _cart = _coords
            _frac = functions.cart2frac(self.lattice6, _cart)

        _magmom = self.__default_magmom(_atom_index.flatten())
        _block2 = np.hstack((_atom_index, _frac, _cart, _dyn, _magmom))
        self.position = np.vstack((self.position, _block2))

    def __parse_qeout(self, lines):
        """Parse the final geometry (lattice + atomic positions) from a QE output file."""
        _latticeunit = 1.0
        _nat = 0
        _crystal_axes_idx = 0
        _cartesian_axes_idx = 0
        _atomic_position_idx = []
        _cell_parameter_idx = []

        for i, line in enumerate(lines):
            if "lattice parameter" in line:
                _latticeunit = float(line.split()[-2])
            elif "number of atoms/cell" in line:
                _nat = int(line.split()[-1])
            elif "crystal axes" in line:
                _crystal_axes_idx = i
            elif "Cartesian axes" in line:
                _cartesian_axes_idx = i
            elif "ATOMIC_POSITIONS" in line:
                _atomic_position_idx.append(i)
            elif "CELL_PARAMETERS" in line:
                _cell_parameter_idx.append(i)

        self.total_atom = _nat

        # lattice: prefer the last relaxed CELL_PARAMETERS block, otherwise
        # fall back to the initial crystal axes block (units of alat, bohr)
        if _cell_parameter_idx:
            i = _cell_parameter_idx[-1]
            self.lattice9 = np.array(
                [[
                    float(v) for v in lines[i + 1 + row].split()[:3]
                ] for row in range(3)]
            ) * _latticeunit * database.uc_bohr2ang
        else:
            self.lattice9 = np.array(
                [[
                    float(v) for v in lines[_crystal_axes_idx + 1 + row].split()[3:6]
                ] for row in range(3)]
            ) * _latticeunit * database.uc_bohr2ang
        self.lattice6 = functions.lattice_conversion(self.lattice9)

        # atomic positions: prefer the last relaxed ATOMIC_POSITIONS block,
        # otherwise the initial Cartesian axes block (units of alat, bohr)
        _symbols = []
        _coords = []
        if _atomic_position_idx:
            i = _atomic_position_idx[-1]
            _is_crystal = "crystal" in lines[i]
            _is_bohr = "bohr" in lines[i]
            for line in lines[i + 1: i + 1 + _nat]:
                _parts = line.split()
                _symbols.append(_parts[0])
                _coords.append([float(v) for v in _parts[1:4]])
            _coords = np.array(_coords)
            if _is_crystal:
                _frac = _coords
                _cart = functions.frac2cart(self.lattice9, _frac)
            else:
                _cart = _coords * database.uc_bohr2ang if _is_bohr else _coords
                _frac = functions.cart2frac(self.lattice6, _cart)
        else:
            for line in lines[_cartesian_axes_idx + 3: _cartesian_axes_idx + 3 + _nat]:
                _parts = line.split()
                _symbols.append(_parts[1])
                _coords.append([float(v) for v in _parts[6:9]])
            _cart = np.array(_coords) * _latticeunit * database.uc_bohr2ang
            _frac = functions.cart2frac(self.lattice6, _cart)

        _atom_index = np.array([database.symbol_map[s] for s in _symbols]).reshape(-1, 1)
        _dyn = np.ones((_nat, 3))
        _magmom = self.__default_magmom(_atom_index.flatten())

        _block = np.hstack((_atom_index, _frac, _cart, _dyn, _magmom))
        self.position = np.vstack((self.position, _block))

    # ------------------------------------------------------------------ #
    # LAMMPS
    # ------------------------------------------------------------------ #
    def __parse_lammps(self, lines):
        """
        Parse a LAMMPS data file (atomic/charge atom style, orthogonal or
        triclinic box). Element symbols are recovered from the Masses
        section by matching each type's mass to the closest entry in
        mkits.database.atom_data.
        """
        self.title = lines[0].strip()

        _xlo = _xhi = _ylo = _yhi = _zlo = _zhi = 0.0
        _xy = _xz = _yz = 0.0
        _masses = {}
        _atoms_start = None
        _n_atoms = 0

        for i, line in enumerate(lines):
            _stripped = line.split("#")[0].strip()
            if re.match(r"^\d+\s+atoms$", _stripped):
                _n_atoms = int(_stripped.split()[0])
            elif "xlo" in _stripped and "xhi" in _stripped:
                _xlo, _xhi = (float(v) for v in _stripped.split()[:2])
            elif "ylo" in _stripped and "yhi" in _stripped:
                _ylo, _yhi = (float(v) for v in _stripped.split()[:2])
            elif "zlo" in _stripped and "zhi" in _stripped:
                _zlo, _zhi = (float(v) for v in _stripped.split()[:2])
            elif "xy" in _stripped and "xz" in _stripped and "yz" in _stripped:
                _xy, _xz, _yz = (float(v) for v in _stripped.split()[:3])
            elif _stripped == "Masses":
                j = i + 2
                while j < len(lines) and lines[j].strip():
                    _parts = lines[j].split("#")[0].split()
                    _masses[int(_parts[0])] = float(_parts[1])
                    j += 1
            elif _stripped == "Atoms":
                _atoms_start = i + 2

        if _atoms_start is None:
            raise functions.MkitsError("LAMMPS data file does not contain an Atoms section.")

        self.lattice9 = np.array([
            [_xhi - _xlo, 0.0, 0.0],
            [_xy, _yhi - _ylo, 0.0],
            [_xz, _yz, _zhi - _zlo]
        ])
        self.lattice6 = functions.lattice_conversion(self.lattice9)

        _type_symbol = {}
        for _atype, _mass in _masses.items():
            _best_z = min(range(1, 119), key=lambda z: abs((database.atom_data[z][3] or 1e9) - _mass))
            _type_symbol[_atype] = _best_z

        _atom_indices = []
        _cart_coords = []
        for line in lines[_atoms_start: _atoms_start + _n_atoms]:
            _parts = line.split("#")[0].split()
            if len(_parts) < 5:
                continue
            _atype = int(_parts[1])
            # supports "atomic" (id type x y z) and "charge" (id type q x y z) styles
            _xyz = _parts[-3:]
            _atom_indices.append(_type_symbol.get(_atype, 0))
            _cart_coords.append([float(v) for v in _xyz])

        self.total_atom = len(_atom_indices)
        _atom_indices = np.array(_atom_indices).reshape(-1, 1)
        _cart_coords = np.array(_cart_coords)
        _frac_coords = functions.cart2frac(self.lattice9, _cart_coords)
        _dyn = np.ones((self.total_atom, 3))
        _magmom = self.__default_magmom(_atom_indices.flatten())

        _block = np.hstack((_atom_indices, _frac_coords, _cart_coords, _dyn, _magmom))
        self.position = np.vstack((self.position, _block))

    # ------------------------------------------------------------------ #
    # WIEN2k
    # ------------------------------------------------------------------ #
    def __parse_wien2k(self, lines):
        """
        Parse a WIEN2k case.struct file. If an atom's declared multiplicity
        (MULT) exceeds the number of equivalent positions explicitly listed
        in the file, the missing ones are completed using the symmetry
        operations block at the end of the file.
        """
        _is_bohr = "ang" not in lines[2].lower()
        _a, _b, _c, _alpha, _beta, _gamma = (float(v) for v in lines[3].split()[:6])
        if _is_bohr:
            _a *= database.uc_bohr2ang
            _b *= database.uc_bohr2ang
            _c *= database.uc_bohr2ang
        self.lattice6 = np.array([_a, _b, _c, _alpha, _beta, _gamma])
        self.lattice9 = functions.lattice_conversion(self.lattice6)

        # trailing symmetry-operation block, used as a fallback to complete
        # atom groups whose equivalent positions are not all listed explicitly
        _sym_start = None
        _n_symops = 0
        for i, line in enumerate(lines):
            if "NUMBER OF SYMMETRY OPERATIONS" in line:
                _sym_start = i
                _n_symops = int(line.split()[0])
                break

        if _sym_start is not None:
            _rotations, _translations = [], []
            i = _sym_start + 1
            for _ in range(_n_symops):
                _rows, _trans = [], []
                for _r in range(3):
                    _parts = lines[i].split()
                    _rows.append([float(_parts[0]), float(_parts[1]), float(_parts[2])])
                    _trans.append(float(_parts[3]))
                    i += 1
                _rotations.append(np.array(_rows))
                _translations.append(np.array(_trans))
                i += 1  # skip the operation-index line
        else:
            _rotations, _translations = [np.eye(3)], [np.zeros(3)]

        _symbols = []
        _frac_coords = []
        i = 4
        while i < len(lines) and lines[i].startswith("ATOM"):
            _coord = [float(v) for v in re.findall(r"[XYZ]=\s*(-?\d+\.\d+)", lines[i])]
            _mult = int(re.findall(r"MULT=\s*(\d+)", lines[i + 1])[0])
            i += 2

            _rep_frac = [_coord]
            for _ in range(_mult - 1):
                _rep_frac.append([float(v) for v in re.findall(r"[XYZ]=\s*(-?\d+\.\d+)", lines[i])])
                i += 1

            _sym = re.sub(r"[^A-Za-z]", "", lines[i].split()[0])
            i += 1
            i += 3  # skip the 3-line local rotation matrix block

            if len(_rep_frac) < _mult:
                _exp_syms, _exp_coords = self.__expand_by_symmetry([_rep_frac[0]], [_sym], _rotations, _translations)
                _symbols += _exp_syms
                _frac_coords += _exp_coords.tolist()
            else:
                _symbols += [_sym] * _mult
                _frac_coords += _rep_frac

        self.total_atom = len(_symbols)
        _atom_index = np.array([database.symbol_map[s] for s in _symbols]).reshape(-1, 1)
        _frac_coords = np.array(_frac_coords)
        _cart_coords = functions.frac2cart(self.lattice9, _frac_coords)
        _dyn = np.ones((self.total_atom, 3))
        _magmom = self.__default_magmom(_atom_index.flatten())

        _block = np.hstack((_atom_index, _frac_coords, _cart_coords, _dyn, _magmom))
        self.position = np.vstack((self.position, _block))

    # ================================================================== #
    # in-place structure edits
    # ================================================================== #
    def sort_atoms(self):
        """Sort the atomic coordinates by atomic index. Updates self.position."""
        self.position = self.position[np.argsort(self.position[:, 0]), :]

    def add_dyn(self,
                xmin=-1e8, xmax=1e8,
                ymin=-1e8, ymax=1e8,
                zmin=-1e8, zmax=1e8,
                fix="none",
                move="none",
                frac=True):
        """
        Set atoms enclosed by the given box to fixed (dyn=0) and the rest
        to free (dyn=1); `fix`/`move` additionally override by element
        symbol (fix -> 0, move -> 1). Updates the dynrange columns of
        self.position.
        """
        _coord = self.position[1:, 1:4] if frac else self.position[1:, 4:7]

        for _i in range(len(_coord)):
            if xmax > _coord[_i, 0] > xmin and ymax > _coord[_i, 1] > ymin and zmax > _coord[_i, 2] > zmin:
                self.position[1 + _i, 7:10] = 0.0
            else:
                self.position[1 + _i, 7:10] = 1.0

        if fix != "none" and fix in database.symbol_map:
            _z = database.symbol_map[fix]
            for _i in range(1, self.total_atom + 1):
                if self.position[_i, 0] == _z:
                    self.position[_i, 7:10] = 0.0
        elif move != "none" and move in database.symbol_map:
            _z = database.symbol_map[move]
            for _i in range(1, self.total_atom + 1):
                if self.position[_i, 0] == _z:
                    self.position[_i, 7:10] = 1.0

    def set_magmom(self,
                    magmom=0.0,
                    xmin=-1e8, xmax=1e8,
                    ymin=-1e8, ymax=1e8,
                    zmin=-1e8, zmax=1e8,
                    element="none",
                    frac=True):
        """
        Set the initial magnetic moment of the atoms enclosed by the given
        box, or of every atom of a given element. Updates the magmom column
        of self.position.
        """
        if element != "none" and element in database.symbol_map:
            _z = database.symbol_map[element]
            for _i in range(1, self.total_atom + 1):
                if self.position[_i, 0] == _z:
                    self.position[_i, 10] = magmom
            return

        _coord = self.position[1:, 1:4] if frac else self.position[1:, 4:7]
        for _i in range(len(_coord)):
            if xmax > _coord[_i, 0] > xmin and ymax > _coord[_i, 1] > ymin and zmax > _coord[_i, 2] > zmin:
                self.position[1 + _i, 10] = magmom

    def scale_lattice9(
            self,
            scale=np.array([
                [1, 1, 1],
                [1, 1, 1],
                [1, 1, 1]
            ])
    ):
        """Scale the Cartesian lattice matrix element-wise by `scale`."""
        self.lattice9 = self.lattice9 * scale

    def apply_strain(self, exx=0.0, eyy=0.0, ezz=0.0, exy=0.0, exz=0.0, eyz=0.0):
        """
        Deform the lattice by the given strain tensor components
        (dimensionless, eg 0.01 = 1% strain), keeping fractional
        coordinates fixed -- the clamped-ion strain convention used for
        deformation-potential and elastic-constant calculations (a
        subsequent fixed-cell ionic relaxation, if wanted, is layered on
        top by the caller, not by this method). Shear components are the
        plain tensor strain, not the engineering 2*strain convention.

        :param exx, eyy, ezz: normal strain along a, b, c
        :param exy, exz, eyz: shear strain components
        """
        _eps = np.array([
            [exx, exy, exz],
            [exy, eyy, eyz],
            [exz, eyz, ezz],
        ])
        self.lattice9 = self.lattice9 @ (np.eye(3) + _eps)
        self.lattice6 = functions.lattice_conversion(self.lattice9)
        self.position[1:, 4:7] = functions.frac2cart(self.lattice9, self.position[1:, 1:4])

    # ================================================================== #
    # writers
    # ================================================================== #
    def write_struct(self,
                      fpath: str = "./",
                      fname: str = "pwscf.in",
                      calculator: str = "qein",
                      dyn: bool = False,
                      frac: bool = True,
                      write2file: bool = True):
        """
        Write the structure to a file.

        Parameters
        ----------
        calculator: str
            Output format: poscar, qein, cif, xyz, lammps, wien2k.
        dyn: bool
            Whether to include selective-dynamics / constraint information
            (poscar, qein only).
        frac: bool
            Write fractional (True) or Cartesian (False) coordinates
            (poscar, qein only; cif is always fractional, xyz/lammps/wien2k
            always use their own native convention).
        write2file: bool
            If False, return the generated lines instead of writing them.
        """
        if calculator == "none":
            calculator = self.calculator

        _qesystemblock = {}
        if calculator == "poscar":
            _lines = self.__write_poscar_lines(dyn, frac)
        elif calculator == "qein":
            _qesystemblock, _lines = self.__write_qein_lines(dyn, frac)
        elif calculator == "cif":
            _lines = self.__write_cif_lines()
        elif calculator == "xyz":
            _lines = self.__write_xyz_lines()
        elif calculator == "lammps":
            _lines = self.__write_lammps_lines()
        elif calculator == "wien2k":
            _lines = self.__write_wien2k_lines()
        else:
            raise functions.MkitsError("Unsupported output format: %s" % calculator)

        if write2file:
            with open(fpath + "/" + fname, "w", newline="\n") as f:
                f.writelines(_lines)
        else:
            return _qesystemblock, _lines

    def __write_poscar_lines(self, dyn, frac):
        """Build the lines of a VASP POSCAR file."""
        _coord = self.position[1:, 1:4] if frac else self.position[1:, 4:7]
        _coord = functions.convert_array2strlist(_coord, fmt="{:20.10f}")

        _atomic_index = [int(_) for _ in self.position[1:, 0]]
        _atomic_type = np.sort(np.array(list(set(_atomic_index))))
        _atomic_num = [_atomic_index.count(_) for _ in _atomic_type]
        _atomic_symbols = [database.atom_data[_][1] for _ in _atomic_type]

        _lines = ["%s\n" % self.title, "1.0\n"]
        _lines += functions.convert_high2writeablelist(functions.convert_array2strlist(self.lattice9))
        _lines.append("   ".join(_atomic_symbols) + "\n")
        _lines += functions.convert_high2writeablelist(functions.convert_array2strlist(_atomic_num, fmt="{:<6d}"))

        if dyn:
            _lines.append("Selective dynamics\n")
            _dyn_flags = self.position[1:, 7:10]
            _dyn_flags = np.where(_dyn_flags > 0.5, "T", "F")
            _dyn_flags = functions.convert_array2strlist(_dyn_flags, fmt="{:>5s}")
            _lines.append("Direct\n" if frac else "Cartesian\n")
            _lines += functions.convert_high2writeablelist(functions.hstack_append_list(_coord, _dyn_flags))
        else:
            _lines.append("Direct\n" if frac else "Cartesian\n")
            _lines += functions.convert_high2writeablelist(_coord)

        return _lines

    def __write_qein_lines(self, dyn, frac):
        """Build the lines of a Quantum ESPRESSO pw.x input file."""
        _atomic_index = [int(_) for _ in self.position[1:, 0]]
        _atomic_type = np.sort(np.array(list(set(_atomic_index))))
        _atomic_symbols = [database.atom_data[_][1] for _ in _atomic_index]

        _coord = self.position[1:, 1:4] if frac else self.position[1:, 4:7]
        _coord = functions.convert_array2strlist(_coord, fmt="{:20.10f}")

        _qesystemblock = {"ibrav": 0, "nat": self.total_atom, "ntyp": len(_atomic_type)}

        _lines = ["&SYSTEM\n"]
        for _key in ("ibrav", "nat", "ntyp"):
            _lines.append("%s=%s\n" % (_key, _qesystemblock[_key]))
        _lines.append("/\n\n")

        _symbol_col = functions.convert_array2strlist(
            np.array(_atomic_symbols).reshape(self.total_atom, -1), fmt="{:>6s}"
        )

        _lines.append("ATOMIC_POSITIONS crystal\n" if frac else "ATOMIC_POSITIONS angstrom\n")
        if dyn:
            _dyn_flags = self.position[1:, 7:10]
            _dyn_flags = np.where(_dyn_flags > 0.5, "1.0", "0.0")
            _dyn_flags = functions.convert_array2strlist(_dyn_flags, fmt="{:>5s}")
            _rows = functions.hstack_append_list(_symbol_col, _coord)
            _rows = functions.hstack_append_list(_rows, _dyn_flags)
            _lines += functions.convert_high2writeablelist(_rows)
        else:
            _lines += functions.convert_high2writeablelist(functions.hstack_append_list(_symbol_col, _coord))

        _lines.append("CELL_PARAMETERS angstrom\n")
        _lines += functions.convert_high2writeablelist(functions.convert_array2strlist(self.lattice9))

        return _qesystemblock, _lines

    def __write_cif_lines(self):
        """Build the lines of a CIF file (P1, fractional coordinates)."""
        _safe_title = self.title.replace(" ", "_") if self.title != "none" else "generated_structure"
        _lines = [
            "data_%s\n" % _safe_title,
            "_audit_creation_method   'Generated by mkits.structure'\n\n",
            "_symmetry_space_group_name_H-M    'P 1'\n",
            "_symmetry_Int_Tables_number       1\n",
            "_symmetry_cell_setting            triclinic\n\n",
        ]

        if len(self.lattice6) != 6:
            self.lattice6 = functions.lattice_conversion(self.lattice9)
        _a, _b, _c, _alpha, _beta, _gamma = self.lattice6

        _lines.append("_cell_length_a        %12.6f\n" % _a)
        _lines.append("_cell_length_b        %12.6f\n" % _b)
        _lines.append("_cell_length_c        %12.6f\n" % _c)
        _lines.append("_cell_angle_alpha     %12.6f\n" % _alpha)
        _lines.append("_cell_angle_beta      %12.6f\n" % _beta)
        _lines.append("_cell_angle_gamma     %12.6f\n\n" % _gamma)

        _lines += [
            "loop_\n",
            " _atom_site_label\n",
            " _atom_site_type_symbol\n",
            " _atom_site_fract_x\n",
            " _atom_site_fract_y\n",
            " _atom_site_fract_z\n",
        ]

        _atomic_index = [int(_) for _ in self.position[1:, 0]]
        _atomic_symbols = [database.atom_data[_][1] for _ in _atomic_index]
        _frac_coords = self.position[1:, 1:4]

        _element_counter = {}
        for i, _sym in enumerate(_atomic_symbols):
            _element_counter[_sym] = _element_counter.get(_sym, 0) + 1
            _label = "%s%d" % (_sym, _element_counter[_sym])
            _fx, _fy, _fz = _frac_coords[i]
            _lines.append(" %-8s %-4s %12.8f %12.8f %12.8f\n" % (_label, _sym, _fx, _fy, _fz))

        return _lines

    def __write_xyz_lines(self):
        """Build the lines of a plain XYZ file (Cartesian coordinates)."""
        _lines = ["%d\n" % self.total_atom, "%s\n" % str(self.title).strip()]

        _indices = self.position[1:, 0]
        _cart_coords = self.position[1:, 4:7]
        for i in range(self.total_atom):
            _idx = int(_indices[i])
            try:
                _sym = database.atom_data[_idx][1]
            except IndexError:
                _sym = "X"
            _x, _y, _z = _cart_coords[i]
            _lines.append("%-4s %12.6f %12.6f %12.6f\n" % (_sym, _x, _y, _z))

        return _lines

    def __write_lammps_lines(self):
        """Build the lines of a LAMMPS data file (atomic style, triclinic box)."""
        _atomic_index = [int(_) for _ in self.position[1:, 0]]
        _atomic_type = sorted(set(_atomic_index))
        _type_of = {z: i + 1 for i, z in enumerate(_atomic_type)}

        _cart = self.position[1:, 4:7].astype(float)
        _cart = _cart - _cart.min(axis=0)

        _a, _b, _c = self.lattice9
        _xlo, _ylo, _zlo = 0.0, 0.0, 0.0
        _xhi, _xy, _yhi, _xz, _yz, _zhi = _a[0], _b[0], _b[1], _c[0], _c[1], _c[2]

        _lines = [
            "%s\n\n" % self.title,
            "%d atoms\n" % self.total_atom,
            "%d atom types\n\n" % len(_atomic_type),
            "%.10f %.10f xlo xhi\n" % (_xlo, _xhi),
            "%.10f %.10f ylo yhi\n" % (_ylo, _yhi),
            "%.10f %.10f zlo zhi\n" % (_zlo, _zhi),
            "%.10f %.10f %.10f xy xz yz\n\n" % (_xy, _xz, _yz),
            "Masses\n\n",
        ]
        for z in _atomic_type:
            _lines.append("%d %.6f\n" % (_type_of[z], database.atom_data[z][3] or 0.0))

        _lines.append("\nAtoms\n\n")
        for i, z in enumerate(_atomic_index):
            _x, _y, _w = _cart[i]
            _lines.append("%d %d %.10f %.10f %.10f\n" % (i + 1, _type_of[z], _x, _y, _w))

        return _lines

    def __write_wien2k_lines(self):
        """
        Build the lines of a WIEN2k case.struct file. Symmetry-equivalent
        atoms and the symmetry-operations block are obtained from spglib;
        local rotation matrices are written as identity (see
        mkits.database.local_rot_matrix).
        """
        _dataset = self.get_symmetry(symprec=1e-3)
        _equiv = np.array(_dataset["equivalent_atoms"])
        _rotations = _dataset["rotations"]
        _translations = _dataset["translations"]
        _sg_number = _dataset["number"]

        _atomic_index = self.position[1:, 0].astype(int)
        _frac = self.position[1:, 1:4]

        _unique_reps = sorted(set(_equiv.tolist()))
        _lines = [
            "mkits generated structure\n",
            "P   LATTICE,NONEQUIV.ATOMS: %3d %3d\n" % (len(_unique_reps), _sg_number),
            "MODE OF CALC=RELA unit=bohr\n",
        ]

        _a, _b, _c, _alpha, _beta, _gamma = self.lattice6
        _lines.append("%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n" % (
            _a / database.uc_bohr2ang, _b / database.uc_bohr2ang, _c / database.uc_bohr2ang,
            _alpha, _beta, _gamma
        ))

        for _rep in _unique_reps:
            _members = np.where(_equiv == _rep)[0]
            _z = _atomic_index[_rep]
            _sym = database.atom_data[_z][1]
            _x, _y, _w = _frac[_rep]
            _lines.append("ATOM %3d: X=%10.8f Y=%10.8f Z=%10.8f\n" % (-(_rep + 1), _x, _y, _w))
            _lines.append("          MULT= %2d          ISPLIT= 8\n" % len(_members))
            for _m in _members[1:]:
                _x, _y, _w = _frac[_m]
                _lines.append("    %d: X=%10.8f Y=%10.8f Z=%10.8f\n" % (_rep + 1, _x, _y, _w))
            _lines.append("%-10s NPT=  781  R0=0.00010000 RMT=   2.0000   Z: %6.2f\n" % (_sym, float(_z)))
            _lines.append(database.local_rot_matrix)

        _lines.append("%4d      NUMBER OF SYMMETRY OPERATIONS\n" % len(_rotations))
        for _op_idx, (_rot, _trans) in enumerate(zip(_rotations, _translations)):
            for _row in range(3):
                # columns are space-separated (rather than WIEN2k's native
                # fixed-width, no-space I2 columns) so mkits' own whitespace-
                # based reader can parse the file it just wrote
                _lines.append("%3d %3d %3d %11.8f\n" % (_rot[_row, 0], _rot[_row, 1], _rot[_row, 2], _trans[_row]))
            _lines.append("      %d\n" % (_op_idx + 1))

        return _lines

    # ================================================================== #
    # geometry utilities
    # ================================================================== #
    def get_rec_basis(self, with_2pi=False):
        """Return the reciprocal lattice basis (rows are reciprocal vectors)."""
        if with_2pi:
            return np.linalg.inv(self.lattice9).T * 2 * np.pi
        else:
            return np.linalg.inv(self.lattice9).T

    def is_2d(self, vacuum_ratio=2.0, vacuum_min=12.0):
        """
        Heuristically flag the structure as 2D when the c-axis is both
        longer than an absolute threshold and much longer than a and b
        (the conventional slab-model convention: vacuum stacked along c).
        Shared by mkits.structure.seekpath and by the VASP/QE/... k-mesh
        generators, so every caller agrees on the same definition of "2D".

        :param vacuum_ratio: c is treated as a vacuum axis if c > vacuum_ratio * max(a, b)
        :param vacuum_min: minimum absolute length (angstrom) of the vacuum axis
        """
        _a, _b, _c = self.lattice6[:3]
        return _c > vacuum_min and _c > vacuum_ratio * max(_a, _b)

    def kmesh(self, kspacing=0.3, oddeven="none", kfix=None):
        """
        Compute a Gamma-centered Monkhorst-Pack mesh (n1, n2, n3) from the
        true reciprocal lattice, honoring 2D structures (the vacuum axis,
        assumed to be c, is forced to 1 point). Code-agnostic: VASP, QE and
        other front-ends format the returned tuple into their own KPOINTS
        card.

        :param kspacing: target spacing between k-points (1/angstrom, 2*pi convention)
        :param oddeven: "odd", "even" or "none": force the mesh to an odd/even
            number of points along each non-fixed direction
        :param kfix: optional list of 3 int; -1 means "not fixed", otherwise
            overrides the computed mesh count for that direction
        """
        if oddeven == "odd":
            _oddeven = 1
        elif oddeven == "even":
            _oddeven = 0
        else:
            _oddeven = -1

        _recip_lengths = np.linalg.norm(self.get_rec_basis(with_2pi=True), axis=1)
        _n = [
            int(max(1, functions.round_even_odd(_b / kspacing, _oddeven)))
            for _b in _recip_lengths
        ]

        if self.is_2d():
            _n[2] = 1

        if kfix is not None:
            for _i in range(3):
                if int(kfix[_i]) != -1:
                    _n[_i] = int(kfix[_i])

        return tuple(_n)

    def local_kbox(self, k0_frac, delta=0.01, npoints=5):
        """
        Build an explicit grid of fractional k-points centered on
        `k0_frac`, spaced `delta` apart (fractional reciprocal-lattice
        coordinates) along each reciprocal lattice direction -- used to
        sample a small neighborhood around a band extremum for a full
        effective-mass-tensor fit (mkits.mobility.effective_mass_tensor),
        unlike kmesh() (BZ coverage) or seekpath (a high-symmetry path).
        Honors is_2d() the same way kmesh() does: the vacuum axis (c*) is
        forced to 1 point, since a 2D structure's effective mass is a 2x2
        in-plane tensor, not a 3x3 one.

        :param k0_frac: (3,) fractional coordinates of the extremum
        :param delta: spacing between adjacent grid points, fractional
            reciprocal-lattice coordinates -- convert via get_rec_basis()
            for the absolute 1/angstrom scale this corresponds to
        :param npoints: number of points per sampled dimension; use an odd
            number so k0_frac itself is included exactly at the grid center

        Return
        ------
        (n, 3) array of fractional k-points, n = npoints^3 (3D) or
        npoints^2 (2D, is_2d()==True)
        """
        _offsets = (np.arange(npoints) - npoints // 2) * delta
        _dims = [_offsets, _offsets, np.array([0.0]) if self.is_2d() else _offsets]
        _grid = np.array(np.meshgrid(*_dims, indexing="ij")).reshape(3, -1).T
        return np.array(k0_frac, dtype=float) + _grid

    def add_atom(self, atomic_symbo, position, is_frac=True, magmom=None):
        """
        Append a single atom to the structure.

        Parameters
        ----------
        atomic_symbo: str
            Element symbol.
        position: (3,) array-like
        is_frac: bool
        magmom: float or None
            Initial magnetic moment; defaults to the element's database value.
        """
        _z = database.symbol_map[atomic_symbo]
        if is_frac:
            _frac_pos = np.array(position, dtype=float)
            _cart_pos = functions.frac2cart(self.lattice9, _frac_pos)
        else:
            _cart_pos = np.array(position, dtype=float)
            _frac_pos = functions.cart2frac(self.lattice9, _cart_pos)

        if magmom is None:
            magmom = database.atom_data[_z][4]

        self.total_atom += 1
        _new_line = np.array([
            _z,
            _frac_pos[0], _frac_pos[1], _frac_pos[2],
            _cart_pos[0], _cart_pos[1], _cart_pos[2],
            1, 1, 1,
            magmom
        ])
        self.position = np.vstack((self.position, _new_line))
        self.sort_atoms()

    def replace_atom(
            self,
            ranges="xmin=-1,xmax=-1,ymin=-1,ymax=-1,zmin=-1,zmax=-1",
            atomic_symbo="Ti",
            is_frac=True
    ):
        """
        Replace the element of every atom enclosed by the given box.

        :param ranges: box boundaries, "xmin=..,xmax=..,ymin=..,ymax=..,zmin=..,zmax=.."
        :param atomic_symbo: the new element symbol
        :param is_frac: whether `ranges` is expressed in fractional coordinates
        """
        if not is_frac:
            raise functions.MkitsError("replace_atom currently only supports fractional ranges.")

        _replace_range = {
            "xmin": -1e8, "ymin": -1e8, "zmin": -1e8,
            "xmax": 1e8, "ymax": 1e8, "zmax": 1e8,
        }
        _given = functions.parser_inputpara(ranges)
        for _key in _given:
            _replace_range[_key] = float(_given[_key])

        for _i in range(self.total_atom):
            if _replace_range["xmin"] < self.position[_i + 1, 1] < _replace_range["xmax"] and \
                    _replace_range["ymin"] < self.position[_i + 1, 2] < _replace_range["ymax"] and \
                    _replace_range["zmin"] < self.position[_i + 1, 3] < _replace_range["zmax"]:
                self.position[_i + 1, 0] = database.symbol_map[atomic_symbo]

        self.sort_atoms()

    def supercell(self, super_matrix=[2, 2, 1]):
        """
        Build a supercell by repeating the cell along each lattice vector.

        :param super_matrix: repetition counts along a, b, c.
        """
        _na, _nb, _nc = super_matrix
        _base = self.position[1:].astype(float)
        _n_base = _base.shape[0]

        _shifts = np.array([
            [i, j, k]
            for i in range(_na) for j in range(_nb) for k in range(_nc)
        ])
        _n_shift = _shifts.shape[0]

        _expanded = np.repeat(_base, _n_shift, axis=0)
        _tiled_shifts = np.tile(_shifts, (_n_base, 1))
        _expanded[:, 1:4] = (_expanded[:, 1:4] + _tiled_shifts) / np.array([_na, _nb, _nc])

        self.lattice6 = np.hstack((self.lattice6[:3] * np.array([_na, _nb, _nc]), self.lattice6[3:]))
        self.lattice9 = functions.lattice_conversion(self.lattice6)
        _expanded[:, 4:7] = functions.frac2cart(self.lattice9, _expanded[:, 1:4])

        self.position = np.vstack((np.zeros((1, _expanded.shape[1])), _expanded))
        self.total_atom = _expanded.shape[0]
        self.sort_atoms()

    def transit_axis(
            self,
            new_axis=np.array([
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 1]
            ])
    ):
        """
        Change the lattice to `new_axis` and update fractional coordinates
        accordingly. Atoms falling outside the new cell are dropped.

        :param new_axis: 3x3 numpy array in angstrom
        """
        _cart_pos = self.position[1:, 4:7].astype(float)

        _new_lattice9 = np.array(new_axis, dtype=float)
        _new_frac_pos = functions.cart2frac(_new_lattice9, _cart_pos)

        # fractional tolerance corresponding to 0.01 angstrom
        _tol_angstrom = 0.01
        _inv_new_lattice = np.linalg.inv(_new_lattice9)
        _recip_lengths = np.linalg.norm(_inv_new_lattice, axis=0)
        _frac_tols = _tol_angstrom * _recip_lengths

        _inside_mask = np.all(
            (_new_frac_pos >= -_frac_tols) & (_new_frac_pos <= 1.0 + _frac_tols),
            axis=1
        )
        _new_frac_pos[_inside_mask] = np.clip(_new_frac_pos[_inside_mask], 0.0, 1.0)

        _buffer_row = self.position[0:1, :]
        _kept_atoms_rows = self.position[1:][_inside_mask].copy()
        _kept_atoms_rows[:, 1:4] = _new_frac_pos[_inside_mask]

        self.position = np.vstack((_buffer_row, _kept_atoms_rows))
        self.lattice9 = _new_lattice9
        self.lattice6 = functions.lattice_conversion(self.lattice9)
        self.total_atom = len(_kept_atoms_rows)

        self.sort_atoms()

    # ================================================================== #
    # symmetry (spglib)
    # ================================================================== #
    def __to_spglib_cell(self):
        """Build the (lattice, positions, numbers) tuple required by spglib."""
        _frac_pos = self.position[1:, 1:4].astype(float)
        _numbers = self.position[1:, 0].astype(int).tolist()
        return (self.lattice9.tolist(), _frac_pos.tolist(), _numbers)

    def get_spacegroup(self, symprec=1e-3):
        """Return the space-group symbol and number detected by spglib, eg 'Fm-3m (225)'."""
        return spg.get_spacegroup(self.__to_spglib_cell(), symprec=symprec)

    def get_symmetry(self, symprec=1e-3):
        """
        Return the symmetry dataset detected by spglib as a plain dict with
        the keys: number, international, hall, rotations, translations,
        equivalent_atoms, wyckoffs, std_lattice.
        """
        _dataset = spg.get_symmetry_dataset(self.__to_spglib_cell(), symprec=symprec)
        _keys = ["number", "international", "hall", "rotations", "translations",
                 "equivalent_atoms", "wyckoffs", "std_lattice"]
        return {_k: _spglib_field(_dataset, _k) for _k in _keys}

    def get_primitive_cell(self, symprec=1e-3):
        """Return a new struct instance reduced to the primitive cell found by spglib."""
        _result = spg.find_primitive(self.__to_spglib_cell(), symprec=symprec)
        if _result is None or _result[0] is None:
            raise functions.MkitsError("spglib failed to find a primitive cell for this structure.")
        _lattice, _positions, _numbers = _result

        _prim = struct("none")
        _prim.calculator = self.calculator
        _prim.title = self.title
        _prim.lattice9 = np.array(_lattice)
        _prim.lattice6 = functions.lattice_conversion(_prim.lattice9)

        _numbers = np.array(_numbers).reshape(-1, 1)
        _frac = np.array(_positions)
        _cart = functions.frac2cart(_prim.lattice9, _frac)
        _dyn = np.ones((len(_numbers), 3))
        _magmom = self.__default_magmom(_numbers.flatten())

        _prim.total_atom = len(_numbers)
        _prim.position = np.vstack((
            np.zeros((1, 11)),
            np.hstack((_numbers, _frac, _cart, _dyn, _magmom))
        ))
        _prim.sort_atoms()
        return _prim


class volumetric(struct):
    """
    Volumetric data object (eg VASP CHGCAR-like files: a structure followed
    by a scalar field on a regular grid).

    :param inp: the name of the input structure/volumetric file
    """

    def __init__(self, inp):
        super().__init__(inp)
        self.grid = np.array([0, 0, 0])
        self.voldata = np.array([])
        self.__parse_volumetric()

    def __parse_volumetric(self):
        """Parse the grid data following the structure block."""
        if self.calculator == "poscar":
            self.grid = np.loadtxt(
                self.inp,
                skiprows=self.total_atom + 9,
                max_rows=1
            )
            _total_grid = int(self.grid[0] * self.grid[1] * self.grid[2])
            self.voldata = np.loadtxt(
                self.inp,
                skiprows=self.total_atom + 10,
                max_rows=_total_grid // 5
            ).flatten()
            if _total_grid % 5 != 0:
                self.voldata = np.hstack((
                    self.voldata,
                    np.loadtxt(
                        self.inp,
                        skiprows=self.total_atom + 9 + _total_grid // 5,
                        max_rows=1
                    )
                ))
            self.voldata = self.voldata.reshape((
                int(self.grid[2]),
                int(self.grid[1]),
                int(self.grid[0])
            ))
        else:
            pass

    def write_volumetric(
        self,
        fpath="./",
        fname="pwscf.in",
        calculator="qein",
        dyn=False,
        frac=True
    ):
        """Write the structure followed by the volumetric grid data."""
        _, _lines = super().write_struct(
            fpath, fname, calculator, dyn, frac, write2file=False
        )
        _total_grid = int(self.grid[0] * self.grid[1] * self.grid[2])

        if calculator == "poscar":
            _lines.append("\n")
            with open(fpath + "/" + fname, "w", newline="\n") as f:
                f.writelines(_lines)
                np.savetxt(f, np.array([self.grid]), fmt="%5d")
                _voldata = self.voldata.flatten()[:_total_grid // 5 * 5]
                np.savetxt(f, _voldata.reshape((-1, 5)), fmt="%10.11e")
                if _total_grid % 5 != 0:
                    np.savetxt(f, self.voldata.flatten()[_total_grid // 5 * 5:], fmt="%10.11e")

    def set_value(
        self,
        value,
        xmin=0, xmax=1,
        ymin=0, ymax=1,
        zmin=0, zmax=1,
        frac=True
    ):
        """Set the grid values inside the given fractional box to `value`."""
        if frac:
            xmin = int(xmin * self.grid[0])
            xmax = int(xmax * self.grid[0])
            ymin = int(ymin * self.grid[1])
            ymax = int(ymax * self.grid[1])
            zmin = int(zmin * self.grid[2])
            zmax = int(zmax * self.grid[2])
        self.voldata[zmin:zmax, ymin:ymax, xmin:xmax] = value


class seekpath(struct):
    """
    Generate a standard high-symmetry k-path for band-structure calculations.

    NOTE: the path tables (mkits.database.kpath_3d / kpath_2d) are a first,
    self-built approximation - they distinguish crystal systems (3D) / 2D
    Bravais lattice types (2D) only, without accounting for lattice
    centering (P/I/F/C). They are meant to unblock band-structure workflows
    and are expected to be replaced with literature-verified paths later.

    :param inp: input structure file
    :param is2d: force 2D/3D path selection; auto-detected if None
    :param vacuum_ratio: c is treated as a vacuum axis if c > vacuum_ratio * max(a, b)
    :param vacuum_min: minimum absolute length (angstrom) of the vacuum axis
    :param symprec: symmetry-detection tolerance passed to spglib
    """

    def __init__(self, inp, is2d=None, vacuum_ratio=2.0, vacuum_min=12.0, symprec=1e-3):
        super().__init__(inp)
        self.symprec = symprec
        self.is2d = self.is_2d(vacuum_ratio, vacuum_min) if is2d is None else is2d
        self.lattice_type = self.__classify_2d() if self.is2d else self.__classify_3d()
        self.kpath_points = {}
        self.kpath_segments = []
        self.__build_kpath()

    def __classify_2d(self):
        """Classify the in-plane (a, b) lattice into a 2D Bravais lattice type."""
        _a, _b, _, _, _, _gamma = self.lattice6
        if abs(_a - _b) < 1e-2 * _a:
            if abs(_gamma - 90.0) < 1.0:
                return "square"
            if abs(_gamma - 120.0) < 1.0 or abs(_gamma - 60.0) < 1.0:
                return "hexagonal"
        if abs(_gamma - 90.0) < 1.0:
            return "rectangular"
        return "oblique"

    def __classify_3d(self):
        """Classify the structure into a crystal system via its spglib space-group number."""
        _number = self.get_symmetry(symprec=self.symprec)["number"]
        if _number is None or _number <= 2:
            return "triclinic"
        elif _number <= 15:
            return "monoclinic"
        elif _number <= 74:
            return "orthorhombic"
        elif _number <= 142:
            return "tetragonal"
        elif _number <= 167:
            return "trigonal"
        elif _number <= 194:
            return "hexagonal"
        else:
            return "cubic"

    def __build_kpath(self):
        """Populate self.kpath_points / self.kpath_segments from the placeholder tables."""
        _table = database.kpath_2d if self.is2d else database.kpath_3d
        _entry = _table[self.lattice_type]
        self.kpath_points = _entry["points"]
        self.kpath_segments = _entry["path"]

    def write_kpath(self, fpath="./", fname="KPOINTS_band", code="vasp", kpoints_per_segment=20,
                     write2file=True):
        """
        Write the high-symmetry k-path to a band-structure k-points file, or
        return its lines (`write2file=False`) so a caller (eg mkits.qe) can
        embed them directly into a larger input file instead of a standalone
        KPOINTS/K_POINTS file.

        :param code: "vasp" (KPOINTS, line mode) or "qe" (K_POINTS crystal_b)
        :param kpoints_per_segment: number of k-points interpolated per segment
        """
        if code == "vasp":
            _lines = ["k-path generated by mkits.structure.seekpath\n",
                      "%d\n" % kpoints_per_segment, "Line-mode\n", "Reciprocal\n"]
            for _segment in self.kpath_segments:
                for _label1, _label2 in zip(_segment[:-1], _segment[1:]):
                    _p1 = self.kpath_points[_label1]
                    _p2 = self.kpath_points[_label2]
                    _lines.append("%12.8f%12.8f%12.8f  ! %s\n" % (_p1[0], _p1[1], _p1[2], _label1))
                    _lines.append("%12.8f%12.8f%12.8f  ! %s\n\n" % (_p2[0], _p2[1], _p2[2], _label2))
        elif code == "qe":
            _lines = []
            _all_labels = [label for _segment in self.kpath_segments for label in _segment]
            _lines.append("K_POINTS crystal_b\n")
            _lines.append("%d\n" % len(_all_labels))
            for _segment in self.kpath_segments:
                for _label in _segment[:-1]:
                    _p = self.kpath_points[_label]
                    _lines.append("%12.8f%12.8f%12.8f %d  ! %s\n" % (_p[0], _p[1], _p[2], kpoints_per_segment, _label))
                _last = _segment[-1]
                _p = self.kpath_points[_last]
                _lines.append("%12.8f%12.8f%12.8f %d  ! %s\n" % (_p[0], _p[1], _p[2], 1, _last))
        else:
            raise functions.MkitsError("Unsupported band-structure k-points format: %s" % code)
            return

        if write2file:
            with open(fpath + "/" + fname, "w", newline="\n") as f:
                f.writelines(_lines)
        else:
            return _lines
