# -*- coding: utf-8 -*-
"""Structure building: slab/heterostructure stacking, molecule adsorption,
and NGT-style quasi-crystal tiling canvases.

Class
-----
adsorpt_molecue    : adsorb a molecule onto a slab surface at detected sites.
canvas             : build a 2D quasi-crystal tiling canvas and export a slab.

Functions
---------
stack_struct       : stack two slabs (substrate + support) and add vacuum.
ngt_til, sigma_til, approx1_til, approx2_til, approx3_til, bigapp_til,
hexapp_til, honeycomb_til
                   : preset NGT/approximant tiling canvases.
"""

import os

import numpy as np

from mkits import database
from mkits import functions
from mkits import structure


# ================================================================== #
# heterostructure stacking
# ================================================================== #
def stack_struct(substrate, support, distance=1.7, vacuum=15):
    """
    Stack `support` on top of `substrate` along the c-axis and add a
    vacuum gap above -- the mkits.structure equivalent of the original
    ase.build/ase.io-based interface builder (ase.build.add_vacuum +
    ase.Atoms concatenation).

    Assumes both structures already share compatible in-plane (a, b)
    lattice vectors (eg after independent surface/supercell construction)
    and a c-axis aligned along Cartesian z, the standard slab convention:
    the combined cell keeps substrate's a/b vectors and rebuilds c as
    [0, 0, top_of(support) + vacuum]. Neither assumption is checked.

    :param substrate: bottom layer -- a structure file path (any
        mkits.structure format) or an already-loaded mkits.structure.struct
    :param support: layer stacked on top of substrate -- file path or
        mkits.structure.struct
    :param distance: vertical gap (angstrom) between substrate's top and
        support's bottom
    :param vacuum: vacuum thickness (angstrom) added above support's
        (shifted) top

    Return
    ------
    A new mkits.structure.struct: substrate's atoms followed by support's
    (z-shifted) atoms, sharing one cell.
    """
    _substrate = substrate if isinstance(substrate, structure.struct) else structure.struct(substrate)
    _support = support if isinstance(support, structure.struct) else structure.struct(support)

    _substrate_top = _substrate.position[1:, 6].max()
    _support_bottom = _support.position[1:, 6].min()
    _shift = _substrate_top + distance - _support_bottom
    _support_top = _support.position[1:, 6].max() + _shift

    _interface = structure.struct("none")
    _interface.calculator = _substrate.calculator
    _interface.title = "interface"
    _interface.lattice9 = _substrate.lattice9.copy()
    _interface.lattice9[2] = [0.0, 0.0, _support_top + vacuum]
    _interface.lattice6 = functions.lattice_conversion(_interface.lattice9)

    for _i in range(1, _substrate.total_atom + 1):
        _symbol = database.atom_data[int(_substrate.position[_i, 0])][1]
        _interface.add_atom(
            _symbol, _substrate.position[_i, 4:7], is_frac=False, magmom=_substrate.position[_i, 10]
        )
    for _i in range(1, _support.total_atom + 1):
        _symbol = database.atom_data[int(_support.position[_i, 0])][1]
        _shifted_cart = _support.position[_i, 4:7] + np.array([0.0, 0.0, _shift])
        _interface.add_atom(_symbol, _shifted_cart, is_frac=False, magmom=_support.position[_i, 10])

    return _interface


# ================================================================== #
# molecule adsorption
# ================================================================== #
class adsorpt_molecue(structure.struct):
    """
    Adsorb a molecule onto a slab surface at detected adsorption sites,
    writing one output structure file per site.

    Detects surface atoms of `element` within `lateral_radius` of
    `center_coord` (periodic lateral images considered) and within
    `direction_thickness` of it along the surface normal, then places
    `molecule` (from mkits.database's mol_* tables) `distance` above (or
    below, for a "-a"/"-b"/"-c" direction) each detected site.

    Attributes
    ----------
    ads_coords_list: list of (3,) arrays
        Cartesian coordinates of the detected adsorption sites.
    """

    def __init__(
            self,
            inp,
            direction="c",
            center_coord=[0.5, 0.5, 0.5],
            distance=2.0,
            lateral_radius=2.0,
            direction_thickness=0.1,
            element="O",
            molecule="h2o",
            write_to="./",
            write_name="tio2",
            calculator="poscar",
            **kwargs
    ):
        """
        :param inp: path to the input structure file (any mkits.structure format)
        :param direction: surface normal axis and sign, eg "c" or "-c"
        :param center_coord: fractional coordinates of the surface search center
        :param distance: vertical bond length (angstrom) between the surface and the molecule
        :param lateral_radius: lateral search radius around center_coord --
            a fractional length if < 1.0, a Cartesian distance (angstrom) if >= 1.0
        :param direction_thickness: fractional-coordinate tolerance along the surface normal
        :param element: element symbol(s) to search for as adsorption sites ("O" or "O,N")
        :param molecule: molecule name, matching one of mkits.database's
            mol_* tables (eg "h2o" -> database.mol_h2o)
        :param write_to: output directory
        :param write_name: output filename prefix
        :param calculator: output format, passed through to write_struct
        """
        super().__init__(inp)
        self.direction = direction.strip().lower()
        self.center_coord = np.array(center_coord, dtype=float)

        self.distance = float(distance)
        self.lateral_radius = float(lateral_radius)
        self.direction_thickness = float(direction_thickness)

        self.element = element.split(",") if isinstance(element, str) else element

        self.molecule = molecule
        self.write_to = write_to
        self.write_name = write_name
        self.calculator = calculator

        self.ads_coords_list = []

        self.__parse_direction()
        self.__calc_ads_coords()
        self.__adsorpt()

    def __parse_direction(self):
        """Parse self.direction into self.normal_axis/self.lateral_axes/self.sign."""
        _direction = self.direction
        if _direction.startswith("-"):
            self.sign = -1.0
            _axis = _direction[1:]
        elif _direction.startswith("+"):
            self.sign = 1.0
            _axis = _direction[1:]
        else:
            self.sign = 1.0
            _axis = _direction

        _axis_map = {"a": (0, [1, 2]), "b": (1, [0, 2]), "c": (2, [0, 1])}
        if _axis not in _axis_map:
            functions.write2log("Unknown direction '%s', defaulting to 'c'." % self.direction)
            _axis = "c"
            self.sign = 1.0
        self.normal_axis, self.lateral_axes = _axis_map[_axis]

    def __calc_ads_coords(self):
        """
        Find surface atoms of `element` within `lateral_radius` of
        `center_coord` (periodic lateral images included) and within
        `direction_thickness` along the surface normal; populates
        self.ads_coords_list with their Cartesian positions.
        """
        _center_frac = self.center_coord
        _center_cart = functions.frac2cart(self.lattice9, _center_frac)

        _raw_vac_vec = self.lattice9[self.normal_axis]
        _vac_unit = (_raw_vac_vec / np.linalg.norm(_raw_vac_vec)) * self.sign

        if self.lateral_radius < 1.0:
            _len_a = np.linalg.norm(self.lattice9[self.lateral_axes[0]])
            _len_b = np.linalg.norm(self.lattice9[self.lateral_axes[1]])
            _search_radius = self.lateral_radius * (_len_a + _len_b) / 2.0
        else:
            _search_radius = self.lateral_radius
        _search_sq = _search_radius ** 2

        _shifts = []
        for _i in [-1, 0, 1]:
            for _j in [-1, 0, 1]:
                _vec = np.zeros(3)
                _vec[self.lateral_axes[0]] = _i
                _vec[self.lateral_axes[1]] = _j
                _shifts.append(np.dot(_vec, self.lattice9))

        self.ads_coords_list = []
        for _i in range(1, self.total_atom + 1):
            _atom_frac = self.position[_i, 1:4]
            _atom_cart = self.position[_i, 4:7]

            _symbol = database.atom_data[int(self.position[_i, 0])][1]
            if _symbol not in self.element:
                continue

            _delta_frac = _atom_frac[self.normal_axis] - _center_frac[self.normal_axis]
            _delta_frac -= np.round(_delta_frac)
            if abs(_delta_frac) > self.direction_thickness:
                continue

            _best_pos = None
            _min_lat_dist_sq = np.inf
            for _shift in _shifts:
                _shifted_pos = _atom_cart + _shift
                _diff = _shifted_pos - _center_cart
                _h_diff = np.dot(_diff, _vac_unit)
                _lat_vec = _diff - _h_diff * _vac_unit
                _lat_dist_sq = np.dot(_lat_vec, _lat_vec)
                if _lat_dist_sq < _search_sq and _lat_dist_sq < _min_lat_dist_sq:
                    _min_lat_dist_sq = _lat_dist_sq
                    _best_pos = _shifted_pos

            if _best_pos is not None:
                self.ads_coords_list.append(_best_pos)

        functions.write2log("Total adsorption sites generated: %d" % len(self.ads_coords_list))

    def __adsorpt(self):
        """
        For each detected adsorption site, attach `molecule` and write a
        separate output structure file; self.position/self.total_atom are
        restored to the clean slab once every site has been written.
        """
        _mol_name = "mol_" + self.molecule
        if not hasattr(database, _mol_name):
            raise functions.MkitsError("Molecule %s not found in mkits.database." % self.molecule)
        _mol_data = getattr(database, _mol_name)

        _original_position = self.position.copy()
        _original_total_atom = self.total_atom

        if os.path.isdir(self.write_to):
            _base_fpath = self.write_to
        else:
            _base_fpath = os.path.dirname(self.write_to) or "./"

        _mol_xyz = np.array([_row[1:4] for _row in _mol_data[1:]], dtype=float)
        _mol_elements = [database.atom_data[int(_row[0])][1] for _row in _mol_data[1:]]

        _raw_vac_vec = self.lattice9[self.normal_axis]
        _vac_unit = _raw_vac_vec / np.linalg.norm(_raw_vac_vec)

        for _idx, _base_ads_coord in enumerate(self.ads_coords_list):
            self.position = _original_position.copy()
            self.total_atom = _original_total_atom

            if self.direction.startswith("-"):
                # face downwards: mirror the molecule across x, displace down
                _oriented_xyz = _mol_xyz * np.array([1.0, -1.0, -1.0])
                _displacement = -self.distance * _vac_unit
            else:
                _oriented_xyz = _mol_xyz
                _displacement = self.distance * _vac_unit

            for _i, _symbol in enumerate(_mol_elements):
                _abs_pos = _base_ads_coord + _oriented_xyz[_i] + _displacement
                self.add_atom(_symbol, _abs_pos, is_frac=False)

            _site_fname = "%s_%s_ads_%d" % (self.write_name, self.molecule, _idx + 1)
            self.write_struct(fpath=_base_fpath, fname=_site_fname, calculator=self.calculator)

        self.position = _original_position.copy()
        self.total_atom = _original_total_atom


# ================================================================== #
# NGT-style quasi-crystal tiling canvas
# ================================================================== #
class canvas(object):
    """
    A 2D tiling canvas built by attaching polygon patterns (square,
    triangle, rhombus, ...) edge-to-edge, used to construct NGT-style
    quasi-crystal approximant tilings (see ngt_til, sigma_til, ... below),
    which build_crystal then truncates and exports as a periodic slab.

    Attributes
    ----------
    coord: (n, 3) numpy array
        Atom coordinates (2D tiling plane in x/y; z used for buckling).
    atom: (n,) numpy array
        Element symbols, same order as coord.
    edge_coordinates: (m, 6) numpy array
        [x0, y0, z0, x1, y1, z1] per open edge, available for add_pattern.
    truncated_edge: (k, 2) numpy array
        The polygon boundary set by truncate(), used by build_crystal.
    export_value: dict
        Populated by export_val()/build_crystal() for external inspection.
    """

    def __init__(self, center_pattern):
        """
        :param center_pattern: dict with "corner" (int), "coord" (Nx3
            array), "atoms" (N array of element symbols) -- the starting
            polygon pattern, eg a square/triangle/rhombus unit read from
            an xyz file:
                square_pos = np.loadtxt("square.xyz", skiprows=2, usecols=[1, 2, 3])
                square_atom = np.loadtxt("square.xyz", skiprows=2, usecols=[0], dtype="str")
                pattern_square = {"corner": 4, "coord": square_pos, "atoms": square_atom}
        """
        self.coord = center_pattern["coord"]
        self.atom = center_pattern["atoms"]
        self.edge_list = np.array([0])
        self.edge_used = np.array([0])
        self.edge_max = 0
        self.flip = np.array([0])
        self.pattern_index = np.array([0])
        self.pattern_center = np.zeros(3)
        self.edge_coordinates = np.zeros(6)
        _pos_head_tail = np.vstack((self.coord[:center_pattern["corner"]], self.coord[0]))
        for _i in range(int(center_pattern["corner"])):
            self.edge_max += 1
            self.edge_list = np.append(self.edge_list, self.edge_max)
            self.edge_coordinates = np.vstack((
                self.edge_coordinates, np.hstack((_pos_head_tail[_i], _pos_head_tail[_i + 1]))
            ))
            self.flip = np.append(self.flip, 0)
        self.truncated = False
        self.truncated_edge = np.zeros(2)
        self.export_value = {}

    def add_pattern(self, pattern, edge_index):
        """
        Attach `pattern` to the open edge at `edge_index`, mirroring/
        rotating it so its first edge matches that edge's direction, then
        appends its atoms and registers its own open edges.

        :param pattern: dict, same shape as canvas()'s center_pattern
        :param edge_index: index into self.edge_coordinates to attach to
        """
        _edge_num = pattern["corner"]
        _edge_init_position = self.edge_coordinates[edge_index]

        _flip = self.flip[edge_index]
        _flip = 0 if _flip else 1

        if _flip:
            _reverse = True
            _clockwise = _edge_init_position[1] <= _edge_init_position[4]
        else:
            _reverse = False
            _clockwise = not (_edge_init_position[1] <= _edge_init_position[4])

        if _reverse:
            _vr = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
            _coord = (_vr @ pattern["coord"].T).T
        else:
            _coord = pattern["coord"]

        _vector_init = (_coord[1] - _coord[0])[0:2]
        _vector_fina = _edge_init_position[3:5] - _edge_init_position[0:2]
        _angle = functions.vector_angle(_vector_init, _vector_fina, "rad")
        if _clockwise:
            _vx = np.array([[np.cos(_angle), np.sin(_angle), 0],
                            [-np.sin(_angle), np.cos(_angle), 0],
                            [0, 0, 1]])
        else:
            _vx = np.array([[np.cos(_angle), -np.sin(_angle), 0],
                            [np.sin(_angle), np.cos(_angle), 0],
                            [0, 0, 1]])
        _coord_rotated = (_vx @ _coord.T).T

        _coord_translated = _coord_rotated + np.append(_edge_init_position[0:2], 0)
        self.coord = np.vstack((self.coord, _coord_translated))
        self.atom = np.hstack((self.atom, pattern["atoms"]))

        _pos_head_tail = np.vstack((_coord_translated[0:_edge_num], _coord_translated[0]))
        self.edge_used = np.append(self.edge_used, edge_index)

        for _i in range(int(pattern["corner"])):
            self.edge_max += 1
            self.edge_list = np.append(self.edge_list, self.edge_max)
            self.edge_coordinates = np.vstack((
                self.edge_coordinates, np.hstack((_pos_head_tail[_i], _pos_head_tail[_i + 1]))
            ))
            self.flip = np.append(self.flip, _flip)

        self.pattern_index = np.append(self.pattern_index, self.pattern_index[-1] + 1)
        self.pattern_center = np.vstack((
            self.pattern_center, np.average(_coord_translated[:pattern["corner"]], axis=0)
        ))

    def plot_canvas(self,
                    show_pattern_index=0,
                    show_edge_index=True,
                    save_fig="none",
                    fig_size="none"):
        """
        Visualize the tiling: open edges as arrows, optional edge/pattern
        index labels, and the truncation boundary if truncate() was called.

        :param show_pattern_index: 0 -> don't show, 1 -> ascending order,
            -1 -> descending order
        :param fig_size: eg [20, 20]
        :param save_fig: if given, save the figure under this filename
        """
        import matplotlib.pyplot as plt

        plt.figure(figsize=[10, 10])
        if fig_size != "none":
            plt.figure(figsize=fig_size)

        for _i in range(len(self.edge_coordinates) - 1):
            if _i + 1 in self.edge_used:
                continue
            plt.arrow(
                self.edge_coordinates[_i + 1, 0], self.edge_coordinates[_i + 1, 1],
                self.edge_coordinates[_i + 1, 3] - self.edge_coordinates[_i + 1, 0],
                self.edge_coordinates[_i + 1, 4] - self.edge_coordinates[_i + 1, 1],
                width=0.1
            )
            if show_edge_index:
                plt.text(
                    (self.edge_coordinates[_i + 1, 3] + self.edge_coordinates[_i + 1, 0]) / 2,
                    (self.edge_coordinates[_i + 1, 4] + self.edge_coordinates[_i + 1, 1]) / 2,
                    str(self.edge_list[_i + 1])
                )
        plt.axis("equal")

        if self.truncated:
            plt.plot(self.truncated_edge[:, 0], self.truncated_edge[:, 1],
                     "--", color="tab:red", linewidth=2.0)
            for _i in range(len(self.truncated_edge) - 1):
                plt.text(self.truncated_edge[_i, 0], self.truncated_edge[_i, 1],
                         str(_i), color="tab:red", fontsize=18)

        if show_pattern_index in (1, -1):
            for _i in range(1, len(self.pattern_index)):
                plt.text(
                    self.pattern_center[_i, 0], self.pattern_center[_i, 1],
                    str(self.pattern_index[_i * show_pattern_index]),
                    color="tab:red", fontsize=16
                )

        if save_fig != "none":
            plt.savefig(save_fig)

    def savexyz(self, out):
        """Write self.coord/self.atom as a plain xyz file."""
        with open(out, "w", newline="\n") as f:
            f.write("%d\n\n" % len(self.coord))
            for _i in range(len(self.coord)):
                f.write("%s%23.15f%23.15f%23.15f\n" % (
                    "{:<2}".format(self.atom[_i]), self.coord[_i, 0], self.coord[_i, 1], self.coord[_i, 2]
                ))

    def points_in_polygon(self, polygon):
        """Drop every atom (from index 1 onward) outside the given 2D polygon."""
        from shapely.geometry import Point, Polygon

        _polygon = Polygon(polygon)
        _atom_to_del = []
        for _i in range(1, len(self.atom)):
            if not Point(self.coord[_i, 0:2]).within(_polygon):
                _atom_to_del.append(_i)
        self.atom = np.delete(self.atom, list(set(_atom_to_del)))
        self.coord = np.delete(self.coord, list(set(_atom_to_del)), axis=0)

    def truncate(self, method, inpara):
        """
        Truncate the tiling boundary, keeping only atoms inside it.

        :param method: "cube" (inpara = [xmin, xmax, ymin, ymax]), "edge"
            (inpara = (edge_indices, [0|1 per edge: head/tail], offsetx,
            offsety)), or "midp" (inpara = (edge_indices, offsetx,
            offsety), using each edge's midpoint)
        """
        self.truncated = True
        if method == "midp":
            _truncated_edge = np.zeros(2)
            _offsetx = np.zeros(1)
            _offsety = np.zeros(1)
            for _i in range(len(inpara[0])):
                _edge_points = np.average(self.edge_coordinates[inpara[0][_i]].reshape((-1, 3)), axis=0)
                _truncated_edge = np.vstack((_truncated_edge, _edge_points[0:2]))
                _offsetx = np.append(_offsetx, inpara[1][_i])
                _offsety = np.append(_offsety, inpara[2][_i])

            _truncated_edge = _truncated_edge + (np.vstack((_offsetx, _offsety))).T
            self.points_in_polygon(_truncated_edge[1:])
            self.truncated_edge = np.vstack((_truncated_edge[1:], _truncated_edge[1]))

        elif method == "cube":
            _xmin, _xmax, _ymin, _ymax = inpara
            _truncated_edge = np.array([
                [_xmin, _ymin], [_xmax, _ymin], [_xmax, _ymax], [_xmin, _ymax]
            ])
            self.points_in_polygon(_truncated_edge)
            self.truncated_edge = np.vstack((_truncated_edge, _truncated_edge[0]))

        elif method == "edge":
            _truncated_edge = np.zeros(2)
            _offsetx = np.zeros(1)
            _offsety = np.zeros(1)
            for _i in range(len(inpara[0])):
                if inpara[1][_i] == 0:
                    _edge_points = self.edge_coordinates[inpara[0][_i]][0:2]
                elif inpara[1][_i] == 1:
                    _edge_points = self.edge_coordinates[inpara[0][_i]][3:5]
                _truncated_edge = np.vstack((_truncated_edge, _edge_points))
                _offsetx = np.append(_offsetx, inpara[2][_i])
                _offsety = np.append(_offsety, inpara[3][_i])

            _truncated_edge = _truncated_edge + (np.vstack((_offsetx, _offsety))).T
            self.points_in_polygon(_truncated_edge[1:])
            self.truncated_edge = np.vstack((_truncated_edge[1:], _truncated_edge[1]))
        else:
            raise functions.MkitsError("Unknown truncate method: %s (expected cube, edge, or midp)." % method)

    def get_edge_point(self, edge, initend):
        """
        NOTE: relies on functions.trans_reflection_xy, which is not
        implemented anywhere in mkits yet -- this method (and
        reflection_transform, which calls it) is ported from the original
        but is not currently usable.
        """
        _edge_point_coord = np.zeros(2)
        for _i in range(len(edge)):
            if initend[_i] == 0:
                _edge_point_coord = np.vstack((_edge_point_coord, self.edge_coordinates[0:2]))
            elif initend[_i] == 1:
                _edge_point_coord = np.vstack((_edge_point_coord, self.edge_coordinates[3:5]))
            else:
                raise functions.MkitsError("initend entries must be 0 (initial point) or 1 (end point).")
        return _edge_point_coord

    def reflection_transform(self, edge, initend):
        """See get_edge_point's note -- not currently usable."""
        _edge_points_coord = self.get_edge_point(edge, initend)
        _edge_coordinates_reflected = functions.trans_reflection_xy(
            self.edge_coordinates, "2ps", tuple(_edge_points_coord)
        )
        _coord_reflected = functions.trans_reflection_xy(
            self.coord, "2ps", tuple(_edge_points_coord)
        )
        self.atom = np.hstack((self.atom, self.atom))
        self.edge_coordinates = np.vstack((self.edge_coordinates, _edge_coordinates_reflected))
        self.coord = np.vstack((self.coord, _coord_reflected))

    def build_crystal(self, out, vacuum, a_index, b_index, fmt="vasp",
                      sort=True, scale_a=1.0, scale_b=1.0):
        """
        Cut a periodic 2D slab out of the tiling along the a_index/b_index
        edges of the truncation boundary and export it as a
        mkits.structure struct (c-axis is pure vacuum spacing, not a real
        periodic direction).

        :param out: output file path (directory + filename)
        :param vacuum: c-axis cell height (angstrom) -- vacuum spacing,
            since the tiling itself is 2D
        :param a_index, b_index: (start, end) truncated_edge point
            indices defining the a/b lattice vectors -- must share the
            same start point
        :param fmt: output format, currently only "vasp" (POSCAR) is supported
        """
        _vector_a = self.truncated_edge[a_index[1], 0:2] - self.truncated_edge[a_index[0], 0:2]
        _vector_b = self.truncated_edge[b_index[1], 0:2] - self.truncated_edge[b_index[0], 0:2]

        if a_index[0] != b_index[0]:
            raise functions.MkitsError("Lattice a and b must start from the same point.")
        _angle = functions.vector_angle(_vector_a, np.array([1, 0]), "rad")

        # translate so a_index[0] sits at the origin
        self.coord = self.coord - np.append(self.truncated_edge[a_index[0]], 0)

        if _vector_a[1] > 0:
            _vx = np.array([[np.cos(_angle), np.sin(_angle), 0],
                            [-np.sin(_angle), np.cos(_angle), 0],
                            [0, 0, 1]])
        else:
            _vx = np.array([[np.cos(_angle), -np.sin(_angle), 0],
                            [np.sin(_angle), np.cos(_angle), 0],
                            [0, 0, 1]])
        self.coord = (_vx @ self.coord.T).T

        _vector_a = (_vx @ np.append(_vector_a, 0).T).T
        _vector_b = (_vx @ np.append(_vector_b, 0).T).T

        # drop periodic-image duplicates at the a/b cell boundaries
        _too_close = []
        for _i in range(len(self.atom)):
            if abs(self.coord[_i, 0]) < 1:
                for _j in range(len(self.atom)):
                    if np.linalg.norm(self.coord[_i] + (_vector_a - self.coord[_j])) < 1:
                        _too_close.append(_j)
        for _i in range(len(self.atom)):
            if abs(self.coord[_i, 1]) < 1:
                for _j in range(len(self.atom)):
                    if np.linalg.norm(self.coord[_i] + (_vector_b - self.coord[_j])) < 1:
                        _too_close.append(_j)
        self.atom = np.delete(self.atom, list(set(_too_close)))
        self.coord = np.delete(self.coord, list(set(_too_close)), axis=0)

        if sort:
            self.sort_atom()

        _crystal = structure.struct("none")
        _crystal.title = "canvas_tiling"
        _crystal.lattice9 = np.array([_vector_a.tolist(), _vector_b.tolist(), [0.0, 0.0, vacuum]])
        _crystal.lattice6 = functions.lattice_conversion(_crystal.lattice9)
        for _i in range(len(self.atom)):
            _crystal.add_atom(str(self.atom[_i]), self.coord[_i], is_frac=False)

        self.export_value["crystal_lattice_a"] = _vector_a
        self.export_value["crystal_lattice_b"] = _vector_b

        _calculator = {"vasp": "poscar"}.get(fmt)
        if _calculator is None:
            raise functions.MkitsError("Unsupported build_crystal output format: %s (expected 'vasp')." % fmt)
        _fpath, _fname = os.path.split(out)
        _crystal.write_struct(fpath=_fpath or ".", fname=_fname, calculator=_calculator)

    def sort_atom(self):
        """Sort self.atom/self.coord by atomic number."""
        _atom_index = np.array([database.symbol_map[_a] for _a in self.atom])
        _order = np.argsort(_atom_index)
        self.atom = np.array(self.atom)[_order]
        self.coord = self.coord[_order]

    def export_val(self):
        """Return self.export_value (atom/coord/edge, plus the crystal lattice vectors after build_crystal)."""
        self.export_value["atom"] = self.atom
        self.export_value["coord"] = self.coord
        self.export_value["edge"] = self.edge_coordinates
        return self.export_value

    def del_neighbor(self, threshold=1):
        """Merge atom pairs of the same element closer than `threshold` into their midpoint."""
        _too_close_pairs = [[-1, -1]]
        for _i in range(len(self.coord)):
            for _j in range(len(self.coord)):
                if _i != _j:
                    _dis = np.linalg.norm(self.coord[_j] - self.coord[_i])
                    if _dis < threshold and self.atom[_i] == self.atom[_j]:
                        _too_close_pairs += [[_i, _j]]

        _seen = set()
        _too_close_pairs = [
            _x for _x in _too_close_pairs
            if tuple(_x[::-1]) not in _seen and not _seen.add(tuple(_x))
        ]

        _new_atom = np.array(["x"])
        _new_coord = np.zeros(3)
        _old_atom_index = []
        for _i, _j in _too_close_pairs:
            if self.atom[_i] != self.atom[_j]:
                raise functions.MkitsError("Atoms in pair %d/%d have different elements." % (_i, _j))
            _old_atom_index += [_i, _j]
            _new_atom = np.append(_new_atom, self.atom[_i])
            _new_coord = np.vstack((_new_coord, (self.coord[_i] + self.coord[_j]) / 2))

        self.atom = np.delete(self.atom, list(set(_old_atom_index)))
        self.coord = np.delete(self.coord, list(set(_old_atom_index)), axis=0)

        for _i in range(1, len(_new_coord)):
            _put_atom = all(
                np.linalg.norm(
                    _new_coord[_i] - _j
                ) >= threshold for _j in self.coord
            )
            if _put_atom:
                self.atom = np.hstack((self.atom, _new_atom[_i]))
                self.coord = np.vstack((self.coord, _new_coord[_i]))


# ================================================================== #
# preset NGT / approximant tiling canvases
# ================================================================== #
def ngt_til(pattern_square, pattern_triangle, pattern_rhombus):
    """A possible NGT ideal tiling canvas."""
    ngt_tiling = canvas(pattern_square)
    ngt_tiling.add_pattern(pattern_triangle, 1)
    ngt_tiling.add_pattern(pattern_triangle, 2)
    ngt_tiling.add_pattern(pattern_triangle, 3)
    ngt_tiling.add_pattern(pattern_triangle, 4)
    ngt_tiling.add_pattern(pattern_triangle, 9)
    ngt_tiling.add_pattern(pattern_triangle, 16)
    ngt_tiling.add_pattern(pattern_square, 13)
    ngt_tiling.add_pattern(pattern_square, 12)
    ngt_tiling.add_pattern(pattern_triangle, 26)
    ngt_tiling.add_pattern(pattern_triangle, 28)
    ngt_tiling.add_pattern(pattern_square, 10)
    ngt_tiling.add_pattern(pattern_triangle, 6)
    ngt_tiling.add_pattern(pattern_square, 15)
    ngt_tiling.add_pattern(pattern_triangle, 7)
    ngt_tiling.add_pattern(pattern_triangle, 40)
    ngt_tiling.add_pattern(pattern_triangle, 52)
    ngt_tiling.add_pattern(pattern_rhombus, 18)
    ngt_tiling.add_pattern(pattern_triangle, 59)
    ngt_tiling.add_pattern(pattern_square, 62)
    ngt_tiling.add_pattern(pattern_triangle, 58)
    ngt_tiling.add_pattern(pattern_triangle, 69)
    ngt_tiling.add_pattern(pattern_square, 72)
    ngt_tiling.add_pattern(pattern_triangle, 65)
    ngt_tiling.add_pattern(pattern_triangle, 63)
    ngt_tiling.add_pattern(pattern_triangle, 82)
    ngt_tiling.add_pattern(pattern_triangle, 75)
    ngt_tiling.add_pattern(pattern_rhombus, 89)
    ngt_tiling.add_pattern(pattern_square, 80)
    ngt_tiling.add_pattern(pattern_triangle, 96)
    ngt_tiling.add_pattern(pattern_square, 100)
    ngt_tiling.add_pattern(pattern_triangle, 103)
    ngt_tiling.add_pattern(pattern_triangle, 66)
    ngt_tiling.add_pattern(pattern_triangle, 79)
    ngt_tiling.add_pattern(pattern_triangle, 97)
    ngt_tiling.add_pattern(pattern_triangle, 99)
    ngt_tiling.add_pattern(pattern_triangle, 104)
    ngt_tiling.add_pattern(pattern_triangle, 106)
    ngt_tiling.add_pattern(pattern_triangle, 88)
    ngt_tiling.add_pattern(pattern_triangle, 76)
    ngt_tiling.add_pattern(pattern_square, 109)
    ngt_tiling.add_pattern(pattern_square, 113)
    ngt_tiling.add_pattern(pattern_square, 115)
    ngt_tiling.add_pattern(pattern_square, 119)
    ngt_tiling.add_pattern(pattern_square, 121)
    ngt_tiling.add_pattern(pattern_square, 125)
    ngt_tiling.add_pattern(pattern_square, 127)
    ngt_tiling.add_pattern(pattern_square, 131)
    ngt_tiling.add_pattern(pattern_square, 35)
    ngt_tiling.add_pattern(pattern_triangle, 39)
    ngt_tiling.add_pattern(pattern_triangle, 134)
    ngt_tiling.add_pattern(pattern_triangle, 133)
    ngt_tiling.add_pattern(pattern_triangle, 138)
    ngt_tiling.add_pattern(pattern_triangle, 139)
    ngt_tiling.add_pattern(pattern_triangle, 142)
    ngt_tiling.add_pattern(pattern_triangle, 141)
    ngt_tiling.add_pattern(pattern_triangle, 146)
    ngt_tiling.add_pattern(pattern_triangle, 147)
    ngt_tiling.add_pattern(pattern_triangle, 150)
    ngt_tiling.add_pattern(pattern_triangle, 149)
    ngt_tiling.add_pattern(pattern_triangle, 154)
    ngt_tiling.add_pattern(pattern_triangle, 155)
    ngt_tiling.add_pattern(pattern_triangle, 158)
    ngt_tiling.add_pattern(pattern_triangle, 157)
    ngt_tiling.add_pattern(pattern_triangle, 162)
    ngt_tiling.add_pattern(pattern_triangle, 163)
    ngt_tiling.add_pattern(pattern_triangle, 166)
    ngt_tiling.add_pattern(pattern_triangle, 167)
    ngt_tiling.add_pattern(pattern_triangle, 29)
    ngt_tiling.add_pattern(pattern_rhombus, 224)
    ngt_tiling.add_pattern(pattern_triangle, 45)
    ngt_tiling.add_pattern(pattern_triangle, 234)
    ngt_tiling.add_pattern(pattern_rhombus, 236)
    ngt_tiling.add_pattern(pattern_triangle, 46)
    ngt_tiling.add_pattern(pattern_rhombus, 170)
    ngt_tiling.add_pattern(pattern_triangle, 247)
    ngt_tiling.add_pattern(pattern_triangle, 251)
    ngt_tiling.add_pattern(pattern_rhombus, 253)
    ngt_tiling.add_pattern(pattern_triangle, 256)
    ngt_tiling.add_pattern(pattern_triangle, 244)
    ngt_tiling.add_pattern(pattern_triangle, 248)
    ngt_tiling.add_pattern(pattern_triangle, 173)
    ngt_tiling.add_pattern(pattern_square, 237)
    ngt_tiling.add_pattern(pattern_square, 233)
    ngt_tiling.add_pattern(pattern_square, 254)
    ngt_tiling.add_pattern(pattern_square, 250)
    ngt_tiling.add_pattern(pattern_square, 175)
    ngt_tiling.add_pattern(pattern_triangle, 274)
    ngt_tiling.add_pattern(pattern_triangle, 293)
    ngt_tiling.add_pattern(pattern_triangle, 296)
    ngt_tiling.add_pattern(pattern_triangle, 272)
    ngt_tiling.add_pattern(pattern_square, 301)
    ngt_tiling.add_pattern(pattern_square, 295)
    ngt_tiling.add_pattern(pattern_square, 299)
    ngt_tiling.add_pattern(pattern_square, 218)
    ngt_tiling.add_pattern(pattern_rhombus, 209)
    ngt_tiling.add_pattern(pattern_rhombus, 197)
    ngt_tiling.add_pattern(pattern_rhombus, 194)
    ngt_tiling.add_pattern(pattern_rhombus, 185)
    ngt_tiling.add_pattern(pattern_rhombus, 182)
    ngt_tiling.add_pattern(pattern_triangle, 264)
    ngt_tiling.add_pattern(pattern_triangle, 277)
    ngt_tiling.add_pattern(pattern_triangle, 304)
    ngt_tiling.add_pattern(pattern_triangle, 305)
    ngt_tiling.add_pattern(pattern_triangle, 306)
    ngt_tiling.add_pattern(pattern_triangle, 310)
    ngt_tiling.add_pattern(pattern_triangle, 309)
    ngt_tiling.add_pattern(pattern_triangle, 308)
    ngt_tiling.add_pattern(pattern_triangle, 313)
    ngt_tiling.add_pattern(pattern_triangle, 314)
    ngt_tiling.add_pattern(pattern_triangle, 318)
    ngt_tiling.add_pattern(pattern_triangle, 317)
    ngt_tiling.add_pattern(pattern_triangle, 316)
    ngt_tiling.add_pattern(pattern_triangle, 321)
    ngt_tiling.add_pattern(pattern_triangle, 322)
    ngt_tiling.add_pattern(pattern_triangle, 206)
    ngt_tiling.add_pattern(pattern_triangle, 280)
    ngt_tiling.add_pattern(pattern_triangle, 266)
    ngt_tiling.add_pattern(pattern_square, 343)
    ngt_tiling.add_pattern(pattern_square, 353)
    ngt_tiling.add_pattern(pattern_square, 364)
    ngt_tiling.add_pattern(pattern_rhombus, 374)
    ngt_tiling.add_pattern(pattern_triangle, 373)
    ngt_tiling.add_pattern(pattern_triangle, 380)
    ngt_tiling.add_pattern(pattern_triangle, 382)
    ngt_tiling.add_pattern(pattern_triangle, 408)
    ngt_tiling.add_pattern(pattern_triangle, 410)
    ngt_tiling.add_pattern(pattern_triangle, 402)
    ngt_tiling.add_pattern(pattern_triangle, 398)
    ngt_tiling.add_pattern(pattern_triangle, 419)
    ngt_tiling.add_pattern(pattern_square, 385)
    ngt_tiling.add_pattern(pattern_triangle, 434)
    ngt_tiling.add_pattern(pattern_triangle, 435)
    ngt_tiling.add_pattern(pattern_triangle, 325)
    ngt_tiling.add_pattern(pattern_triangle, 326)
    ngt_tiling.add_pattern(pattern_triangle, 330)
    ngt_tiling.add_pattern(pattern_triangle, 329)
    ngt_tiling.add_pattern(pattern_triangle, 333)
    ngt_tiling.add_pattern(pattern_triangle, 334)
    ngt_tiling.add_pattern(pattern_triangle, 338)
    ngt_tiling.add_pattern(pattern_triangle, 337)
    ngt_tiling.add_pattern(pattern_triangle, 288)
    ngt_tiling.add_pattern(pattern_triangle, 289)
    ngt_tiling.add_pattern(pattern_square, 391)
    ngt_tiling.add_pattern(pattern_triangle, 476)
    ngt_tiling.add_pattern(pattern_triangle, 475)
    ngt_tiling.add_pattern(pattern_triangle, 474)
    ngt_tiling.add_pattern(pattern_rhombus, 485)
    ngt_tiling.add_pattern(pattern_square, 482)
    ngt_tiling.add_pattern(pattern_square, 444)
    ngt_tiling.add_pattern(pattern_square, 451)
    ngt_tiling.add_pattern(pattern_square, 460)
    ngt_tiling.add_pattern(pattern_square, 465)
    ngt_tiling.add_pattern(pattern_square, 442)
    ngt_tiling.add_pattern(pattern_triangle, 438)
    ngt_tiling.add_pattern(pattern_rhombus, 439)
    ngt_tiling.add_pattern(pattern_square, 420)
    ngt_tiling.add_pattern(pattern_square, 379)
    ngt_tiling.add_pattern(pattern_square, 416)
    ngt_tiling.add_pattern(pattern_triangle, 523)
    ngt_tiling.add_pattern(pattern_triangle, 527)
    ngt_tiling.add_pattern(pattern_triangle, 530)
    ngt_tiling.add_pattern(pattern_triangle, 515)
    ngt_tiling.add_pattern(pattern_triangle, 540)
    ngt_tiling.add_pattern(pattern_triangle, 534)
    ngt_tiling.add_pattern(pattern_triangle, 538)
    ngt_tiling.add_pattern(pattern_triangle, 531)
    ngt_tiling.add_pattern(pattern_triangle, 512)
    ngt_tiling.add_pattern(pattern_triangle, 513)
    ngt_tiling.add_pattern(pattern_triangle, 496)
    ngt_tiling.add_pattern(pattern_triangle, 447)
    ngt_tiling.add_pattern(pattern_triangle, 457)
    ngt_tiling.add_pattern(pattern_triangle, 462)
    ngt_tiling.add_pattern(pattern_square, 543)
    ngt_tiling.add_pattern(pattern_triangle, 558)
    ngt_tiling.add_pattern(pattern_rhombus, 546)
    ngt_tiling.add_pattern(pattern_square, 556)
    ngt_tiling.add_pattern(pattern_square, 581)
    ngt_tiling.add_pattern(pattern_square, 567)
    ngt_tiling.add_pattern(pattern_square, 571)
    ngt_tiling.add_pattern(pattern_square, 573)
    ngt_tiling.add_pattern(pattern_triangle, 508)
    ngt_tiling.add_pattern(pattern_triangle, 507)
    ngt_tiling.add_pattern(pattern_triangle, 604)
    ngt_tiling.add_pattern(pattern_triangle, 603)
    ngt_tiling.add_pattern(pattern_square, 611)
    ngt_tiling.add_pattern(pattern_triangle, 620)
    ngt_tiling.add_pattern(pattern_triangle, 623)
    ngt_tiling.add_pattern(pattern_triangle, 491)
    ngt_tiling.add_pattern(pattern_rhombus, 626)
    ngt_tiling.add_pattern(pattern_triangle, 492)
    ngt_tiling.add_pattern(pattern_square, 627)

    ngt_tiling.add_pattern(pattern_triangle, 641)
    ngt_tiling.add_pattern(pattern_triangle, 640)
    ngt_tiling.add_pattern(pattern_triangle, 639)
    ngt_tiling.add_pattern(pattern_square, 610)
    ngt_tiling.add_pattern(pattern_triangle, 652)
    ngt_tiling.add_pattern(pattern_triangle, 653)
    ngt_tiling.add_pattern(pattern_triangle, 654)
    ngt_tiling.add_pattern(pattern_square, 663)
    ngt_tiling.add_pattern(pattern_triangle, 665)
    ngt_tiling.add_pattern(pattern_triangle, 666)
    ngt_tiling.add_pattern(pattern_rhombus, 660)
    ngt_tiling.add_pattern(pattern_triangle, 647)
    ngt_tiling.add_pattern(pattern_square, 680)
    ngt_tiling.add_pattern(pattern_triangle, 682)
    ngt_tiling.add_pattern(pattern_triangle, 683)
    ngt_tiling.add_pattern(pattern_triangle, 684)
    ngt_tiling.add_pattern(pattern_square, 614)
    ngt_tiling.add_pattern(pattern_triangle, 599)
    ngt_tiling.add_pattern(pattern_triangle, 600)
    ngt_tiling.add_pattern(pattern_triangle, 601)
    ngt_tiling.add_pattern(pattern_triangle, 616)
    ngt_tiling.add_pattern(pattern_triangle, 699)
    ngt_tiling.add_pattern(pattern_rhombus, 708)
    ngt_tiling.add_pattern(pattern_square, 712)
    ngt_tiling.add_pattern(pattern_square, 365)
    ngt_tiling.add_pattern(pattern_triangle, 592)
    ngt_tiling.add_pattern(pattern_triangle, 596)
    ngt_tiling.add_pattern(pattern_triangle, 595)
    ngt_tiling.add_pattern(pattern_triangle, 706)
    ngt_tiling.add_pattern(pattern_triangle, 593)
    ngt_tiling.add_pattern(pattern_triangle, 669)
    ngt_tiling.add_pattern(pattern_triangle, 693)
    ngt_tiling.add_pattern(pattern_triangle, 690)
    ngt_tiling.add_pattern(pattern_square, 646)
    ngt_tiling.add_pattern(pattern_square, 726)
    ngt_tiling.add_pattern(pattern_square, 730)
    ngt_tiling.add_pattern(pattern_square, 550)
    ngt_tiling.add_pattern(pattern_square, 422)
    ngt_tiling.add_pattern(pattern_triangle, 756)
    ngt_tiling.add_pattern(pattern_triangle, 755)
    ngt_tiling.add_pattern(pattern_triangle, 760)
    ngt_tiling.add_pattern(pattern_triangle, 732)
    ngt_tiling.add_pattern(pattern_triangle, 762)
    ngt_tiling.add_pattern(pattern_square, 553)
    ngt_tiling.add_pattern(pattern_square, 431)
    ngt_tiling.add_pattern(pattern_square, 689)
    ngt_tiling.add_pattern(pattern_rhombus, 779)
    ngt_tiling.add_pattern(pattern_rhombus, 552)
    ngt_tiling.add_pattern(pattern_rhombus, 741)
    ngt_tiling.add_pattern(pattern_triangle, 759)
    ngt_tiling.add_pattern(pattern_triangle, 720)
    ngt_tiling.add_pattern(pattern_triangle, 696)
    ngt_tiling.add_pattern(pattern_triangle, 672)
    ngt_tiling.add_pattern(pattern_triangle, 794)
    ngt_tiling.add_pattern(pattern_square, 742)
    ngt_tiling.add_pattern(pattern_triangle, 787)
    ngt_tiling.add_pattern(pattern_triangle, 588)
    ngt_tiling.add_pattern(pattern_triangle, 576)
    ngt_tiling.add_pattern(pattern_triangle, 766)
    ngt_tiling.add_pattern(pattern_triangle, 767)
    ngt_tiling.add_pattern(pattern_triangle, 724)
    ngt_tiling.add_pattern(pattern_triangle, 361)
    ngt_tiling.add_pattern(pattern_rhombus, 838)
    ngt_tiling.add_pattern(pattern_rhombus, 340)
    ngt_tiling.add_pattern(pattern_rhombus, 580)
    ngt_tiling.add_pattern(pattern_triangle, 764)
    ngt_tiling.add_pattern(pattern_triangle, 789)
    ngt_tiling.add_pattern(pattern_triangle, 723)
    ngt_tiling.add_pattern(pattern_triangle, 260)
    ngt_tiling.add_pattern(pattern_triangle, 394)
    ngt_tiling.add_pattern(pattern_triangle, 395)
    ngt_tiling.add_pattern(pattern_triangle, 752)
    ngt_tiling.add_pattern(pattern_triangle, 751)
    ngt_tiling.add_pattern(pattern_triangle, 826)
    ngt_tiling.add_pattern(pattern_square, 841)
    ngt_tiling.add_pattern(pattern_square, 359)
    ngt_tiling.add_pattern(pattern_square, 870)
    ngt_tiling.add_pattern(pattern_triangle, 790)
    ngt_tiling.add_pattern(pattern_triangle, 763)
    ngt_tiling.add_pattern(pattern_triangle, 864)
    ngt_tiling.add_pattern(pattern_triangle, 888)
    ngt_tiling.add_pattern(pattern_triangle, 840)
    ngt_tiling.add_pattern(pattern_triangle, 893)
    ngt_tiling.add_pattern(pattern_triangle, 894)
    ngt_tiling.add_pattern(pattern_triangle, 399)
    ngt_tiling.add_pattern(pattern_triangle, 896)
    ngt_tiling.add_pattern(pattern_triangle, 897)
    ngt_tiling.add_pattern(pattern_triangle, 898)
    ngt_tiling.add_pattern(pattern_triangle, 388)
    ngt_tiling.add_pattern(pattern_triangle, 479)
    ngt_tiling.add_pattern(pattern_square, 907)
    ngt_tiling.add_pattern(pattern_square, 481)
    ngt_tiling.add_pattern(pattern_square, 931)
    ngt_tiling.add_pattern(pattern_square, 702)
    ngt_tiling.add_pattern(pattern_square, 747)
    ngt_tiling.add_pattern(pattern_square, 748)
    ngt_tiling.add_pattern(pattern_rhombus, 844)
    ngt_tiling.add_pattern(pattern_rhombus, 350)
    ngt_tiling.add_pattern(pattern_triangle, 874)
    ngt_tiling.add_pattern(pattern_square, 912)
    ngt_tiling.add_pattern(pattern_triangle, 964)
    ngt_tiling.add_pattern(pattern_triangle, 890)
    ngt_tiling.add_pattern(pattern_triangle, 889)
    ngt_tiling.add_pattern(pattern_triangle, 909)
    ngt_tiling.add_pattern(pattern_triangle, 956)
    ngt_tiling.add_pattern(pattern_triangle, 947)
    ngt_tiling.add_pattern(pattern_rhombus, 906)
    ngt_tiling.add_pattern(pattern_square, 985)
    ngt_tiling.add_pattern(pattern_triangle, 1002)
    ngt_tiling.add_pattern(pattern_triangle, 984)
    ngt_tiling.add_pattern(pattern_triangle, 975)
    ngt_tiling.add_pattern(pattern_triangle, 979)
    ngt_tiling.add_pattern(pattern_triangle, 915)
    ngt_tiling.add_pattern(pattern_square, 982)
    ngt_tiling.add_pattern(pattern_square, 1011)
    ngt_tiling.add_pattern(pattern_triangle, 955)
    ngt_tiling.add_pattern(pattern_triangle, 577)
    ngt_tiling.add_pattern(pattern_triangle, 859)
    ngt_tiling.add_pattern(pattern_triangle, 771)
    ngt_tiling.add_pattern(pattern_triangle, 952)
    ngt_tiling.add_pattern(pattern_triangle, 858)
    ngt_tiling.add_pattern(pattern_square, 1042)
    ngt_tiling.add_pattern(pattern_triangle, 1043)
    ngt_tiling.add_pattern(pattern_square, 1036)
    ngt_tiling.add_pattern(pattern_rhombus, 1049)
    ngt_tiling.add_pattern(pattern_triangle, 1052)
    ngt_tiling.add_pattern(pattern_square, 928)
    ngt_tiling.add_pattern(pattern_square, 972)
    ngt_tiling.add_pattern(pattern_square, 919)
    ngt_tiling.add_pattern(pattern_square, 922)
    ngt_tiling.add_pattern(pattern_square, 1007)
    ngt_tiling.add_pattern(pattern_square, 882)
    ngt_tiling.add_pattern(pattern_triangle, 1001)
    ngt_tiling.add_pattern(pattern_triangle, 941)
    ngt_tiling.add_pattern(pattern_triangle, 1079)
    ngt_tiling.add_pattern(pattern_triangle, 1020)
    ngt_tiling.add_pattern(pattern_triangle, 1021)
    ngt_tiling.add_pattern(pattern_triangle, 1073)
    ngt_tiling.add_pattern(pattern_triangle, 1069)
    ngt_tiling.add_pattern(pattern_triangle, 1068)
    ngt_tiling.add_pattern(pattern_triangle, 1067)
    ngt_tiling.add_pattern(pattern_triangle, 428)
    ngt_tiling.add_pattern(pattern_triangle, 877)
    ngt_tiling.add_pattern(pattern_triangle, 924)
    ngt_tiling.add_pattern(pattern_triangle, 948)
    ngt_tiling.add_pattern(pattern_triangle, 880)
    ngt_tiling.add_pattern(pattern_triangle, 636)
    ngt_tiling.add_pattern(pattern_triangle, 943)
    ngt_tiling.add_pattern(pattern_triangle, 493)
    ngt_tiling.add_pattern(pattern_triangle, 1025)
    ngt_tiling.add_pattern(pattern_triangle, 1024)
    ngt_tiling.add_pattern(pattern_square, 1100)
    ngt_tiling.add_pattern(pattern_rhombus, 1142)
    ngt_tiling.add_pattern(pattern_triangle, 944)
    ngt_tiling.add_pattern(pattern_rhombus, 933)
    ngt_tiling.add_pattern(pattern_triangle, 1096)
    ngt_tiling.add_pattern(pattern_triangle, 1072)
    ngt_tiling.add_pattern(pattern_triangle, 1076)
    ngt_tiling.add_pattern(pattern_triangle, 1106)
    ngt_tiling.add_pattern(pattern_rhombus, 1166)
    ngt_tiling.add_pattern(pattern_square, 1103)
    ngt_tiling.add_pattern(pattern_square, 1108)
    ngt_tiling.add_pattern(pattern_triangle, 1175)
    ngt_tiling.add_pattern(pattern_triangle, 1177)
    ngt_tiling.add_pattern(pattern_triangle, 1180)
    ngt_tiling.add_pattern(pattern_triangle, 1179)
    ngt_tiling.add_pattern(pattern_triangle, 1109)
    ngt_tiling.add_pattern(pattern_triangle, 959)
    ngt_tiling.add_pattern(pattern_square, 1111)
    ngt_tiling.add_pattern(pattern_square, 1112)
    ngt_tiling.add_pattern(pattern_square, 1196)
    ngt_tiling.add_pattern(pattern_square, 993)
    ngt_tiling.add_pattern(pattern_square, 1132)
    ngt_tiling.add_pattern(pattern_square, 1136)
    ngt_tiling.add_pattern(pattern_triangle, 1211)
    ngt_tiling.add_pattern(pattern_triangle, 1202)
    ngt_tiling.add_pattern(pattern_triangle, 1201)
    ngt_tiling.add_pattern(pattern_triangle, 1207)
    ngt_tiling.add_pattern(pattern_triangle, 1118)
    ngt_tiling.add_pattern(pattern_square, 1235)
    ngt_tiling.add_pattern(pattern_rhombus, 1229)
    ngt_tiling.add_pattern(pattern_triangle, 1246)
    ngt_tiling.add_pattern(pattern_triangle, 1245)
    ngt_tiling.add_pattern(pattern_triangle, 1242)
    ngt_tiling.add_pattern(pattern_triangle, 1241)
    ngt_tiling.add_pattern(pattern_square, 1251)
    ngt_tiling.add_pattern(pattern_square, 1255)
    ngt_tiling.add_pattern(pattern_rhombus, 1258)
    ngt_tiling.add_pattern(pattern_triangle, 1261)
    ngt_tiling.add_pattern(pattern_triangle, 1266)
    ngt_tiling.add_pattern(pattern_triangle, 1264)
    ngt_tiling.add_pattern(pattern_triangle, 1265)
    ngt_tiling.add_pattern(pattern_triangle, 1064)
    ngt_tiling.add_pattern(pattern_triangle, 1213)
    ngt_tiling.add_pattern(pattern_triangle, 1214)
    ngt_tiling.add_pattern(pattern_triangle, 1215)
    ngt_tiling.add_pattern(pattern_triangle, 1219)
    ngt_tiling.add_pattern(pattern_triangle, 1217)
    ngt_tiling.add_pattern(pattern_triangle, 1028)
    ngt_tiling.add_pattern(pattern_triangle, 1284)
    ngt_tiling.add_pattern(pattern_square, 1276)
    ngt_tiling.add_pattern(pattern_square, 1282)
    ngt_tiling.add_pattern(pattern_triangle, 1278)
    ngt_tiling.add_pattern(pattern_square, 1306)
    ngt_tiling.add_pattern(pattern_triangle, 1291)
    ngt_tiling.add_pattern(pattern_triangle, 1308)
    ngt_tiling.add_pattern(pattern_triangle, 1309)
    ngt_tiling.add_pattern(pattern_triangle, 1310)
    ngt_tiling.add_pattern(pattern_triangle, 1297)
    ngt_tiling.add_pattern(pattern_triangle, 1222)
    ngt_tiling.add_pattern(pattern_triangle, 1313)
    ngt_tiling.add_pattern(pattern_triangle, 1314)
    ngt_tiling.add_pattern(pattern_rhombus, 1330)
    ngt_tiling.add_pattern(pattern_rhombus, 1342)
    ngt_tiling.add_pattern(pattern_rhombus, 1294)
    ngt_tiling.add_pattern(pattern_rhombus, 1129)
    ngt_tiling.add_pattern(pattern_triangle, 1321)
    ngt_tiling.add_pattern(pattern_triangle, 1344)
    ngt_tiling.add_pattern(pattern_square, 1324)
    ngt_tiling.add_pattern(pattern_triangle, 1356)
    ngt_tiling.add_pattern(pattern_triangle, 1369)
    ngt_tiling.add_pattern(pattern_triangle, 1370)
    ngt_tiling.add_pattern(pattern_triangle, 1371)
    ngt_tiling.add_pattern(pattern_triangle, 1084)
    ngt_tiling.add_pattern(pattern_square, 1335)
    ngt_tiling.add_pattern(pattern_triangle, 1218)
    ngt_tiling.add_pattern(pattern_square, 1299)
    ngt_tiling.add_pattern(pattern_triangle, 1395)
    ngt_tiling.add_pattern(pattern_triangle, 1396)
    ngt_tiling.add_pattern(pattern_triangle, 1397)
    ngt_tiling.add_pattern(pattern_triangle, 1083)
    ngt_tiling.add_pattern(pattern_triangle, 1409)
    ngt_tiling.add_pattern(pattern_rhombus, 1411)
    ngt_tiling.add_pattern(pattern_triangle, 822)
    ngt_tiling.add_pattern(pattern_square, 1418)
    return ngt_tiling


def sigma_til(pattern_square, pattern_triangle, pattern_rhombus):
    """A sigma-approximant tiling canvas."""
    sigma = canvas(pattern_triangle)
    sigma.add_pattern(pattern_triangle, 3)
    sigma.add_pattern(pattern_square, 2)
    sigma.add_pattern(pattern_square, 1)
    sigma.add_pattern(pattern_triangle, 10)
    sigma.add_pattern(pattern_triangle, 12)
    sigma.add_pattern(pattern_square, 5)
    sigma.add_pattern(pattern_triangle, 14)
    sigma.add_pattern(pattern_square, 6)
    sigma.add_pattern(pattern_triangle, 8)
    return sigma


def approx1_til(pattern_square, pattern_triangle, pattern_rhombus):
    """A first-approximant tiling canvas."""
    approx1 = canvas(pattern_triangle)
    approx1.add_pattern(pattern_square, 2)
    approx1.add_pattern(pattern_square, 3)
    approx1.add_pattern(pattern_triangle, 5)
    approx1.add_pattern(pattern_triangle, 11)
    approx1.add_pattern(pattern_square, 13)
    approx1.add_pattern(pattern_square, 17)
    approx1.add_pattern(pattern_triangle, 25)
    approx1.add_pattern(pattern_square, 28)
    approx1.add_pattern(pattern_triangle, 20)
    approx1.add_pattern(pattern_triangle, 35)
    approx1.add_pattern(pattern_triangle, 24)
    approx1.add_pattern(pattern_triangle, 30)
    approx1.add_pattern(pattern_triangle, 23)
    approx1.add_pattern(pattern_triangle, 10)
    approx1.add_pattern(pattern_triangle, 6)
    approx1.add_pattern(pattern_triangle, 53)
    approx1.add_pattern(pattern_square, 1)
    approx1.add_pattern(pattern_triangle, 7)
    approx1.add_pattern(pattern_triangle, 62)
    approx1.add_pattern(pattern_triangle, 9)
    approx1.add_pattern(pattern_triangle, 60)
    return approx1


def approx2_til(pattern_square, pattern_triangle, pattern_rhombus):
    """A second-approximant tiling canvas."""
    approx2 = canvas(pattern_square)

    approx2.add_pattern(pattern_triangle, 1)
    approx2.add_pattern(pattern_triangle, 2)
    approx2.add_pattern(pattern_triangle, 3)
    approx2.add_pattern(pattern_triangle, 4)
    approx2.add_pattern(pattern_triangle, 13)
    approx2.add_pattern(pattern_triangle, 6)
    approx2.add_pattern(pattern_square, 12)
    approx2.add_pattern(pattern_square, 19)
    approx2.add_pattern(pattern_square, 9)
    approx2.add_pattern(pattern_square, 10)
    approx2.add_pattern(pattern_triangle, 30)
    approx2.add_pattern(pattern_triangle, 33)
    approx2.add_pattern(pattern_triangle, 32)
    approx2.add_pattern(pattern_triangle, 38)
    approx2.add_pattern(pattern_rhombus, 7)
    approx2.add_pattern(pattern_triangle, 24)
    approx2.add_pattern(pattern_triangle, 25)
    approx2.add_pattern(pattern_triangle, 26)
    approx2.add_pattern(pattern_rhombus, 63)
    approx2.add_pattern(pattern_triangle, 53)
    approx2.add_pattern(pattern_triangle, 67)
    approx2.add_pattern(pattern_triangle, 72)
    approx2.add_pattern(pattern_square, 69)
    approx2.add_pattern(pattern_square, 21)
    approx2.add_pattern(pattern_triangle, 54)
    approx2.add_pattern(pattern_triangle, 86)
    approx2.add_pattern(pattern_triangle, 37)
    approx2.add_pattern(pattern_triangle, 82)
    return approx2


def approx3_til(pattern_square, pattern_triangle, pattern_rhombus):
    """A third-approximant tiling canvas."""
    approx3 = canvas(pattern_square)
    approx3.add_pattern(pattern_triangle, 1)
    approx3.add_pattern(pattern_triangle, 2)
    approx3.add_pattern(pattern_triangle, 3)
    approx3.add_pattern(pattern_triangle, 4)
    approx3.add_pattern(pattern_triangle, 12)
    approx3.add_pattern(pattern_triangle, 13)
    approx3.add_pattern(pattern_square, 21)
    approx3.add_pattern(pattern_triangle, 24)
    approx3.add_pattern(pattern_rhombus, 29)
    approx3.add_pattern(pattern_rhombus, 18)
    approx3.add_pattern(pattern_triangle, 37)
    approx3.add_pattern(pattern_triangle, 32)
    approx3.add_pattern(pattern_square, 19)
    approx3.add_pattern(pattern_triangle, 35)
    approx3.add_pattern(pattern_triangle, 15)
    approx3.add_pattern(pattern_square, 7)
    approx3.add_pattern(pattern_square, 6)
    approx3.add_pattern(pattern_triangle, 46)
    approx3.add_pattern(pattern_triangle, 64)
    approx3.add_pattern(pattern_rhombus, 66)
    approx3.add_pattern(pattern_triangle, 56)
    approx3.add_pattern(pattern_triangle, 57)
    approx3.add_pattern(pattern_triangle, 59)
    approx3.add_pattern(pattern_triangle, 61)
    approx3.add_pattern(pattern_triangle, 25)
    approx3.add_pattern(pattern_triangle, 85)
    approx3.add_pattern(pattern_rhombus, 83)
    approx3.add_pattern(pattern_triangle, 91)
    return approx3


def bigapp_til(pattern_square, pattern_triangle, pattern_rhombus):
    """A larger approximant tiling canvas."""
    bigapprox = canvas(pattern_square)
    bigapprox.add_pattern(pattern_triangle, 1)
    bigapprox.add_pattern(pattern_triangle, 2)
    bigapprox.add_pattern(pattern_triangle, 3)
    bigapprox.add_pattern(pattern_triangle, 4)
    bigapprox.add_pattern(pattern_square, 6)
    bigapprox.add_pattern(pattern_square, 7)
    bigapprox.add_pattern(pattern_square, 9)
    bigapprox.add_pattern(pattern_triangle, 22)
    bigapprox.add_pattern(pattern_triangle, 23)
    bigapprox.add_pattern(pattern_triangle, 24)
    bigapprox.add_pattern(pattern_triangle, 18)
    bigapprox.add_pattern(pattern_triangle, 19)
    bigapprox.add_pattern(pattern_triangle, 20)
    bigapprox.add_pattern(pattern_triangle, 26)
    bigapprox.add_pattern(pattern_triangle, 27)
    bigapprox.add_pattern(pattern_triangle, 28)
    bigapprox.add_pattern(pattern_rhombus, 16)
    bigapprox.add_pattern(pattern_rhombus, 34)
    bigapprox.add_pattern(pattern_rhombus, 46)
    bigapprox.add_pattern(pattern_triangle, 55)
    bigapprox.add_pattern(pattern_triangle, 58)
    bigapprox.add_pattern(pattern_triangle, 59)
    bigapprox.add_pattern(pattern_triangle, 62)
    bigapprox.add_pattern(pattern_triangle, 63)
    bigapprox.add_pattern(pattern_triangle, 33)
    bigapprox.add_pattern(pattern_triangle, 43)
    bigapprox.add_pattern(pattern_triangle, 66)
    bigapprox.add_pattern(pattern_triangle, 67)
    bigapprox.add_pattern(pattern_triangle, 52)
    bigapprox.add_pattern(pattern_triangle, 75)
    bigapprox.add_pattern(pattern_triangle, 73)
    bigapprox.add_pattern(pattern_square, 69)
    bigapprox.add_pattern(pattern_square, 76)
    bigapprox.add_pattern(pattern_square, 78)
    bigapprox.add_pattern(pattern_rhombus, 102)
    bigapprox.add_pattern(pattern_square, 99)
    bigapprox.add_pattern(pattern_triangle, 88)
    bigapprox.add_pattern(pattern_triangle, 91)
    bigapprox.add_pattern(pattern_rhombus, 125)
    bigapprox.add_pattern(pattern_triangle, 118)
    bigapprox.add_pattern(pattern_triangle, 117)
    bigapprox.add_pattern(pattern_triangle, 110)
    bigapprox.add_pattern(pattern_triangle, 111)
    bigapprox.add_pattern(pattern_triangle, 121)
    bigapprox.add_pattern(pattern_triangle, 122)
    bigapprox.add_pattern(pattern_triangle, 114)
    bigapprox.add_pattern(pattern_square, 142)
    bigapprox.add_pattern(pattern_square, 141)
    bigapprox.add_pattern(pattern_square, 96)
    bigapprox.add_pattern(pattern_square, 97)
    bigapprox.add_pattern(pattern_triangle, 106)
    bigapprox.add_pattern(pattern_triangle, 169)
    bigapprox.add_pattern(pattern_triangle, 170)
    bigapprox.add_pattern(pattern_triangle, 164)
    bigapprox.add_pattern(pattern_triangle, 165)
    bigapprox.add_pattern(pattern_square, 37)
    bigapprox.add_pattern(pattern_square, 39)
    bigapprox.add_pattern(pattern_square, 90)
    bigapprox.add_pattern(pattern_square, 129)
    bigapprox.add_pattern(pattern_triangle, 162)
    bigapprox.add_pattern(pattern_triangle, 199)
    bigapprox.add_pattern(pattern_triangle, 200)
    bigapprox.add_pattern(pattern_triangle, 201)
    bigapprox.add_pattern(pattern_triangle, 189)
    bigapprox.add_pattern(pattern_triangle, 81)
    bigapprox.add_pattern(pattern_triangle, 150)
    bigapprox.add_pattern(pattern_triangle, 158)
    bigapprox.add_pattern(pattern_triangle, 160)
    bigapprox.add_pattern(pattern_rhombus, 204)
    bigapprox.add_pattern(pattern_triangle, 161)
    bigapprox.add_pattern(pattern_square, 84)
    bigapprox.add_pattern(pattern_square, 218)
    bigapprox.add_pattern(pattern_square, 173)
    bigapprox.add_pattern(pattern_triangle, 247)
    bigapprox.add_pattern(pattern_triangle, 245)
    bigapprox.add_pattern(pattern_triangle, 237)
    bigapprox.add_pattern(pattern_triangle, 238)
    bigapprox.add_pattern(pattern_triangle, 188)
    bigapprox.add_pattern(pattern_square, 126)
    bigapprox.add_pattern(pattern_triangle, 264)
    bigapprox.add_pattern(pattern_triangle, 266)
    bigapprox.add_pattern(pattern_triangle, 241)
    return bigapprox


def hexapp_til(pattern_square, pattern_triangle, pattern_rhombus):
    """A hexagonal-approximant tiling canvas."""
    hexapprox = canvas(pattern_rhombus)
    hexapprox.add_pattern(pattern_triangle, 1)
    hexapprox.add_pattern(pattern_triangle, 2)
    hexapprox.add_pattern(pattern_triangle, 3)
    hexapprox.add_pattern(pattern_triangle, 4)
    hexapprox.add_pattern(pattern_triangle, 7)
    hexapprox.add_pattern(pattern_triangle, 13)
    hexapprox.add_pattern(pattern_triangle, 18)
    hexapprox.add_pattern(pattern_triangle, 21)
    hexapprox.add_pattern(pattern_square, 15)
    hexapprox.add_pattern(pattern_square, 16)
    hexapprox.add_pattern(pattern_square, 9)
    hexapprox.add_pattern(pattern_square, 10)
    hexapprox.add_pattern(pattern_triangle, 22)
    hexapprox.add_pattern(pattern_triangle, 19)
    hexapprox.add_pattern(pattern_triangle, 35)
    hexapprox.add_pattern(pattern_triangle, 36)
    hexapprox.add_pattern(pattern_triangle, 30)
    hexapprox.add_pattern(pattern_triangle, 31)
    hexapprox.add_pattern(pattern_triangle, 43)
    hexapprox.add_pattern(pattern_triangle, 44)
    hexapprox.add_pattern(pattern_triangle, 38)
    hexapprox.add_pattern(pattern_triangle, 39)
    hexapprox.add_pattern(pattern_square, 47)
    hexapprox.add_pattern(pattern_square, 50)
    hexapprox.add_pattern(pattern_triangle, 76)
    hexapprox.add_pattern(pattern_triangle, 77)
    hexapprox.add_pattern(pattern_triangle, 78)
    hexapprox.add_pattern(pattern_triangle, 80)
    hexapprox.add_pattern(pattern_triangle, 81)
    hexapprox.add_pattern(pattern_triangle, 82)
    hexapprox.add_pattern(pattern_rhombus, 27)
    hexapprox.add_pattern(pattern_rhombus, 46)
    hexapprox.add_pattern(pattern_rhombus, 24)
    hexapprox.add_pattern(pattern_rhombus, 49)
    hexapprox.add_pattern(pattern_square, 93)
    hexapprox.add_pattern(pattern_square, 84)
    hexapprox.add_pattern(pattern_triangle, 122)
    hexapprox.add_pattern(pattern_triangle, 118)
    hexapprox.add_pattern(pattern_triangle, 120)
    hexapprox.add_pattern(pattern_triangle, 124)
    hexapprox.add_pattern(pattern_square, 52)
    hexapprox.add_pattern(pattern_square, 64)
    hexapprox.add_pattern(pattern_triangle, 56)
    hexapprox.add_pattern(pattern_triangle, 68)
    return hexapprox


def honeycomb_til(pattern_triangle):
    """A honeycomb approximant tiling canvas."""
    honeycomb = canvas(center_pattern=pattern_triangle)
    honeycomb.add_pattern(pattern_triangle, 1)
    honeycomb.add_pattern(pattern_triangle, 5)
    honeycomb.add_pattern(pattern_triangle, 9)
    honeycomb.add_pattern(pattern_triangle, 2)
    honeycomb.add_pattern(pattern_triangle, 11)
    honeycomb.add_pattern(pattern_triangle, 14)
    honeycomb.add_pattern(pattern_triangle, 17)
    honeycomb.add_pattern(pattern_triangle, 20)
    honeycomb.add_pattern(pattern_triangle, 23)
    honeycomb.add_pattern(pattern_triangle, 29)
    honeycomb.add_pattern(pattern_triangle, 6)
    honeycomb.add_pattern(pattern_triangle, 12)
    honeycomb.add_pattern(pattern_triangle, 24)
    honeycomb.add_pattern(pattern_triangle, 33)
    honeycomb.add_pattern(pattern_triangle, 39)
    honeycomb.add_pattern(pattern_triangle, 42)
    honeycomb.add_pattern(pattern_triangle, 45)
    honeycomb.add_pattern(pattern_triangle, 36)
    honeycomb.add_pattern(pattern_triangle, 8)
    honeycomb.add_pattern(pattern_triangle, 38)
    honeycomb.add_pattern(pattern_triangle, 48)
    honeycomb.add_pattern(pattern_triangle, 51)
    honeycomb.add_pattern(pattern_triangle, 54)
    honeycomb.add_pattern(pattern_triangle, 56)
    honeycomb.add_pattern(pattern_triangle, 59)
    honeycomb.add_pattern(pattern_triangle, 62)
    honeycomb.add_pattern(pattern_triangle, 66)
    honeycomb.add_pattern(pattern_triangle, 69)
    honeycomb.add_pattern(pattern_triangle, 72)
    honeycomb.add_pattern(pattern_triangle, 57)
    honeycomb.add_pattern(pattern_triangle, 81)
    return honeycomb
