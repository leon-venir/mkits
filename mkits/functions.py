# -*- coding: utf-8 -*-
"""General-purpose utilities shared across mkits modules.

mkits is used as an imported library (Jupyter, a future web backend, and
eventually a CLI), so it must never call sys.exit() itself -- that would
be fatal to a long-running host process (eg a web worker). Error
conditions raise MkitsError instead; whichever front-end wraps mkits
decides what to do with it (print a traceback, return an HTTP error,
convert to a CLI exit code, ...). Likewise, this module only obtains a
logger -- it never calls logging.basicConfig(), since configuring
handlers/output location is the host application's decision, not a
library's.
"""

import os
import math
import logging

import numpy as np


logger = logging.getLogger(__name__)


class MkitsError(Exception):
    """Raised for user-facing input/parsing errors across mkits."""


def write2log(message):
    """Log an informational message."""
    logger.info(message)


def rmspace(string):
    """Remove all space characters from a string."""
    return "".join(c for c in string if c != " ")


def listcross(list1, list2):
    """
    Expand list1 by repeating each entry according to the counts in list2.
    len(list1) must equal len(list2).

    Example
    -------
    listcross(["Ge", "Te"], [2, 5])
        -> ["Ge", "Ge", "Te", "Te", "Te", "Te", "Te"]
    """
    if not isinstance(list1, list) or not isinstance(list2, (list, np.ndarray)):
        raise MkitsError("The input lists are not list type.")
    try:
        num_list = [int(i) for i in list2]
    except ValueError:
        raise MkitsError("Wrong numeric list.")
    crosslist = []
    for i in range(len(num_list)):
        crosslist += [list1[i]] * num_list[i]
    return crosslist


def hstack_append_list(list1, list2):
    """
    Append list2 to list1 column-wise, ie horizontally stack two nested lists.

    list1                  list2               result
    [["1","2"],        +  [["8", "9"],    ->  [["1","2","8","9"],
     ["3","4"],             ["10","11"]]        ["3","4","10","11"],
     ["5","7"]]                                 ["5","7"]]
    """
    if isinstance(list1[0], list) and isinstance(list2[0], list):
        if len(list1) >= len(list2):
            for i in range(len(list2)):
                list1[i] = list1[i] + list2[i]
        else:
            for i in range(len(list1)):
                list1[i] = list1[i] + list2[i]
            for i in range(len(list1), len(list2)):
                list1.append(list2[i])
    else:
        raise MkitsError("Input is not a nested list.")
    return list1


def convert_high2writeablelist(inp, endstr="\n"):
    """Flatten a 2D array of strings into a list of joined, newline-terminated lines."""
    if isinstance(inp, list):
        inp = np.array(inp)

    if inp.ndim == 1:
        inp = inp.reshape(1, len(inp))
    elif inp.ndim != 2:
        raise MkitsError("convert_high2writeablelist does not support arrays with more than 2 dimensions.")

    return ["".join(inp[i, :]) + endstr for i in range(inp.shape[0])]


def convert_array2strlist(inp, fmt="{:20.10f}"):
    """Convert a numeric array or list into a nested list of formatted strings."""
    if isinstance(inp, list):
        inp = np.array(inp)

    if inp.ndim == 1:
        inp = inp.reshape(1, len(inp))
    elif inp.ndim != 2:
        raise MkitsError("convert_array2strlist does not support arrays with more than 2 dimensions.")

    return [[fmt.format(inp[i, j]) for j in range(inp.shape[1])] for i in range(inp.shape[0])]


def parser_inputpara(inputstring):
    """
    Parse a compact "key:value,key:value" (or "key=value key=value") string
    into a dictionary.

    Parameters
    ----------
    inputstring: str
        Key-value pairs separated by "," (or ";" / " "), with the key and
        the value separated by ":" (or "=").

    Returns
    -------
    dict
    """
    input_dict = {}

    if "," in inputstring:
        separator_outkey = ","
    elif ";" in inputstring:
        separator_outkey = ";"
    else:
        separator_outkey = " "

    if inputstring and inputstring[0] == separator_outkey:
        inputstring = inputstring[1:]
    if inputstring and inputstring[-1] == separator_outkey:
        inputstring = inputstring[:-1]

    if ":" in inputstring:
        separator_inkey = ":"
    elif "=" in inputstring:
        separator_inkey = "="
    else:
        raise MkitsError("Cannot find a key-value separator (':' or '=') in the input string.")

    try:
        for item in inputstring.split(separator_outkey):
            key, value = item.split(separator_inkey)
            input_dict[key] = value
    except ValueError:
        raise MkitsError("Error: make sure the input parameter string follows: key1:para1,key2:para2")

    return input_dict


def apply_dynrange(struct_obj, dynrange_str):
    """
    Parse a "xmin=..,xmax=..,fix=O,move=Ti"-style range string and apply it
    to `struct_obj` via its add_dyn method. Shared by mkits.vasp and
    mkits.qe so both accept the same selective-dynamics / atom-constraint
    syntax.
    """
    _range = parser_inputpara(dynrange_str)
    _fix = _range.pop("fix", "none")
    _move = _range.pop("move", "none")
    _bounds = {"xmin": -1e8, "xmax": 1e8, "ymin": -1e8, "ymax": 1e8, "zmin": -1e8, "zmax": 1e8}
    _bounds.update({_k: float(_v) for _k, _v in _range.items()})
    struct_obj.add_dyn(fix=_fix, move=_move, **_bounds)


def write_runsh(runsh, cmd):
    """Write a shell run script and make it executable."""
    with open(runsh, "w", newline="\n") as f:
        f.write(cmd)
    os.chmod(runsh, 0o775)


def round_even_odd(num, even_odd):
    """
    Round a number to the nearest even integer, nearest odd integer, or
    just the nearest integer.

    Parameters
    ----------
    num: float or int
    even_odd: int
        0 -> nearest even, 1 -> nearest odd, -1 -> plain round.

    Returns
    -------
    int
    """
    if even_odd == 0:
        return int(num) + 1 if math.modf(num / 2)[0] >= 0.5 else int(num)
    elif even_odd == 1:
        num -= 1
        return int(num) + 2 if math.modf(num / 2)[0] >= 0.5 else int(num) + 1
    elif even_odd == -1:
        return round(num)
    else:
        raise MkitsError("Error in even_odd, use one of [0 (even), 1 (odd), -1 (round only)].")


def parse_inputfile(lines, comment_sym="#", assign_sym="=", seprate_sym=","):
    """
    Parse a list of "key<assign_sym>value" lines (eg a Fortran namelist body)
    into a dictionary, stripping comments and expanding multiple
    assignments per line.

    Parameters
    ----------
    lines: list of str
    comment_sym: str
        Comment marker; everything after it on a line is discarded.
    assign_sym: str
        Key/value separator.
    seprate_sym: str
        Separator between multiple assignments on the same line.

    Returns
    -------
    dict
    """
    keywords = {}

    for i in range(len(lines)):
        lines[i] = lines[i].strip()
        if comment_sym in lines[i]:
            lines[i] = lines[i][:lines[i].index(comment_sym)]

    newlines = []
    for line in lines:
        if seprate_sym in line:
            newlines += line.split(seprate_sym)
        else:
            newlines.append(line)

    for line in newlines:
        if not line or assign_sym not in line:
            continue
        assign_sym_index = line.index(assign_sym)
        key = line[:assign_sym_index].strip()
        value = line[assign_sym_index + 1:].strip()
        if key:
            keywords[key] = value

    return keywords


def lattice_conversion(give_lattice):
    """
    Convert between a 3x3 Cartesian lattice matrix and 6 lattice parameters.

    Mode 1: 3x3 matrix -> [a, b, c, alpha, beta, gamma]
    Mode 2: [a, b, c, alpha, beta, gamma] -> 3x3 matrix
            (standard convention: a along x, b in the xy-plane)

    Parameters
    ----------
    give_lattice: (3, 3) or (6,) array-like

    Returns
    -------
    numpy array
    """
    give_lattice = np.array(give_lattice, dtype=float)

    if give_lattice.shape == (3, 3):
        lengths = np.linalg.norm(give_lattice, axis=1)
        a, b, c = lengths[0], lengths[1], lengths[2]

        cos_alpha = np.clip(np.dot(give_lattice[1], give_lattice[2]) / (b * c), -1.0, 1.0)
        cos_beta = np.clip(np.dot(give_lattice[0], give_lattice[2]) / (a * c), -1.0, 1.0)
        cos_gamma = np.clip(np.dot(give_lattice[0], give_lattice[1]) / (a * b), -1.0, 1.0)

        alpha = np.degrees(np.arccos(cos_alpha))
        beta = np.degrees(np.arccos(cos_beta))
        gamma = np.degrees(np.arccos(cos_gamma))

        return np.array([a, b, c, alpha, beta, gamma])

    elif give_lattice.size == 6:
        params = give_lattice.flatten()
        a, b, c = params[0], params[1], params[2]
        alpha, beta, gamma = np.radians(params[3:6])

        # standard convention: a along x-axis, b in the xy-plane
        val = (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
        v_a = [a, 0.0, 0.0]
        v_b = [b * np.cos(gamma), b * np.sin(gamma), 0.0]
        cz = np.sqrt(max(1.0 - np.cos(beta) ** 2 - val ** 2, 0.0))
        v_c = [c * np.cos(beta), c * val, c * cz]

        return np.array([v_a, v_b, v_c])

    else:
        raise MkitsError("Lattice must be a 3x3 matrix or a 6-parameter array.")


def frac2cart(cart_lattice, fraction_pos):
    """
    Convert fractional coordinates to Cartesian coordinates.
    r_cart = frac_pos @ lattice_matrix (lattice vectors stored as rows).

    Parameters
    ----------
    cart_lattice: (3, 3) array-like
    fraction_pos: (3,) or (n, 3) array-like

    Returns
    -------
    numpy array with the same shape as fraction_pos.
    """
    lattice = np.array(cart_lattice, dtype=float)
    pos = np.array(fraction_pos, dtype=float)
    return pos @ lattice


def cart2frac(lattice, cart_pos):
    """
    Convert Cartesian coordinates to fractional coordinates.
    r_frac = cart_pos @ inverse(lattice_matrix).

    Parameters
    ----------
    lattice: (3, 3) array-like or 6 lattice parameters
    cart_pos: (3,) or (n, 3) array-like

    Returns
    -------
    numpy array with the same shape as cart_pos.
    """
    lattice = np.array(lattice, dtype=float)

    if lattice.shape == (6,) or lattice.size == 6:
        lattice = lattice_conversion(lattice)

    try:
        inv_lattice = np.linalg.inv(lattice)
    except np.linalg.LinAlgError:
        raise MkitsError("Lattice matrix is singular (zero volume), check the lattice parameters.")

    return np.asarray(cart_pos) @ inv_lattice


def vector_angle(vector1, vector2, unit="deg"):
    """
    Angle between two vectors, via arccos of their normalized dot product.

    Parameters
    ----------
    vector1, vector2: array-like
    unit: str
        "deg" or "rad".

    Returns
    -------
    float
    """
    _unit1 = np.asarray(vector1, dtype=float) / np.linalg.norm(vector1)
    _unit2 = np.asarray(vector2, dtype=float) / np.linalg.norm(vector2)
    _angle = np.arccos(np.clip(np.dot(_unit1, _unit2), -1.0, 1.0))
    return np.rad2deg(_angle) if unit == "deg" else _angle


def vector_angle_cclockwise(vector1, vector2, unit="deg"):
    """
    Counterclockwise angle (in [0, 2*pi) / [0, 360)) from vector2 to
    vector1, via the difference of their atan2 angles -- unlike
    vector_angle's undirected arccos-based angle (always in [0, pi]),
    this is signed/directional and only defined for 2D vectors.

    Parameters
    ----------
    vector1, vector2: (2,) array-like
    unit: str
        "deg" or "rad".

    Returns
    -------
    float
    """
    _v1 = list(vector1)
    _v2 = list(vector2)
    if len(_v1) != 2 or len(_v2) != 2:
        raise MkitsError("vector_angle_cclockwise only supports 2D vectors.")

    _angle1 = np.arctan2(_v1[1], _v1[0])
    _angle2 = np.arctan2(_v2[1], _v2[0])
    _angle = (_angle1 - _angle2) % (2 * np.pi)
    return np.rad2deg(_angle) if unit == "deg" else _angle


def group_atoms_by_symbol(symbols, positions):
    """
    Group atoms by element symbol in alphabetical order (NOT atomic-number
    order -- matches how mkits.builder's layered-structure generator sorts
    a POSCAR species block when it has no mkits.database.symbol_map
    ordering to follow), reordering `positions`'s rows to match.

    Parameters
    ----------
    symbols: (n,) array-like of str
    positions: (n, k) array-like
        Rows reordered to match the grouped symbols.

    Returns
    -------
    (symbols_sorted, positions_sorted, unique_symbols, counts)
        symbols_sorted: (n,) grouped symbols
        positions_sorted: (n, k) positions reordered to match
        unique_symbols, counts: (n_unique,) alphabetically sorted
    """
    _symbols = np.asarray(symbols)
    _sort_idx = _symbols.argsort()
    _symbols_sorted = _symbols[_sort_idx]
    _positions_sorted = np.asarray(positions)[_sort_idx]

    _unique_symbols = np.array(sorted(set(_symbols_sorted.tolist())))
    _counts = np.array([np.count_nonzero(_symbols_sorted == _s) for _s in _unique_symbols])
    return _symbols_sorted, _positions_sorted, _unique_symbols, _counts
