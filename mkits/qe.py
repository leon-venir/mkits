from mkits.globle import *
import ase.io
import numpy as np
import sys


"""
:func qe_input_lines_parser: get 
:func qe_get_struct: return a standard output for QE input file
"""


def qe_input_lines_parser(lines:str):
    """
    &CONTROL
    &SYSTEM
    &ELECTRONS
    &IONS
    &CELL
    ATOMIC_SPECIES
    K_POINTS
    CELL_PARAMETERS
    ATOMIC_POSITIONS
    """

    lines_control = {
        "CONTROL_beg": -1,
        "CONTROL_end": -1,
        "SYSTEM_beg": -1,
        "SYSTEM_end": -1,
        "ELECTRONS_beg": -1,
        "ELECTRONS_end": -1,
        "IONS_beg": -1,
        "IONS_end": -1,
        "CELL_beg": -1,
        "CELL_end": -1,
        "ATOMIC_SPECIES_beg": -1,
        "ATOMIC_SPECIES_end": -1,
        "K_POINTS_beg": -1,
        "K_POINTS_end": -1,
        "CELL_PARAMETERS_beg": -1,
        "CELL_PARAMETERS_end": -1,
        "ATOMIC_POSITIONS_beg": -1,
        "ATOMIC_POSITIONS_end": -1
    }

    for i in range(len(lines)):
        if "&CONTROL" in lines[i]:
            lines_control["CONTROL_beg"] = i
            for j in range(i, len(lines)):
                if "/" in lines[j]:
                    lines_control["CONTROL_end"] = j
                    break
        elif "&SYSTEM" in lines[i]:
            lines_control["SYSTEM_beg"] = i
            for j in range(i, len(lines)):
                if "/" in lines[j]:
                    lines_control["SYSTEM_end"] = j
                    break
        elif "&ELECTRONS" in lines[i]:
            lines_control["ELECTRONS_beg"] = i
            for j in range(i, len(lines)):
                if "/" in lines[j]:
                    lines_control["ELECTRONS_end"] = j
                    break
        elif "&IONS" in lines[i]:
            lines_control["IONS_beg"] = i
            for j in range(i, len(lines)):
                if "/" in lines[j]:
                    lines_control["IONS_end"] = j
                    break
        elif "&CELL" in lines[i]:
            lines_control["CELL_beg"] = i
            for j in range(i, len(lines)):
                if "/" in lines[j]:
                    lines_control["CELL_end"] = j
                    break
        elif "ATOMIC_SPECIES" in lines[i]:
            lines_control["ATOMIC_SPECIES_beg"] = i
            for j in range(i, len(lines)):
                if lines[j].split():
                    pass
                else:
                    lines_control["ATOMIC_SPECIES_end"] = j
                    break
        elif "K_POINTS" in lines[i]:
            lines_control["K_POINTS_beg"] = i
            for j in range(i, len(lines)):
                if lines[j].split():
                    pass
                else:
                    lines_control["K_POINTS_end"] = j
                    break
        elif "CELL_PARAMETERS" in lines[i]:
            lines_control["CELL_PARAMETERS_beg"] = i
            for j in range(i, len(lines)):
                if lines[j].split():
                    pass
                else:
                    lines_control["CELL_PARAMETERS_end"] = j
                    break
        elif "ATOMIC_POSITIONS" in lines[i]:
            lines_control["ATOMIC_POSITIONS_beg"] = i
            for j in range(i, len(lines)):
                if lines[j].split():
                    pass
                else:
                    lines_control["ATOMIC_POSITIONS_end"] = j
                    break
    return lines_control


def qe_get_struct(file_inp:str):
    crystals = ase.io.read(filename=file_inp)
    ase.io.write(filename="./tmp", images=crystals, format="espresso-in")
    with open("./tmp", "r") as f:
        crystals_lines = f.readlines()
    lines_control = qe_input_lines_parser(crystals_lines)

    qe_struct = {
        "&SYSTEM": crystals_lines[lines_control["SYSTEM_beg"]:lines_control["SYSTEM_end"]+1],
        "ATOMIC_SPECIES": crystals_lines[lines_control["ATOMIC_SPECIES_beg"]:lines_control["ATOMIC_SPECIES_end"]+1],
        "K_POINTS": crystals_lines[lines_control["K_POINTS_beg"]:lines_control["K_POINTS_end"]+1],
        "CELL_PARAMETERS": crystals_lines[lines_control["CELL_PARAMETERS_beg"]:lines_control["CELL_PARAMETERS_end"]+1],
        "ATOMIC_POSITIONS": crystals_lines[lines_control["ATOMIC_POSITIONS_beg"]:lines_control["ATOMIC_POSITIONS_end"]+1],
        "CELL_UNITS": crystals_lines[lines_control["CELL_PARAMETERS_beg"]].split()[1],
        "atoms_num": crystals.get_atomic_numbers()
    }
    return qe_struct


def qe_get_upf(atom_numbers, upfpath:str="", wkdir:str="./"):
    """
    "return : list []
    """
    np.unique(atom_numbers)


def qe_geninput(calculation:str="scf", wpath:str="./", struct_inp:str="cu.cif", dryrun:bool=False):
    func_help = """
    Generate input file for QE calculations
    :param calculation: [scf, relax, vcrelax, conv_kmesh, conv_ecutwfc]
    """
    if dryrun:
        print(func_help)
        lexit("Show help with dryrun.")
    
    #
    if wpath == "./": wpath = os.path.abspath("./")
    wdir = wpath+"/"+calculation+"/" 


    qe_struct = qe_get_struct(file_inp=struct_inp)


