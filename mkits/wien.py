from mkits.globle import *

"""
:func gen_wien_born: generate the structure file for BORN effective charge calculation
"""

def gen_wien_born(structure: str, atom_list_index:list, direction:list, displacement:float = 0.01, fpath:str= "./"):
    """
    generate the 
    :param structure: string, input file
    :param atom_list: integer list from 1, select the atom to calculate, [1,2,3,4]
    :param direction: string list, select the direction to move, optional ["x", "y", "z"]
    :param displacement: float, 
    """
    
    for idx in atom_list_index:
        for direc in direction:
            for disp in [displacement, -displacement]:
                if direc == "x":
                    direc_idx = 0
                elif direc == "y":
                    direc_idx = 1
                elif direc == "z":
                    direc_idx = 2
                structure_class = struct(fpath+"/"+structure)
                structure_dict = structure_class.return_dict()
                new_pos = structure_dict["pos_frac"][idx-1][direc_idx] + disp
                if new_pos < 0:
                    structure_dict["pos_frac"][idx-1][direc_idx] = 1 + new_pos
                else:
                    structure_dict["pos_frac"][idx-1][direc_idx] = new_pos

                structure_class.update_struct_dict(structure_dict)
                structure_class.write_struct(fpath, "%s_%d_%s_%s_%.3f.struct" % (structure.replace(".struct", ""), idx, structure_dict["atoms_type"][idx-1], direc, disp))
                
