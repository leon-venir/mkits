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


def wien_parse_energy(energy:str="case.energy"):
    """
    wien2k case.energy and case.energyso parser
    :param energy: input file
    :return klist(kx3 array), energy(k x states x 2 array): 
    """
    with open(energy, "r") as f:
        enelines = f.readlines()
    
    # find the begin line index
    beg_index = 0
    for line in enelines:
        if line[90:100] != "":
            beg_index += 1
        else:
            break
    
    # line_index_of_1st_kpoints  the_index_of_1st_states_num
    numline = len(enelines)
    kpoint_line_index = []
    for _ in range(beg_index, numline):
        if len(enelines[_]) > 50: 
            kpoint_line_index.append(_)

    kpoint_lines = [enelines[_] for _ in kpoint_line_index]
    kpoint_lines = [_.split() for _ in kpoint_lines]
    kpoint_lines = np.array(kpoint_lines, dtype=float)
    min_states = np.min(kpoint_lines[:,-2])

    # make sure the coordinates of k are always smaller than 1.0
    energy_lines = np.loadtxt(energy, skiprows=beg_index, usecols=[0,1])
    energy_lines = energy_lines[energy_lines[:, 0]<=min_states]
    
    return kpoint_lines[:,:3], energy_lines.reshape(-1,int(min_states)+1, 2)[:, 1:, :]
    

def wien_fermi_surface(structure:str="case.struct", energy:str="ase.energy"):
    pass