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


def wien_parse_energy(structure:str="case.struct", energy:str="case.energy"):
    """
    wien2k case.energy and case.energyso parser
    :param energy: input file
    :return klist(kx4 array: ), energy(k x states x 2 array): 
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
    kpoint_lines = [[_[0:19], _[19:38], _[38:57], _[57:67], _[67:73], _[73:79], _[79:-1]] for _ in kpoint_lines]
    kpoint_lines = np.array(kpoint_lines, dtype=float)
    
    # add kpath 
    struct_recip = struct(structure).return_recip()
    kpath = klist_kpath(kpoint_lines[:,0:3], struct_recip)
    
    min_states = np.min(kpoint_lines[:,-2])

    # make sure the coordinates of k are always smaller than 1.0
    energy_lines = enelines[0:]
    energy_lines_index = list(set([_ for _ in range(beg_index, len(enelines))])^set(kpoint_line_index))
    energy_lines = [enelines[_] for _ in energy_lines_index]

    energy_lines = [[_[:12], _[12:-1]] for _ in energy_lines]
    energy_lines = np.array(energy_lines, dtype=float)
    energy_lines = energy_lines[energy_lines[:, 0]<=min_states]
    
    return np.hstack((kpoint_lines[:,:3], kpath.reshape(-1, 1))), energy_lines.reshape(-1,int(min_states), 2)[:, :, :]


def wien_merge_energy(energy_head:str="case.energy", energy_num:int=4, fpath:str="./"):
    """
    Merge the energy file from case.energy_num or case.energyso_num
    :param energy_head: string
    :param energy_num: int
    :write case.energy: write energy file
    """

    with open("%s/%s_%d" % (fpath, energy_head, 1), "r") as f:
        ene1lines = f.readlines()
    # find the begin line index
    beg_index = 0
    for line in ene1lines:
        if line[90:100] != "":
            beg_index += 1
        else:
            break
    with open("%s/%s" % (fpath, energy_head), "w") as f:
        f.writelines(ene1lines)
        for _ in range(energy_num-1):
            with open("%s/%s_%d" % (fpath, energy_head, _+2), "r") as ff:
                lines = ff.readlines()
            f.writelines(lines[beg_index:])


def wien_kgen(kend:str="0,0,0,0.1,0.1,0", kmesh:str="11-11-1", fpath:str="/Users/", fname:str="gen.klist"):
    """
    Generate arbitrary klist for effective mass calculation, please do not use --band 
    :param kend: string,    "0.5,0.5,0.5,0.01,0.01,0.01" --> the center of the kmesh cuboid with the length of the sides
    :                       "0.0,0.0,0.0,0.0,0.0,0.5  --> two ends of the klist line
    :param kmesh: int, 
    :param fpath: str, 
    :param fname: str, 
    :write fname: a wien2k DOS calculation klist
    """
    kend = [float(_) for _ in kend.split(",")]
    kmesh = [int(_) for _ in kmesh.split("-")] 

    # cube mesh
    if len(kmesh) == 3:
        x_ = np.linspace(kend[0]-kend[3]/2, kend[0]+kend[3]/2, num=kmesh[0])
        y_ = np.linspace(kend[1]-kend[4]/2, kend[1]+kend[4]/2, num=kmesh[1])
        z_ = np.linspace(kend[2]-kend[5]/2, kend[2]+kend[5]/2, num=kmesh[2])
        xx, yy, zz = np.meshgrid(x_, y_, z_, indexing='ij')
        x, y, z = xx.flatten(), yy.flatten(), zz.flatten()
        
        with open(fpath+"/"+fname, "w") as f:
            lineindex = 1
            for _ in range(int(kmesh[0]*kmesh[1]*kmesh[2])):
                xyz = convert_3decimal_to_4integer_fractioin(x[_], y[_], z[_])          
                f.write("%10d%10d%10d%10d%10d%5.1f\n" % (lineindex, xyz[0], xyz[1], xyz[2], xyz[3], 1))
                lineindex += 1
    # line mesh
    if len(kmesh) == 1:
        x_ = np.linspace(kend[0], kend[3], num=kmesh[0])
        y_ = np.linspace(kend[1], kend[4], num=kmesh[0])
        z_ = np.linspace(kend[2], kend[5], num=kmesh[0])

        with open(fpath+"/"+fname, "w") as f:
            lineindex = 1
            for _ in range(kmesh[0]):
                xyz = convert_3decimal_to_4integer_fractioin(x_[_], y_[_], z_[_])
                f.write("%10d%10d%10d%10d%10d%5.1f\n" % (lineindex, xyz[0], xyz[1], xyz[2], xyz[3], 1))
                lineindex += 1


def wien_fermi_surface(structure:str="case.struct", energy:str="ase.energy"):
    pass
