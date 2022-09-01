import numpy as np
from mkits.globle import *
from mkits.sysfc import *
from mkits.database import *

"""
:func fdmnes_gen_inp: generate
"""

def fdmnes_gen_inp(struct_inp, fpath="./", index=0, params=""):
    """
    generate the input file for FDMNES simulation
    :param struct_inp: the xcrystal
    :param param: additional parameters, conpound=TiS2,element=Ti 
                : element       ->  
                : radius        -> 
                : method        -> [FDMNES, GREEN]
    """

    struct_inp = struct(struct_inp).return_dict()
    lattice = lattice_conversion(struct_inp["lattice"])
    compound = "".join([i for j in [list(atom) for atom in zip(struct_inp["atoms_type"], [str(_) for _ in struct_inp["atoms_num"]])] for i in j])
    atom_list = np.array([])
    for _ in range(len(struct_inp["atoms_type"])):
        atom_list = np.hstack((atom_list, [symbol_map[struct_inp["atoms_type"][_]]]*int(struct_inp["atoms_num"][_])))
    #default setting to write
    param_default = {
        "Radius": "7.0",
        "Green": "n",
        "SCF": "y",
        "Quadrupole": "y",
        "Convolution": "y"
    }
    if params:
        params = parser_inputpara(params)
    else:
        params = {}
    param_default.update(params)
    param_default.pop("element", None)
    
    if "element" in params:
        param_default["Z_Absorber"] = symbol_map[params["element"]]
    else:
        param_default["Z_Absorber"] = str(int(atom_list[0]))
        params["element"] = atom_data[int(atom_list[0])][1]
    
    fdmnes_inp_head = """! Fdmnes indata file
! Calculation for the %s K-edge in %s
! Finite difference method calculation with convolution

 Filout
./out/%s_%s

 Range                               ! Energy range of calculation (eV)
  -8. 0.5  10. 1.  18. 2. 60.        ! first energy, step, intermediary energy, step ..., last energy
 !-8. 0.2  13. 0.5 18. 1. 50. 2 120. ! first energy, step, intermediary energy, step ..., last energy    

 Crystal
"""

    with open(fpath+"/"+compound+"_inp_"+str(index)+".txt", "w") as f:
        f.write(fdmnes_inp_head % (params["element"], compound, compound, str(index)))
        f.write("%15.10f%15.10f%15.10f%6.1f%6.1f%6.1f\n" % (lattice[0], lattice[1], lattice[2], lattice[3], lattice[4], lattice[5]))
        np.savetxt(f, np.hstack((atom_list.reshape(-1,1), struct_inp["pos_frac"])), fmt="%3d%15.10f%15.10f%15.10f")
        f.write("\n")
        # write params
        for key in param_default.keys():
            if param_default[key] == "n":
                pass
            elif param_default[key] == "y":
                f.write(key+"\n\n")
            else:
                f.write(key+"\n")
                f.write(str(param_default[key])+"\n\n")
        f.write("END\n")