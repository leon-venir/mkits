
import numpy as np
import itertools as it
import json
from pathlib import Path
import spglib as spg
from mkits.globle import *
from mkits.database import *
from mkits.sysfc import *
from mkits.vasp import *


def get_node_type(node_type):
    if node_type == "trigonal":
        node_chains = np.array([[0.333333333, 0.666666667, 0], [0.666666667, 0.333333333, 0], [2,2,2]])
    return node_chains


def gen_nodes_frac(grid_num, node_type):
    # return the fractional coordinates of struct
    # [
    #   [0.33333, 0.66667, 0, 0.33333, 0.66666, 0] x -> 
    #   [0.66667, 0.33333, 0, 0.66667, 0.33333, 0] y
    #   [2,     , 2      , 2, 2      , 2      , 2] z -> if z is not fixed, then 2
    # ]
    node_chains_single = get_node_type(node_type)
    grid_per_node = len(node_chains_single[0, :])
    node_num = int(grid_num/grid_per_node)
    node_chains = node_chains_single
    for i in range(1, node_num):
        node_chains = np.hstack((node_chains, node_chains_single))
    return node_chains
#print(gen_nodes_frac(9, "trigonal"))


def randomization_block(node_num, block_num):
    # node_num=4 (int)
    # block_num=[2, 3, 5, 7, 9] (int)
    combinations_w_re = list(it.combinations_with_replacement(block_num, node_num))
    #print(combinations_w_re)
    # permitations
    results = []
    for i in range(len(combinations_w_re)):
        tmp = list(it.permutations(combinations_w_re[i], len(combinations_w_re[i])))
        tmp = list(set(tmp))+["d"]
        for j in range(len(tmp)-1):
            # replace the reversed item: (1,2,3) & (3,2,1)
            if tmp[j][::-1] in tmp[j+1:]:
                tmp[j] = "d"
            # replace the translation: (1,2,3) & (2,3,1) (3,1,2)
            translation_tmp = tmp[j]+tmp[j]
            for k in range(1, len(tmp[j])):
                if tmp[j] == "d":
                    break
                elif translation_tmp[k:k+len(tmp[j])] in tmp[j+1:]:
                    tmp[j] = "d"
                    break
            # delete block 3&9 (3,3,3) (9.01,9.01,9.01) (9.02,9.02,9.02)
            if node_num == 3 and tmp[j] != "d":
                if tmp[j][0] == tmp[j][1] and tmp[j][0] == tmp[j][2] and (tmp[j][0] == 3 or tmp[j][0] == 9.01 or tmp[j][0] == 9.02):
                    tmp[j] = "d"
            if node_num == 6 and tmp[j] != "d":
                if (tmp[j][0] == 3 or tmp[j][0] == 9.01 or tmp[j][0] == 9.02) and tmp[j][0] == tmp[j][1] and tmp[j][0] == tmp[j][2] and tmp[j][0] == tmp[j][3] and tmp[j][0] == tmp[j][4] and tmp[j][0] == tmp[j][5]:
                    tmp[j] = "d"
            if node_num == 9 and tmp[j] != "d":
                if (tmp[j][0] == 3 or tmp[j][0] == 9.01 or tmp[j][0] == 9.02) and tmp[j][0] == tmp[j][1] and tmp[j][0] == tmp[j][2] and tmp[j][0] == tmp[j][3] and tmp[j][0] == tmp[j][4] and tmp[j][0] == tmp[j][5] and tmp[j][0] == tmp[j][6] and tmp[j][0] == tmp[j][7] and tmp[j][0] == tmp[j][8]:
                    tmp[j] = "d"

        results += list(filter(("d").__ne__, tmp))
    tmp = [round(sum(i)%3) for i in results]
    tmp = [i for i,val in enumerate(tmp) if val==0]
    #print([results[i] for i in tmp])
    return [results[i] for i in tmp] 


def struct_write_abs(chain, block_num_type, block_bond_c_axis, node_type):
    #                                                                               gap 
    #block_num_type = {9.005: ['Te', 'Ge', 'Te', 'Sb', 'Te', 'Sb', 'Te', 'Ge', 'Te', 'Sb', 'Te']}
    #block_bond_c_axis = {9.005: [1.393991367, 2.199894149, 1.723797617, 2.01772035, 2.01772035, 1.723797617, 2.199894149, 1.393991367, 2.583808499, 4.264490915]}
    """
    {
    #.     bond insid, gap with Sb, gap with Te
    2.001: [1.552877096, 2.362777542, 2.739889686],
    5.002: [1.712154812, 1.994616362, 1.994616362, 1.712154812, 2.632587528, 3.908657866],
    7.003: [1.699051612, 2.035772228, 1.696968734, 1.696968734, 2.035772228, 1.699051612, 2.710683844, 4.224055311],
    9.004: [1.695489891, 2.049167754, 1.669882087, 1.747324809, 1.747324809, 1.669882087, 2.049167754, 1.695489891, 3.13253256, 4.600418979],
    9.005: [1.393991367, 2.199894149, 1.723797617, 2.01772035, 2.01772035, 1.723797617, 2.199894149, 1.393991367, 2.583808499, 4.264490915]
}
    """
    #node_type = "trigonal"


    chain_plus = chain + [chain[0]]
    atoms_chain = []
    for i in chain:
        atoms_chain += block_num_type[i][:-2]
    #print(atoms_chain)

    z_cord_abs_diff = []
    for i in range(len(chain)):
        z_cord_abs_diff += block_bond_c_axis[chain[i]][:-2]
        if chain_plus[i+1] > 2.1:
            z_cord_abs_diff.append(block_bond_c_axis[chain[i]][-1])
        else:
            z_cord_abs_diff.append(block_bond_c_axis[chain[i]][-2])

    #print(z_cord_abs_diff)
    z_cord_abs = [np.sum(z_cord_abs_diff[:i+1]) for i in range(len(z_cord_abs_diff))]
    z_cord_abs = [0] + z_cord_abs
    #print(z_cord_abs)
    z = z_cord_abs[-1]
    pos_frac_z = np.array(z_cord_abs[:-1])/z
    #print(pos_frac_z)

    total_grid = round(sum(chain))
    grid_pos = gen_nodes_frac(total_grid, node_type) # results fractional coordinates
    grid_pos[2,:] = np.array(pos_frac_z[:])
    #print(grid_pos)

    lattice = np.array([[ 4.3005422627252692, 0.0000000000000000, 0.0000000000000000],
                        [-2.1502711314087866, 3.7243788496048982, 0.0000000000000000],
                        [ 0.0000000000000000, 0.0000000000000000, z]])

    atoms_sort_idx = np.array(atoms_chain).argsort()
    atoms_chain = np.array(atoms_chain)[atoms_sort_idx]
    grid_pos = grid_pos[:, atoms_sort_idx]
    atom_type = np.array(list(set(atoms_chain)))
    atoms_num = np.array([])
    for i in list(set(atoms_chain)):
        atoms_num = np.append(atoms_num, atoms_chain.tolist().count(i))
    atom_type_idx = atom_type.argsort()
    atom_type = atom_type[atom_type_idx]
    atoms_num = atoms_num[atom_type_idx]

    poscar_dict = {
        "title"      : "",
        "ratio"      : 1.0,
        "atoms_type" : [],
        "atoms_num"  : [],
            
    }
    poscar_dict["atoms_type"] = atom_type.tolist()
    poscar_dict["pos_frac"] = grid_pos.T
    poscar_dict["atoms_num"] = atoms_num
    poscar_dict["lattice"] = lattice
    #print(poscar_dict)
    return poscar_dict



def struct_write(chain, block_num_type, node_type):
    # chain: (2.01, 5.0, 2.01); node_num: 3; node_type: trigonal
    total_grid = round(sum(chain))
    grid_pos = gen_nodes_frac(total_grid, node_type) # results fractional coordinates
    atoms_chain = block_num_type[chain[0]]
    for i in range(1, len(chain)):
        atoms_chain = block_num_type[chain[i]] + atoms_chain

    # write z coordinates bond:slab=4:5
    z_cord_idx = [int(i) for i in chain] # (2, 5, 2)
    #print(z_cord_idx)
    spli_per = 1/(4*(sum(z_cord_idx)-len(z_cord_idx))+5*len(z_cord_idx))
    #print(spli_per)
    z_cord = [0]
    for i in range(len(z_cord_idx)):
        z_cord = z_cord + [z_cord[-1]+spli_per*5]
        for j in range(z_cord_idx[i]-1):
            z_cord = z_cord + [z_cord[-1]+spli_per*4]
    grid_pos[2,:] = np.array(z_cord[1:])
    # lattice c
    z = 29.2422161102/60*(4*(sum(z_cord_idx)-len(z_cord_idx))+5*len(z_cord_idx))
    lattice = np.array([[ 4.3005422627252692, 0.0000000000000000, 0.0000000000000000],
                        [-2.1502711314087866, 3.7243788496048982, 0.0000000000000000],
                        [ 0.0000000000000000, 0.0000000000000000, z]])

    atoms_sort_idx = np.array(atoms_chain).argsort()
    atoms_chain = np.array(atoms_chain)[atoms_sort_idx]
    grid_pos = grid_pos[:, atoms_sort_idx]
    atom_type = np.array(list(set(atoms_chain)))
    atoms_num = np.array([])
    for i in list(set(atoms_chain)):
        atoms_num = np.append(atoms_num, atoms_chain.tolist().count(i))
    atom_type_idx = atom_type.argsort()
    atom_type = atom_type[atom_type_idx]
    atoms_num = atoms_num[atom_type_idx]

    poscar_dict = {
        "title"      : "",
        "ratio"      : 1.0,
        "atoms_type" : [],
        "atoms_num"  : [],
        
    }
    poscar_dict["atoms_type"] = atom_type.tolist()
    poscar_dict["pos_frac"] = grid_pos.T
    poscar_dict["atoms_num"] = atoms_num
    poscar_dict["lattice"] = lattice
    poscar_dict["calculator"] = "vasp"
    
    #return poscar_dict
    return struct(poscar_dict)


def gen_struct_fix_block(node_num, node_type, block_type, bond_abs_c=False):
    # node_num=2,3,4,5,6(str), node_type=trigonal(str), 
    # block_num=2.01,2.02,3,5,7,9.01,9.02(str)
    # block_type=bi-bi,te-te,te-bi-a,a-k-a-k-a,a-k-a-k-a-k-a,a-k-a-k-a-k-a-k-a
    # kation=Pb,Bi(str)
    # anion=Te(str)

    # get parameter
    node_num = [int(i) for i in node_num.split(",")]
    block_type = block_type.split(",")
    # from block_type -> block_num
    # block_num = [float(i) for i in block_num.split(",")]
    block_num = [float(i.count("-"))+1 for i in block_type]
    block_num = np.array(block_num) + np.linspace(0.001, 0.001*len(block_num), len(block_num))
    block_num_type = {}
    for i in range(len(block_num)):
        block_num_type[block_num[i]] = block_type[i].split("-")
    #print(block_num_type)
    
    chains = randomization_block(node_num[0], block_num)

    for chain in chains:
        if bond_abs_c:
            # add -Sb-Te 
            block_num_type = {2.001: ['Sb', 'Sb', 'Sb', 'Te'], 5.002: ['Te', 'Sb', 'Te', 'Sb', 'Te', 'Sb', 'Te'], 7.003: ['Te', 'Sb', 'Te', 'Ge', 'Te', 'Sb', 'Te', 'Sb', 'Te'], 9.004: ['Te', 'Sb', 'Te', 'Ge', 'Te', 'Ge', 'Te', 'Sb', 'Te', 'Sb', 'Te'], 9.005: ['Te', 'Ge', 'Te', 'Sb', 'Te', 'Sb', 'Te', 'Ge', 'Te', 'Sb', 'Te']}
            poscar = struct_write_abs(list(chain), block_num_type, bond_abs_c, node_type)
        else:
            poscar = struct_write(chain, block_num_type, node_type)
            
        filename1 = "struct_node%s_block" % str(node_num[0])
        filename2 = "_".join([str(i) for i in chain])
        try:
            Path("./layered/chains%s/" % str(node_num[0])).mkdir(parents=True, exist_ok=True)
        except:
            lexit("Cannot create layered folder, pls check it.")
        poscar.get_refined_struct()
        poscar.write_struct("./layered/chains%s/" % str(node_num[0]), filename1+filename2+".vasp")
    
    ch = open("./layered/chains%s/chains%s.data" % (str(node_num[0]), str(node_num[0])), "w")
    ch.write("Total structures: "+str(len(chains))+"\n")
    ch.write("========== Chains ==========\n"+str(chains)+"\n")
    ch.write("========== Blocks ==========\n"+json.dumps(block_num_type)+"\n")
    ch.close()

def gen_struct_random():
    x = 1 

def test_struct():
    #gen_struct_fix_block('3', 'trigonal', "Bi-Bi,Te-Te,Te-Bi-Te,Te-Bi-Te-Bi-Te,Te-Bi-Te-Pb-Te-Bi-Te,Te-Bi-Te-Pb-Te-Pb-Te-Bi-Te,Te-Pb-Te-Bi-Te-Bi-Te-Pb-Te")
    #gen_struct_fix_block('3', 'trigonal', "Bi-Bi,Te-Bi-Te-Pb-Te-Pb-Te-Bi-Te,Te-Pb-Te-Bi-Te-Bi-Te-Pb-Te")
    #gen_struct_fix_block('2', 'trigonal', "Te-Bi-Te-Bi-Te,Te-Bi-Te-Pb-Te-Bi-Te,Te-Bi-Te-Pb-Te-Pb-Te-Bi-Te,Te-Pb-Te-Bi-Te-Bi-Te-Pb-Te")
    #gen_struct_fix_block('3', 'trigonal', "Te-Bi-Te-Bi-Te,Te-Bi-Te-Pb-Te-Bi-Te,Te-Bi-Te-Pb-Te-Pb-Te-Bi-Te")
    gen_struct_fix_block('3', 'trigonal', "Sb-Sb,Te-Sb-Te-Sb-Te,Te-Sb-Te-Ge-Te-Sb-Te,Te-Sb-Te-Ge-Te-Ge-Te-Sb-Te,Te-Ge-Te-Sb-Te-Sb-Te-Ge-Te")
    #gen_struct_fix_block('3', 'trigonal', "Sb-Sb,Te-Sb-Te-Sb-Te,Te-Sb-Te-Sb-Te-Sb-Te,Te-Sb-Te-Ge-Te-Ge-Te-Sb-Te,Te-Ge-Te-Sb-Te-Sb-Te-Ge-Te", {2.001: [1.552877096, 2.362777542, 2.739889686], 5.002: [1.712154812, 1.994616362, 1.994616362, 1.712154812, 2.632587528, 3.908657866], 7.003: [1.699051612, 2.035772228, 1.696968734, 1.696968734, 2.035772228, 1.699051612, 2.710683844, 4.224055311], 9.004: [1.695489891, 2.049167754, 1.669882087, 1.747324809, 1.747324809, 1.669882087, 2.049167754, 1.695489891, 3.13253256, 4.600418979], 9.005: [1.393991367, 2.199894149, 1.723797617, 2.01772035, 2.01772035, 1.723797617, 2.199894149, 1.393991367, 2.583808499, 4.264490915]})
