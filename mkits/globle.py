import numpy as np
from mkits.sysfc import lexit
from mkits.database import *
import ase.geometry
import spglib as spg
import os
import math

"""
:func effective_mass        : calculate effective mass based on parabolic band approximation, unit: m, J
:func best_polyfit_range    : find best fit range
:func globe_polyfit         : poly-fitting with numpy for r square
:func klist_kpath           : convert a list of kpoints to a kpath in units of reciprocal of angstrom
:func convert_3decimal_to_4 : convert 3 decimal numbers to 4 integer fraction: (0.625, 0.15625, 0.125) --> (4, 10, 8, 64)
:func np_ployfit            : 
:func round_even_odd        : round input varieties to nearest even number or odd number
:func vector_angle          : calculate the angle between two vectors
:func lattice_conversion    : convert lattice between cartesian vector and base vector
:func parser_inputpara      : get input parameters from string and return a dictionary, eg, oddkpoints:Ture;key2:attrib2 -> {"oddkpoints": "Ture"}
:func hstack_append_list    : append list2 to list1 in the horizontal direction despite of the equality of the length, [[1,2,3],[2,3,4]]+[[5,6,7],[6,7,8],[7,8,9]]->[[1, 2, 3, 5, 6, 7], [2, 3, 4, 6, 7, 8], [7, 8, 9]]
:func listcross             : get cross product of two list, [a, b, c]*[1,2,3]=[a, b, b, c, c, c]
:func frac2cart             :  
:func cart2frac             : 
:class struct               : 
"""


def effective_mass_calculator(kpath, eigen):
    """
    Calculate effective mass from single-band parabolic approximation
    :param kpath: a list or array like kpath (unit: angstrom^{-1})
    :param eigen: a list or array like eigen values of VBM or CBM band coorresponding to kpath dimension
    :return effective_mass, r square
    """
    if type(kpath) == list:
        kpath = np.array(kpath, dtype=float)
    if type(eigen) == list:
        eigen = np.array(eigen, dtype=float)
    
    coef1 = np.polyfit(kpath, eigen, 2)

    effective_mass = cons_hbar*cons_hbar/coef1[0]/2/cons_emass
    ssreg, sstot, rsquares = rsquare(kpath, eigen, 2)

    return effective_mass, rsquares


def best_polyfit_range(x, y, degree:int, center_index:int, init_range:int=5):
    """
    Search the best polyfitting range, by minimize the residuals
    :param center_index: the index of the chosen center, from 1
    :return array [[slice_left, slice_right, polynomial coefficients, r square],...[]]
    """
    
    leftrange_init = np.vstack((np.arange(center_index-init_range-1, -1, -1), np.arange(center_index-init_range-1, -1, -1))).T.flatten()[1:]
    rightrange_init = np.vstack((np.arange(center_index+init_range, len(y)), np.arange(center_index+init_range, len(y)))).T.flatten()

    # filling the range array
    rangediff = len(rightrange_init) - len(leftrange_init)
    if rangediff > 0:
        fittingrange = np.vstack(( np.hstack((leftrange_init, np.zeros(rangediff))), rightrange_init ))
    elif rangediff < 0:
        fittingrange = np.vstack(( leftrange_init, np.hstack(( rightrange_init, np.full(-rangediff, rightrange_init[-1]) )) ))
    
    fitparam = np.array([])
    for _ in fittingrange.T:
        fitparam = np.hstack((fitparam, np.array(_), rsquare(x[int(_[0]):int(_[1])], y[int(_[0]):int(_[1])], degree)))
    
    return fitparam.reshape(-1, degree+5)


def rsquare(x, y, degree):
    """
    :            (0     , 1:degree+2             , degree+3 or -1)
    :return ssregression, sstot, rsquare
    """
    results = [degree]

    coeffs = np.polyfit(x, y, degree)

    # r-squared
    p = np.poly1d(coeffs)
    yhat = p(x)
    ybar = np.sum(y)/len(y) 
    ssreg = np.sum((y - yhat)**2)
    sstot = np.sum((y - ybar)**2) 
    rsquare = 1 - ssreg/sstot
    

    return ssreg, sstot, rsquare


def klist_kpath(klist, recip):
    """
    convert a list of kpoints to a kpath in units of reciprocal of angstrom for band plotting or analyzing
    :param klist: array, x3 shape array
    :param recip: array, 3x3 shape array containing reciprocal lattice vector
    :return x shape array with kpath
    """
    kpath_abs = [0]
    cartesian_reciprocal = frac2cart(recip, klist)
    for _ in range(1, len(klist)):
        kpath_abs.append(kpath_abs[_-1]+np.sqrt(np.dot(cartesian_reciprocal[_]-cartesian_reciprocal[_-1], cartesian_reciprocal[_]-cartesian_reciprocal[_-1])))
    
    return np.array(kpath_abs)

def convert_3decimal_to_4integer_fractioin(x:float, y:float, z:float):
    """
    convert 3 decimal numbers to 4 integer fraction: (0.625, 0.15625, 0.125) --> (4, 10, 8, 64)
    maximum 10 digits 
    :return tuple
    """
    return (int(x*1e8), int(y*1e8), int(z*1e8), 1e8)


def np_ployfit(x, y, n, d):
    """ polyfit
    :param x: array or list
    :param y: array or list
    :param n: int, total points
    :param d: int, degree
    :return : dict, key = ["coefficient", "max_y", "min_y", "max_y_index", "min_y_index", "x", "y"] 
    """
    coe = np.polyfit(x, y, d)
    x = np.linspace(np.min(x), np.max(x), n)
    y = 0
    for dd in range(d+1):
        y += coe[dd]*x**(d-dd)
    
    return {
            "coefficient": coe, 
            "max_y": np.max(y),
            "min_y": np.min(y),
            "max_y_index": np.argmax(y),
            "min_y_index": np.argmin(y),
            "x": x,
            "y": y
            }


def round_even_odd(num, even_odd):
    """
    round input varieties to nearest even number or odd number
    :param num      : float or int, input 
    :param even_odd : int; 0 for even, 1 for odd, -1 for round only
    :route          :
                        num              3.5         2.5         2.9
                        num/2            1.75        1.25        1.45
                        math.modf(num/2) (0.75, 1)   (0.25, 1)   (0.45, 1)
    """
    if even_odd == 0:
        return int(num)+1 if math.modf(num/2)[0] >= 0.5 else int(num)
    elif even_odd ==1:
        num -= 1
        return int(num)+2 if math.modf(num/2)[0] >= 0.5 else int(num)+1
    elif even_odd == -1:
        return round(num)
    else:
        lexit("Error in even_odd, pls use [0 for even, 1 for odd, -1 for round only]")


def vector_angle(vector1, vector2, unit="deg"):
    """
    calculate the angle between two vectors
    :param vector1  : array, vector1
    :param vector2  : array, vector2
    :param unit     : string, optional [rad, deg]
    """
    unit_vector1 = vector1 / np.linalg.norm(vector1)
    unit_vector2 = vector2 / np.linalg.norm(vector2)
    dot_product = np.dot(unit_vector1, unit_vector2)
    angle = np.arccos(dot_product)
    return angle/np.pi*180 if unit == "deg" else angle


def lattice_conversion(give_lattice):
    """
    1. give cartesian lattice array; return vectors with angle:
    array([a_x, a_y, a_z], [b_x, b_y, b_z], [c_x, c_y, c_z]) -> array([a, b, c, alpha, beta, gamma])
    2. give vectors with angle; return cartisian lattice array:
    array[a, b, c, alpha, beta, gamma] -> array([a_x, a_y, a_z], [b_x, b_y, b_z], [c_x, c_y, c_z])
    """

    # cartesian -> vector
    if len(give_lattice) == 3:
        lattice_a = np.sqrt(np.dot(give_lattice[0, :], give_lattice[0, :]))
        lattice_b = np.sqrt(np.dot(give_lattice[1, :], give_lattice[1, :]))
        lattice_c = np.sqrt(np.dot(give_lattice[2, :], give_lattice[2, :]))
        angle_alpha = vector_angle(give_lattice[0, :], give_lattice[2, :])
        angle_beta  = vector_angle(give_lattice[2, :], give_lattice[1, :])
        angle_gamma = vector_angle(give_lattice[1, :], give_lattice[0, :])
        return(np.array([lattice_a, lattice_b, lattice_c, angle_alpha, angle_beta, angle_gamma]))

    # vector -> cartesian
    elif len(give_lattice) == 6:
        # convert arc to degree
        new_lattice = np.hstack((give_lattice[:3], give_lattice[3:]/180*np.pi))
        # a in x-axis; b in xy-plane
        vector_x = np.array([new_lattice[0], 0, 0])
        vector_y = np.array([np.cos(new_lattice[5]), np.sin(new_lattice[5]), 0])*new_lattice[1]
        x = np.cos(new_lattice[3])
        y = (np.linalg.norm(vector_y)*np.cos(new_lattice[4])-vector_y[0]*np.cos(new_lattice[3]))/vector_y[1]
        z = np.sqrt(1-x**2-y**2)
        vector_z_direction = np.array([x,y,z])
        vector_z = vector_z_direction/ np.linalg.norm(vector_z_direction) * new_lattice[2]
        return np.array([vector_x, vector_y, vector_z])
        
    else:
        lexit("Error, give the lattice with following format: array([a_x, a_y, a_z], [b_x, b_y, b_z], [c_x, c_y, c_z]) or array([a, b, c, alpha, beta, gamma])")


def parser_inputpara(inputstring):
    """
    :param inputstring: the input parameters are seperated by ",", and key and attribution are seperated by ":"
    :return : dictionary
    """
    input_dict = {}
    if "," in inputstring: separator_outkey = ","
    elif ";" in inputstring: separator_outkey = ";"
    else: separator_outkey = " "
    if ":" in inputstring: separator_inkey = ":"
    elif "=" in inputstring: separator_inkey = "="
    try:
        for inp in inputstring.split(separator_outkey):
            inp_key_para = inp.split(separator_inkey)
            input_dict[inp_key_para[0]] = inp_key_para[1]
    except:
        lexit("Error: make sure the string like input parameters with following format: key1:para1,key2:para2")
    return input_dict


def hstack_append_list(list1, list2):
    """
    append list2 to list1 in the horizontal direction: 
    list1
    [["1","2"],        [["8", "9"]             [["1","2","8","9"],
     ["3","4"],    +    ["10","11"]]    ->      ["3","4","10","11"],
     ["5","7"]]                                 ["5","7"]]
    :param list1 and list2
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
        print("Input is not a list")
    return list1


def listcross(list1, list2):
    """
    list2 must be a numeric type list and len(list1) = len(list2)
    example:
    listcross(['Ge', 'Te'], np.array([2.0, 5.0]))
    listcross(['Ge', 'Te'], [2.0, 5.0])
    output: ['Ge', 'Ge', 'Te', 'Te', 'Te', 'Te', 'Te']
    """
    if (type(list1) != list or type(list2) != list) and type(list2) != np.ndarray:
        lexit("The input lists are not list type.")
    if (type(list1) != list or type(list2) != list) and type(list2) != np.ndarray:
        lexit("The input lists are not list type.")
    try:
        num_list = [int(i) for i in list2]
    except:
        lexit("Wrong numeric list")
    crosslist = []
    for i in range(len(num_list)):
        crosslist += [list1[i]]*num_list[i]
    return crosslist


def frac2cart(cart_lattice, fraction_pos):
    """
    fractional coordinates to cartisian coordinates
    cartesian lattice: [[4.45, 0, 0], [-2.22, 3.86, 0], [0, 0, 17.50]]
    fractional coordinates: [[], []]
    :para 
    """
    cart_ = np.array([0,0,0])
    for _ in range(len(fraction_pos)):
        cart_ = np.vstack((cart_, np.sum(cart_lattice*fraction_pos[_], axis=0)))
    return cart_[1:]


def cart2frac(cart_lattice, cart_pos):
    """
    fractional coordinates to cartisian coordinates
    """

    return np.matmul(cart_pos, np.matrix(cart_lattice.T).I)


class struct:
    """
    structure class
    :input struct_file      : strings, case.vasp, POSCAR
    :func return_dict       : return a dictionary of the structure
    :varieties struct_dict  : key = [
                                        "atom_type": a list contain the alphabetic name of atoms
                                        "atom_num": an array contain the number of each atom in atom_type
                                    ]
    """

    def __init__(self, struct_file): 
        """read from dictionary or string or a file name"""
        try:
            if type(struct_file) == dict:
                self.struct_dict = struct_file
                self.calculator = struct_file["calculator"]
            elif type(struct_file) == str and "\n" in struct_file:
                lexit("The string func hasn't  been added yet.")
            else:
                with open(struct_file, "r") as f: 
                    struct_lines = f.readlines() # a list of each line of the structure file
                self.struct_file = struct_file
                self.calculator = self.parse_calculator(struct_lines) # vasp, wien, qe
                self.struct_dict = self.parse_struct(struct_lines, self.calculator)
        except:
            lexit("Cannot find the structure file: ", struct_file)

    def __repr__(self) -> str:
        """"""
        return "This is a Class containing structures information: \n" + "Lattice:\n{}\n".format(str(self.struct_dict["lattice"]))
    
    def return_dict(self):
        return self.struct_dict
    
    def return_recip(self):
        return ase.geometry.Cell.fromcellpar([self.struct_dict["lattice_direct"][0], self.struct_dict["lattice_direct"][1], self.struct_dict["lattice_direct"][2], self.struct_dict["lattice_direct"][3], self.struct_dict["lattice_direct"][4], self.struct_dict["lattice_direct"][5]]).reciprocal()*2*np.pi

    def parse_calculator(self, struct_lines):
        """ 
        return calculator
        :return vasp_dyn:  
        :return vasp_direct:
        :return vasp_cart:
        :return wien:
        :return qeout:
        :return qein:
        """
        # vasp wien qeou
        if "S" in struct_lines[7] and "Select" in struct_lines[7]:
            print("Read a selective dynamic VASP file: "+self.struct_file)
            return("vasp_dyn")
        elif "D" in struct_lines[7] or "Cartesian" in struct_lines[7] or "Cartesian" in struct_lines[8]:
            print("Read a VASP file: "+self.struct_file)
            return "vasp_direct"
        elif "unit=bohr" in struct_lines[2] or "RELA" in struct_lines[2]:
            print("Read a WIEN2k file: "+self.struct_file)
            return "wien"
        elif "Quantum ESPRESSO" in "".join(struct_lines[:]): 
            print("Read a QE output file: "+self.struct_file)
            return "qeout"
        elif "&CONTROL" in "".join(struct_lines[:]):
            print("Read a QE input file: "+self.struct_file)
            return "qein"

    def parse_struct(self, struct_lines, calculator): # fraction position, angstrom
        struct_dict = {
            "title"      : "",
            "atoms_type" : [],
            "atoms_num"  : [],
            "volume"     : 0,
            "calculator" : ""
        }
        if calculator == "vasp_direct":
            struct_dict["title"] = struct_lines[0][:-1]
            struct_dict["ratio"] = float(struct_lines[1][:-1])
            struct_dict["lattice"] = np.array([list(map(float, struct_lines[2].split())), list(map(float, struct_lines[3].split())), list(map(float, struct_lines[4].split()))])
            struct_dict["lattice_direct"] = lattice_conversion(struct_dict["lattice"])
            struct_dict["atoms_type"] = struct_lines[5].split()
            struct_dict["atoms_num"] = np.array(list(map(int, struct_lines[6].split())))
            struct_dict["atoms_index"] = [symbol_map[_] for _ in struct_dict["atoms_type"]]
            struct_dict["pos_frac"] = np.array([list(map(float, (i[:-1].split())[:3])) for i in struct_lines[8:8+np.sum(struct_dict["atoms_num"])]])
            struct_dict["calculator"] = "vasp_direct"
        elif calculator == "vasp_dyn":
            struct_dict["title"] = struct_lines[0][:-1]
            struct_dict["ratio"] = float(struct_lines[1][:-1])
            struct_dict["lattice"] = np.array([list(map(float, struct_lines[2].split())), list(map(float, struct_lines[3].split())), list(map(float, struct_lines[4].split()))])
            struct_dict["lattice_direct"] = lattice_conversion(struct_dict["lattice"])
            struct_dict["atoms_type"] = struct_lines[5].split()
            struct_dict["atoms_num"] = np.array(list(map(int, struct_lines[6].split())))
            struct_dict["atoms_index"] = [symbol_map[_] for _ in struct_dict["atoms_type"]]
            struct_dict["pos_frac"] = np.array([list(map(float, (i[:-1].split())[:3])) for i in struct_lines[9:9+np.sum(struct_dict["atoms_num"])]])
            struct_dict["pos_move"] = [list(map(str, (i.split())[3:])) for i in struct_lines[9:9+np.sum(struct_dict["atoms_num"])]]
            struct_dict["calculator"] = "vasp_dyn"
        elif calculator == "qeout":
            # get data line: celldm(1); crystal axes; site n.     atom                  positions (alat units); number of k points=    12
            dataget_num = 3
            lattice_inx = 0
            atom_position_bgn_inx = 0
            atom_num_inx = 0
            for i in range(1000):
                if "CELL_PARAMETERS" in struct_lines[-i]:
                    lattice_inx = -i
                    dataget_num -= 1
                elif "ATOMIC_POSITIONS" in struct_lines[-i]:
                    atom_position_bgn_inx = -i
                    dataget_num -= 1
                elif "number of atoms/cell" in struct_lines[-i]:
                    atom_num_inx = -i
                    dataget_num -= 1
                else:
                    continue

                if dataget_num > 0 and i < 999:
                    continue
                elif dataget_num == 0 and i < 999:
                    break
                else:
                    lexit("Cannot find expected information in QE structure file")
            #print(struct_lines[celldm_inx], struct_lines[lattice_inx], struct_lines[atom_cartesian_bgn_inx], struct_lines[atom_cartesian_end_inx])
            total_atoms = int(struct_lines[atom_num_inx].split()[-1])
            struct_dict["lattice"] = np.array([list(map(float, i.split()[-3:])) for i in struct_lines[lattice_inx+1:lattice_inx+4]])
            atoms_pos = ([i.split()[:1][0] for i in struct_lines[atom_position_bgn_inx+1:atom_position_bgn_inx+1+total_atoms]])[:]
            struct_dict["pos_frac"] = np.array([list(map(float, i.split()[-3:])) for i in struct_lines[atom_position_bgn_inx+1:atom_position_bgn_inx+1+total_atoms]])[np.argsort(atoms_pos)]
            atoms_pos = np.array(atoms_pos)[np.argsort(atoms_pos)]
            struct_dict["atoms_type"] = list(dict.fromkeys(atoms_pos))
            struct_dict["atoms_num"] = [list(atoms_pos).count(i) for i in struct_dict["atoms_type"]]
            struct_dict["atoms_sequence"] = listcross(struct_dict["atoms_type"], struct_dict["atoms_num"])
        elif calculator == "wien":
            struct_dict["title"] = struct_lines[0][:-1]
            #struct_dict["sym_type"] = struct_lines[1][:1]
            struct_dict["sym_type"] = "P"
            struct_dict["lattice_direct"] = np.array([float(struct_lines[3][0:10])*uc_bohr2ang, float(struct_lines[3][10:20])*uc_bohr2ang, float(struct_lines[3][20:30])*uc_bohr2ang, float(struct_lines[3][30:40]), float(struct_lines[3][40:50]), float(struct_lines[3][50:60])])
            struct_dict["lattice"] = lattice_conversion(struct_dict["lattice_direct"])
            struct_dict["atoms_alias"] = []
            struct_dict["atoms_rmt"] = []
            struct_dict["atoms_index"] = []
            struct_dict["atoms_r0"] = []
            struct_dict["calculator"] = "wien"

            # the file with "ATOM"
            line_idx = 4
            struct_dict["pos_frac"] = np.array([[0,0,0]])
            while not "OPERATIONS" in struct_lines[line_idx]:
                # line line_idx
                struct_dict["pos_frac"] = np.vstack((struct_dict["pos_frac"], np.array([float(struct_lines[line_idx][12:22]), float(struct_lines[line_idx][25:35]), float(struct_lines[line_idx][38:48])])))
                # line_idx +1
                atom_mlti = int(struct_lines[line_idx+1][15:17])
                for _ in range(1, atom_mlti):
                    struct_dict["pos_frac"] = np.vstack((struct_dict["pos_frac"], np.array([float(struct_lines[line_idx+_+1][12:22]), float(struct_lines[line_idx+_+1][25:35]), float(struct_lines[line_idx+_+1][38:48])])))
                struct_dict["atoms_alias"] += [(struct_lines[line_idx+atom_mlti+1][:5]).replace(" ", "")] * atom_mlti
                struct_dict["atoms_r0"] += [float(struct_lines[line_idx+atom_mlti+1][25:34])] * atom_mlti
                struct_dict["atoms_rmt"] += [float(struct_lines[line_idx+atom_mlti+1][40:48])] * atom_mlti
                struct_dict["atoms_index"] += [int(float(struct_lines[line_idx+atom_mlti+1][55:-1]))] * atom_mlti
                line_idx += atom_mlti+5
                if line_idx > len(struct_lines):
                    lexit("Cannot read a correct wien2k structure, make sure you have <SYMMETRY OPERATIONS> at the end of structure file")
            struct_dict["pos_frac"] = struct_dict["pos_frac"][1:, :]
            struct_dict["atoms_type"] = [atom_data[_][1] for _ in struct_dict["atoms_index"]]
            struct_dict["atoms_tot"] = len(struct_dict["atoms_type"])
            
        return struct_dict
    
    def get_refined_struct(self, symprecision=1e-2, primitive=False):
        self.struct_dict["atoms_sequence"] = listcross(self.struct_dict["atoms_type"], self.struct_dict["atoms_num"])
        spg_number = [symbol_map[i] for i in self.struct_dict["atoms_sequence"]]
        spg_tuple = (np.ndarray.tolist(self.struct_dict["lattice"]), np.ndarray.tolist(self.struct_dict["pos_frac"]), spg_number)
        lattice, positions, spgnumber = spg.refine_cell(spg_tuple, symprec=symprecision)

        spg_indx = np.argsort(np.array(spgnumber))
        spgnumber = list(np.array(spgnumber)[spg_indx])
        atom_type = list(set(spgnumber))
        atom_num = np.array([spgnumber.count(i) for i in atom_type])
        
        self.struct_dict["lattice"] = np.array(lattice)
        self.struct_dict["pos_frac"] = np.array(positions)[spg_indx, :]
        self.struct_dict["atoms_type"] = [atom_data[i][1] for i in atom_type]
        self.struct_dict["atoms_num"] = atom_num
    
    def update_struct_dict(self, new_dict):
        """ update struct dictionary with new and changed structure """
        self.struct_dict = new_dict
    
    def write_struct(self, fpath="./", fname="POSCAR"):
        if self.calculator == "vasp_direct":
            with open(fpath+"/"+fname, "w") as f:
                f.write(self.struct_dict["title"]+" generate by mkits\n")
                f.write("{:<.5f}".format(self.struct_dict["ratio"])+"\n")
                np.savetxt(f, self.struct_dict["lattice"], fmt="%20.10f")
                for i in self.struct_dict["atoms_type"]: f.write("{:>5s}".format(i))
                f.write("\n")
                np.savetxt(f, self.struct_dict["atoms_num"].reshape((1,-1)), fmt="%5d")
                f.write("Direct\n")
                np.savetxt(f, self.struct_dict["pos_frac"][:int(sum(self.struct_dict["atoms_num"]))], fmt="%20.16f")
        elif self.calculator == "vasp_dyn":
            with open(fpath+"/"+fname, "w") as f:
                f.write(self.struct_dict["title"]+" generate by mkits\n")
                f.write("{:<.5f}".format(self.struct_dict["ratio"])+"\n")
                np.savetxt(f, self.struct_dict["lattice"]*self.struct_dict["ratio"], fmt="%20.10f")
                for i in self.struct_dict["atoms_type"]: f.write("{:>5s}".format(i))
                f.write("\n")
                np.savetxt(f, self.struct_dict["atoms_num"].reshape((1,-1)), fmt="%5d")
                f.write("Selective dynamics\nDirect\n")
                np.savetxt(f, np.hstack((self.struct_dict["pos_frac"], self.struct_dict["pos_move"])) , fmt="%20s %20s %20s %6s %6s %6s")
        elif self.calculator == "wien":
            with open(fpath+"/"+fname, "w") as f:
                f.write("generate by mkits\n")
                f.write("%-4sLATTICE,NONEQUIV.ATOMS%4d%11s\n" % ("H", self.struct_dict["atoms_tot"], "P1"))
                f.write("MODE OF CALC=RELA unit=bohr                                                    \n")
                f.write("%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n" % (self.struct_dict["lattice_direct"][0]/0.529177, self.struct_dict["lattice_direct"][1]/0.529177, self.struct_dict["lattice_direct"][2]/0.529177, self.struct_dict["lattice_direct"][3], self.struct_dict["lattice_direct"][4], self.struct_dict["lattice_direct"][5]))
                for _ in range(self.struct_dict["atoms_tot"]):
                    f.write("ATOM%4d: X=%-10.8f Y=%-10.8f Z=%-10.8f\n" % (_+1, self.struct_dict["pos_frac"][_,0], self.struct_dict["pos_frac"][_,1], self.struct_dict["pos_frac"][_,2]))
                    f.write("          MULT= 1          ISPLIT= 4                                           \n")
                    f.write("%-5s      NPT=  781  R0=%10.8f RMT=%8.5f     Z:%10.5f\n" % (self.struct_dict["atoms_alias"][_], self.struct_dict["atoms_r0"][_], self.struct_dict["atoms_rmt"][_], self.struct_dict["atoms_index"][_]))
                    f.write("LOCAL ROT MATRIX:    1.0000000 0.0000000 0.0000000                             \n")
                    f.write("                     0.0000000 1.0000000 0.0000000                             \n")
                    f.write("                     0.0000000 0.0000000 1.0000000                             \n")
                f.write("   1      NUMBER OF SYMMETRY OPERATIONS                                        \n 1 0 0 0.00000000                                                              \n 0 1 0 0.00000000                                                              \n 0 0 1 0.00000000                                                              \n       1\n")
                
    def gen_potcar(self, potpath: str = "./", fpath: str ="./", fname: str ="POTCAR"):
        # VASP] generate POTCAR
        atoms = self.struct_dict["atoms_type"]
        potcar = []
        for atom in atoms:
            try:
                potcar += open(potpath+"/"+atom+"/POTCAR", "r").readlines()
            except:
                try:
                    potcar += open(potpath+"/POTCAR_"+atom, "r").readlines()
                except:
                    lexit("The POTCAR of " + atom + " doesn't exit.")
        if os.path.exists(fpath+"/POTCAR"):
            lexit("POTCAT exits, rename it and re-excute the code.")
        else:
            with open(fpath+"/"+fname, "w") as f:
                f.writelines(potcar)