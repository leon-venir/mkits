# -*- coding: utf-8 -*

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mkits.sysfc import lexit
from mkits.database import *
import ase.geometry
import ase.io
import spglib as spg
import os
import math
import sys


"""
Class
-----
struct:
    A class for read, write, operate the material structure.
struct_ase: 
    A new class of structures based on ase library.
about:
    A class for parsering the dft output files.

Functions
---------
list_index_by_list: 
    index a list with another list, the way in numpy
trans_reflection_xy    :
split_list             : 
arith_prog_split_even  : split the arithmetic progression into approximate 
                         equal part.
coord_trans_xy         : coordinates translation in x-y plane
coord_rotate_xy        : oordinates rotation in x-y plane
parser_inputfile       : get input parameters from inputfile
gnudata2numpy          : convert the data blocks of gnuplot to numpy readable text.
progressbar            : draw processing bar
center_array           : generate an array with specific center value,step 
                         and total number
del_list_dupli_neighbor: removing neighboring duplicates in list
effective_mass         : calculate effective mass based on parabolic band
                         approximation, unit: m, J
best_polyfit_range     : find best fit range
globe_polyfit          : poly-fitting with numpy for r square
klist_kpath            : convert a list of kpoints to a kpath in units of 
                         reciprocal of angstrom
convert_3decimal_to_4  : convert 3 decimal numbers to 4 integer fraction: 
                         (0.625, 0.15625, 0.125) --> (4, 10, 8, 64)
np_ployfit             : 
round_even_odd         : round input varieties to nearest even number or 
                         odd number
vector_angle           : calculate the angle between two vectors
lattice_conversion     : convert lattice between cartesian vector and base vector
parser_inputpara       : get input parameters from string and return a dictionary, 
                         eg, oddkpoints:Ture;key2:attrib2 -> {"oddkpoints": "Ture"}
parser_inputlines      : get input parameters from file and return a dictionary
hstack_append_list     : append list2 to list1 in the horizontal direction 
                         despite of the equality of the length, 
                         [[1,2,3],[2,3,4]]+[[5,6,7],[6,7,8],[7,8,9]]->
                         [1, 2, 3, 5, 6, 7], [2, 3, 4, 6, 7, 8], [7, 8, 9]]
listcross              : get cross product of two list, 
                         [a, b, c]*[1,2,3]=[a, b, b, c, c, c]
frac2cart              :  
cart2frac              : 
"""


def list_index_by_list(list1, index):
    """
    T = [L[i] for i in Idx]
    """
    return [list1[i] for i in index]


def trans_reflection_xy(vector, method, inp):
    """
    Parameters
    ----------
    vector: array
        The original vector
    method: str [2ps]
        The form of reflection line
        2ps: 2 points
    inp: tuple
        2ps: ([1, 0], [2,1])
    Returns
    -------
        A reflected vector
    """
    reflected = np.array([0])
    if method == "2ps":
        point1 = np.array(inp[0])
        point2 = np.array(inp[1])

        p1p2xdiff = point1[0]-point2[0]
        p1p2ydiff = point1[1]-point2[1]
        if p1p2xdiff < 1e-8:
            reflected = vector + np.array([-p1p2xdiff, 0])
            reflected *= np.array([1,-1])
            reflected += np.array([p1p2xdiff, 0])


def split_list(a, n):
    """

    """
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))


def arith_prog_split_even(n):
    """ 
    Algorithm
    ---------
    20, 19, 18, ..., 3, 2, 1
        |

    \left\{\begin{matrix}
        x_{10} = n \\ 
        \frac{x_{1} + 0}{2} \frac{x_{1}}{d} \approx 
        \frac{x_{2} + x_{1}}{2} \frac{x_{2} - x_{1}}{d} \approx 
        \frac{x_{3} + x_{2}}{2} \frac{x_{3} - x_{2}}{d} \approx 
        \frac{x_{4} + x_{3}}{2} \frac{x_{4} - x_{3}}{d} 
    \end{matrix}\right.

    Parameters
    """


def vector_trans_xy(vector, mv_vector):
    """
    Parameters
    ----------
    vector: array
        The original vector
    mv_vector: array
        The translation vector
    Returns
    -------
        A translated vector
    """
    return vector + mv_vector


def vector_rotate_xy(vector, angle, clockwise=True):
    """
    Parameters
    ----------
    """
    vx = np.array([])
    if clockwise:
        vx = np.array([[np.cos(angle), np.sin(angle), 0],
                        [-np.sin(angle), np.cos(angle), 0],
                        [0,0,1]])
    else:
        vx = np.array([[np.cos(angle), -np.sin(angle), 0],
                        [np.sin(angle), np.cos(angle), 0],
                        [0,0,1]])
    return (vx @ vector.T).T


def parser_inputfile(inp, comment_sym="#", assign_sym="=", block=False):
    """
    Parameters
    ----------
    inp:str
        the name of input file
    style: str
        option [keywords, block] 
    comment_sym: str
        comment symbol, default #
    assign_sym: str
        value assigned symbol
    
    Returns
    -------
    dictionary
    """
    keywords = {}

    with open(inp, "r") as f:
        lines = f.readlines()

    if block:
        pass
    else:
        blockoutlist = list(range(len(lines)))
        
    # delete comments and \n
    for i in blockoutlist:
        lines[i] = lines[i].strip()
        try:
            harshsym_index = lines[i].index(comment_sym)
            lines[i] = lines[i][:harshsym_index]
        except:
            continue
    
    # delete blank line and find the keys - values
    for i in blockoutlist:
        if lines[i] == "\n":
            blockoutlist.remove(i)
        else:
            try:
                assign_sym_index = lines[i].index(assign_sym)
                key = lines[i][:assign_sym_index].strip()
                value = lines[i][assign_sym_index+1:].strip()
                keywords[key] = value
            except:
                continue
    
    return keywords


def gnudata_split2numpy(gnufile:str, breakline:int=0):
    """Return an array read from gnuplot data 

    Parameters
    ----------
    gnufile: str
        File name
    breakline: int
        line index from 1

    Returns
    -------
    array with 
    """
    gnudata = np.loadtxt(gnufile)

    #f len(gnudata)%breakline != 0:
    #    lexit("Wrong size of the each data block.")
    
    gnudata_block_size = int(len(gnudata[:,0])/breakline)
    gnudata_clumn_init = len(gnudata[0, :])

    # reshape the array
    gnudata = gnudata.reshape(gnudata_block_size, breakline, -1)
    gnudata = gnudata.transpose(1,0,2).reshape(breakline,-1)
    gnudata = np.delete(gnudata, 
                        np.s_[gnudata_clumn_init::gnudata_clumn_init],
                        axis=1)

    return gnudata


def progressbar(total_step:int, current_step:int):
    """
    draw progressing bar
    """
    # percent
    i = int(current_step/total_step*50)

    sys.stdout.write('\r')
    sys.stdout.write("[%-50s] %d%%" % ('='*i, 2*i))
    sys.stdout.flush()

    #for i in range(21):
    #    sys.stdout.write('\r')
    #    sys.stdout.write("[%-100s] %d%%" % ('='*i, 1*i))
    #    sys.stdout.flush()


def center_array(center:float=1.0, step:float=0.005, total_num:int=5):
    """
    """
    if total_num%2 == 1:
        return np.linspace(center-step*(total_num-1)/2, 
                           center+step*(total_num-1)/2, 
                           total_num, endpoint=True)
    elif total_num % 2 == 0:
        return np.linspace(center-step*total_num/2, 
                           center+step*(total_num-1)/2, 
                           total_num, endpoint=True)


def del_list_dupli_neighbor(inp_list:list):
    """ Removing neighboring duplicates in list """
    res = [inp_list[0]]
    for i, c in enumerate(inp_list[1:]):
        if c != inp_list[i]:
            res.append(c)
    return res


def carrier_mobility(cii, ed, T, dimension:str="3d"):
    """
    Calculate relaxation time based on Deformation Potential theory 
    Parameters
    ----------
    cii: float
        elastic constants in Pa
    effect_mass: float
        effective mass
    ed: float 
        deformation potenital
    T: float
        temperature in K
    dimension: str
        option [3d, 2d]
    """
    if dimension == "3d":
        mu = 2 * np.sqrt(2*np.pi) * cons_echarge * cons_hbar**4 * cii
        mu /= 3 * cons_emass**(5/2) * (cons_kb*T)**(3/2) * ed**2
    return mu


def relaxation_time(cii, effect_mass, ed, T, dimension:str="3d"):
    """
    Calculate relaxation time based on Deformation Potential theory 
    Parameters
    ----------
    cii: float
        elastic constants in Pa
    effect_mass: float 
        effective mass 
    ed: float
        deformation potenital in eV
    T: float
        temperature in K
    dimension: str
        optional [3d, 2d]
    """
    if dimension == "3d":
        #tau = 2 * np.sqrt(2*np.pi) * cons_hbar**4 * cii * effect_mass 
        #tau /= 3 * cons_emass**(3/2) * (cons_kb*T)**(3/2) * ed**2
        tau = 2 * np.sqrt(2*np.pi) * cons_hbar**4 * cii * uc_gpa2nm2
        tau /= 3 * (cons_emass*effect_mass)**(3/2) * (cons_kb*T)**(3/2) * (ed * uc_ev2j)**2
    elif dimension == "2d":
        mu = cons_echarge * cons_hbar**3 * cii
        mu /= cons_kb * T * (effect_mass * cons_emass)**2 * (ed * uc_ev2j)**2
        tau = mu * (effect_mass * cons_emass) / cons_echarge
    return tau


def mean_free_path(cii, ed, T, dimension:str="3d"):
    """
    Calculate relaxation time based on Deformation Potential theory 
    :param cii: elastic constants in Pa
    :param effect_mass: effective mass 
    :param ed: deformation potenital
    :param T: temperature in K
    :param dimension: [3d, 2d]
    """
    if dimension == "3d":
        mean_free_path = np.pi * cons_hbar**4 * cii / \
                         ( cons_emass**2 * ed**2 * cons_kb * T)
    return mean_free_path


def effective_mass_calculator(kpath, 
                              eigen, 
                              fpath:str="./", 
                              fname="fitting.png", 
                              plot:bool=True):
    """
    Calculate effective mass from single-band parabolic approximation

    Parameters:
    -----------
    kpath: a list or array like kpath (unit: angstrom^{-1})
    eigen: a list or array like eigen values of VBM or CBM band coorresponding 
           to kpath dimension
    Return:
    -------
    effective_mass, r square
    """
    if type(kpath) == list:
        kpath = np.array(kpath, dtype=float)
    if type(eigen) == list:
        eigen = np.array(eigen, dtype=float)
    
    coef1 = np.polyfit(kpath, eigen, 2)

    effective_mass = cons_hbar*cons_hbar/coef1[0]/2/cons_emass
    ssreg, sstot, rsquares = rsquare(kpath, eigen, 2)

    if plot:
        plt.plot(kpath, eigen, marker="o")
        plt.xlabel("Reciprocal space (agstrom$^{-1}$)")
        plt.ylabel("Energy (J)")
        plt.title("effective mass: %10.5f; rsquare: %15.9f" % (effective_mass, rsquares))
        plt.savefig("%s/%s" %(fpath, fname))
        plt.close()

    return effective_mass, rsquares


def best_polyfit_range(x, y, degree:int, min_points:int=10, max_points:int=20):
    """
    Search the best polyfitting range, by minimize the residuals
    :param x:
    :param y:
    :return array [[best_beg_index, best_end_index, best_rsqure], [all fitting data]]
    """
    if len(x) != len(y):
        lexit("Input x and y must be 1-dimensional array and process same length,")
    if max_points < min_points:
        lexit("max_points mush be larger than min_points.")
    if max_points > len(x):
        lexit("The number of max_points should be smaller than the length of x.")
    ssreg, sstot, rsquares = rsquare(x, y, degree=degree)
    fit_history = np.array([0, len(x)-1, rsquares])

    # to 
    step = 0
    total_step = int((2*len(y) - min_points - max_points)/2*(max_points-min_points+1))
    print("Searching the best fitting range: ")
    for points in range(min_points, max_points+1):
        beg_index = 0
        while beg_index + points <= len(y):
            ssreg, sstot, rsquares = rsquare(x[beg_index:beg_index+points], 
                                             y[beg_index:beg_index+points], 
                                             degree=degree)
            fit_history = np.vstack((fit_history, 
                                     np.array([beg_index, beg_index+points, rsquares])))
            beg_index += 1
            # progressing bar
            progressbar(total_step=total_step, current_step=step)
            step += 1
    max_rsqures_index = np.argmax(fit_history[:, 2])

    return np.vstack((fit_history[max_rsqures_index], fit_history))


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
    Convert a list of kpoints to a kpath in units of reciprocal of angstrom 
    for band plotting or analyzing

    Parameters:
    -----------
    klist: array, x3 shape array
    recip: array, 3x3 shape array containing reciprocal lattice vector
    
    Return:
    -------
    x shape array with kpath
    """
    kpath_abs = [0]
    cartesian_reciprocal = frac2cart(recip, klist)
    for _ in range(1, len(klist)):
        kpath_abs.append(kpath_abs[_-1]+np.sqrt(np.dot(cartesian_reciprocal[_]-cartesian_reciprocal[_-1], 
                                                       cartesian_reciprocal[_]-cartesian_reciprocal[_-1])))
    
    return np.array(kpath_abs)

def convert_3decimal_to_4integer_fractioin(x:float, y:float, z:float):
    """
    convert 3 decimal numbers to 4 integer fraction: (0.625, 0.15625, 0.125) --> (4, 10, 8, 64)
    maximum 10 digits 
    :return tuple
    """
    return (int(x*1e8), int(y*1e8), int(z*1e8), 1e8)


def np_ployfit(x, y, n, d):
    """
    Polyfit
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
    Parameters
    ----------
    vector1: array
        The first vector
    vector2: array
        The second vector
    unit:str
        units, option [rad, deg]
    Returns
    -------
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



def parser_inputlines(input_lines, comment_letter:str="#"):
    """
    :param input_lines: list or filename, input lines from readlines()
    :return : dictionary
    """
    if isinstance(input_lines, list):
        pass
    else:
        with open(input_lines) as f:
            input_lines = f.readlines()

    input_dict = {}
    for line in input_lines:
        line = line.replace(" ", "")
        line = line.replace(",", "")
        if comment_letter in line:
            line = line[:line.index(comment_letter)]
        else:
            line = line.replace("\n", "")
        if "=" in line:
            line = line.split("=")
            input_dict[line[0]] = str(line[1])
        elif ":" in line:
            line = line.split(":")
            input_dict[line[0]] = str(line[1])
            
    return input_dict


def parser_inputpara(inputstring):
    """
    :param inputstring: the input parameters are seperated by ",", and key and attribution are seperated by ":"
    :return : dictionary
    """
    input_dict = {}

    # get the separator
    if "," in inputstring: 
        separator_outkey = ","
    elif ";" in inputstring: 
        separator_outkey = ";"
    else: separator_outkey = " "

    # delete the redundant separator
    if inputstring[0] == separator_outkey:
        inputstring = inputstring[1:]
    if inputstring[-1] == separator_outkey:
        inputstring = inputstring[:-1]

    # get the assignment symbol
    if ":" in inputstring: 
        separator_inkey = ":"
    elif "=" in inputstring: 
        separator_inkey = "="
        
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
    cart_ = np.zeros(3)

    if fraction_pos.shape == (1,3):
        fraction_pos = np.vstack((fraction_pos,
                                  np.array([[0,0,0]])))
        for _ in range(len(fraction_pos)):
            cart_ = np.vstack((cart_, 
                               np.sum(cart_lattice*fraction_pos[_],
                                      axis=0)))
        cart_ = cart_[:-1]
    else:
        for _ in range(len(fraction_pos)):
            cart_ = np.vstack((cart_, np.sum(cart_lattice, axis=0)))
    
    return cart_[1:]


def cart2frac_single(lattice, cart_pos):
    """
    fractional coordinates to cartesian coordinates
    """
    lattice_a = 0
    lattice_b = 0
    lattice_c = 0
    alpha = 0
    beta = 0
    gamma = 0

    if len(lattice) == 3:
        cart_lattice = lattice
        cart_lattice_angle = lattice_conversion(give_lattice=lattice)
        lattice_a = cart_lattice_angle[0]
        lattice_b = cart_lattice_angle[1]
        lattice_c = cart_lattice_angle[2]
        alpha = cart_lattice_angle[3]
        beta = cart_lattice_angle[4]
        gamma = cart_lattice_angle[5]
    elif len(lattice) == 6:
        cart_lattice = lattice_conversion(give_lattice=lattice)
        cart_lattice_angle = lattice
        lattice_a = cart_lattice_angle[0]
        lattice_b = cart_lattice_angle[1]
        lattice_c = cart_lattice_angle[2]
        alpha = cart_lattice_angle[3]
        beta = cart_lattice_angle[4]
        gamma = cart_lattice_angle[5]
    
    # for o
    frac_a, frac_b, frac_c = 0, 0, 0
    if abs(cart_lattice_angle[3]-90) < 0.1 and abs(cart_lattice_angle[4]-90) < 0.1:
        frac_a = (cart_pos[0] + cart_pos[1]*np.tan((gamma-90)*uc_d2a))/lattice_a
        frac_b = (cart_pos[1]/np.cos((gamma-90)*uc_d2a))/lattice_b
        frac_c = cart_pos[2]/lattice_c

    return frac_a, frac_b, frac_c


def cart2frac(cart_lattice, cart_pos):
    """
    fractional coordinates to cartisian coordinates
    """

    #return np.matmul(cart_pos, np.matrix(cart_lattice.T).I)
    frac_pos = np.zeros((len(cart_pos), 3))
    for i in range(len(cart_pos)):
        frac_a, frac_b, frac_c = cart2frac_single(cart_lattice, cart_pos[i])
        frac_pos[i, 0] = frac_a
        frac_pos[i, 1] = frac_b
        frac_pos[i, 2] = frac_c
    
    return frac_pos


class struct(object):
    """
    structure class
    Parameters
    ----------
    inp, string
        The name of the input structure 

    Attributs
    ---------
    atom_continous_num
        The continous indexing of atom number
        atom number [60 18 12  4] --> [60, 78, 90, 94]
    nearest_neighbor, array
        The list of 5 nearest neighbor and their distance
        [[atom1, neighbor1, distance_to_neighbor1, neighbor2, dis...],
         [atom2, neighbor1, distance_to_neighbor1, neighbor2, dis...],
         ...]]
    struct_dict
        A dictionary
    
    Functions
    ---------
    return_dict
        return a dictionary of the structure
    del_atom
        delete atoms by elements name, by index, by coordinates
    """

    def __init__(self, inp): 
        """read from dictionary or string or a file name"""
        self.total_atom_num = 0
        self.title = "none"
        self.ratio = 1.0
        # cartesian lattice xx xy xz yz yy yz zx zy zz
        self.lattice9 = np.array([])
        # direct lattice a b c alpha beta gamma
        self.lattice6 = np.array([])
        self.coord_direct = np.array([])
        self.coord_cart = np.array([])
        self.atom_dyn = []
        self.dyn = False
        self.atom_num = np.array([])
        self.atom_type = []
        self.atom_index = []
        self.nearest_neighbor = np.array([])
        try:
            if type(inp) == dict:
                self.struct_dict = inp              
                self.calculator = inp["calculator"]
            elif type(inp) == str and "\n" in inp:
                lexit("The string func hasn't  been added yet.")
            else:
                with open(inp, "r") as f: 
                    struct_lines = f.readlines()
                self.inp = inp
                self.calculator = "none"
                self.parse_calculator(struct_lines)
                self.parse_struct(struct_lines, self.calculator)
        except:
            lexit("Cannot find the structure file: ", inp)
        self.atoms_sequence = listcross(self.atom_type, self.atom_num)
        self.find_neighbor()
        self.atom_continous_num = self.find_continous_atom_idex()

    def __repr__(self) -> str:
        """"""
        struct_dict = self.return_dict()
        pringinfo = "This is a Class containing structures information:\n"
        pringinfo += "Lattice direct: {}\n".format(struct_dict["lattice"])
        return pringinfo
    
    def find_continous_atom_idex(self):
        return [int(np.sum(np.array(self.atom_num)[:i+1])) for i in range(len(self.atom_type))]
    
    def find_neighbor(self):
        """ """
        neighbor_num = min(5, self.total_atom_num-1)
        self.nearest_neighbor = np.zeros(neighbor_num*2+1)

        for i in range(self.total_atom_num):
            distance = np.array([])
            for j in range(self.total_atom_num):
                if i != j:
                    distance = np.append(distance, 
                                         np.linalg.norm(self.coord_cart[i]-self.coord_cart[j]))
                else:
                    distance = np.append(distance, 1e8)
            distance_sort = np.argsort(distance)

            distance_tmp = np.array([i])
            for j in range(neighbor_num):
                distance_tmp = np.append(distance_tmp, distance_sort[j])
                distance_tmp = np.append(distance_tmp, distance[distance_sort[j]])

            self.nearest_neighbor = np.vstack((self.nearest_neighbor, distance_tmp))
    
    def return_dict(self):
        struct_dict = {}
        struct_dict["title"] = self.title
        struct_dict["atoms_type"] = self.atom_type
        struct_dict["calculator"] = self.calculator
        struct_dict["lattice"] = self.lattice9
        struct_dict["lattice_direct"] = self.lattice6
        struct_dict["atoms_index"] = self.atom_index
        struct_dict["pos_frac"] = self.coord_direct
        struct_dict["pos_cart"] = self.coord_cart
        return struct_dict
    
    def return_recip(self):
        cell = ase.geometry.Cell.fromcellpar([self.lattice6[0], 
                                              self.lattice6[1], 
                                              self.lattice6[2], 
                                              self.lattice6[3], 
                                              self.lattice6[4], 
                                              self.lattice6[5]])
        return cell.reciprocal()*2*np.pi

    def parse_calculator(self, struct_lines):
        """ 
        Parse the calculator from the input file
        Return
        ------
        A string contains calculator. 
        optional [vasp_dyn_direct, vasp_direct, vasp_cart, wien, qein, qeout]
        """
        # vasp
        l7 = struct_lines[7]
        vasp_dyn = False
        if "Selective" in l7 or \
           "selective" in l7 or \
           "S" == l7[0] or \
           "s" == l7[0]:
            vasp_dyn = True
        l8 = struct_lines[8]
        vasp_cart = True
        if "Direct" in l7 or \
           "direct" in l7 or \
           "Direct" in l8 or \
           "direct" in l8 or \
           "D" == l7[0] or \
           "d" == l8[0]:
            vasp_cart = False
        
        if vasp_dyn:
            self.dyn = True
            if vasp_cart:
                self.calculator = "vasp_dyn_cart"
            else:
                self.calculator = "vasp_dyn_direct"
        else:
            if vasp_cart:
                self.calculator = "vasp_cart"
            else:
                self.calculator = "vasp_direct"

        # wien
        if "unit=bohr" in struct_lines[2] or "RELA" in struct_lines[2]:
            self.calculator = "wien"
        # QE
        if "Quantum ESPRESSO" in "".join(struct_lines[:]): 
            self.calculator = "qeout"
        elif "&CONTROL" in "".join(struct_lines[:]):
            self.calculator = "qein"

    def parse_struct(self, struct_lines, calculator):
        struct_dict = {
            "title"      : "",
            "atoms_type" : [],
            "atoms_num"  : [],
            "volume"     : 0,
            "calculator" : ""
        }
        if "vasp" in self.calculator:
            self.title = struct_lines[0][:-1]
            self.ratio = float(struct_lines[1][:-1])
            lattice_x = list(map(float, struct_lines[2].split()))
            lattice_y = list(map(float, struct_lines[3].split()))
            lattice_z = list(map(float, struct_lines[4].split()))
            self.lattice9 = np.array([lattice_x,
                                      lattice_y,
                                      lattice_z])
            self.lattice6 = lattice_conversion(self.lattice9)
            self.atom_type = struct_lines[5].split()
            self.atom_num = np.array(list(map(int, struct_lines[6].split())))
            self.total_atom_num = int(np.sum(self.atom_num))
            self.atom_index = [symbol_map[i] for i in self.atom_type]

            coord_start_line = 8
            if "dyn" in self.calculator:
                coord_start_line = 9
                # check if the dyn symbols present, if not, use "T T T" for all
                if len(struct_lines[9].split()) == 6:
                    self.atom_dyn = np.array([list(map(str, (i[:-1].split())[3:])) \
                                    for i in struct_lines[coord_start_line:coord_start_line+self.total_atom_num]])
                else:
                    self.atom_dyn = np.array([["T", "T", "T"]] * self.total_atom_num)
                self.atom_dyn = self.one_atom_reshape(self.atom_dyn)
                    
            if "direct" in self.calculator:
                self.coord_direct = np.array([list(map(float, 
                                                       (i[:-1].split())[:3])) \
                                                        for i in struct_lines[coord_start_line:coord_start_line+self.total_atom_num]])
                self.coord_direct = self.one_atom_reshape(self.coord_direct)
                self.coord_cart = frac2cart(self.lattice9,
                                            self.coord_direct)
                #self.coord_cart = self.one_atom_reshape(self.coord_cart)
            else:
                self.coord_cart = np.array([list(map(float, 
                                                       (i[:-1].split())[:3])) \
                                                        for i in struct_lines[coord_start_line:coord_start_line+self.total_atom_num]])
                self.coord_cart = self.one_atom_reshape(self.coord_cart)
                self.coord_direct = cart2frac(self.lattice6,
                                              self.coord_cart)
                #self.coord_direct = self.one_atom_reshape(self.cood_direct)
                   
        if calculator == "qeout":
            # get data line: celldm(1); crystal axes; site n.     atom       
            # positions (alat units); number of k points=    12
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

            """total_atoms = int(struct_lines[atom_num_inx].split()[-1])
            struct_dict["lattice"] = np.array([list(map(float, 
                                                        i.split()[-3:])) \
                                                        for i in struct_lines[lattice_inx+1:lattice_inx+4]])
            atoms_pos = ([i.split()[:1][0] for i in struct_lines[atom_position_bgn_inx+1:atom_position_bgn_inx+1+total_atoms]])[:]
            struct_dict["pos_frac"] = np.array([list(map(float, i.split()[-3:])) for i in struct_lines[atom_position_bgn_inx+1:atom_position_bgn_inx+1+total_atoms]])[np.argsort(atoms_pos)]
            atoms_pos = np.array(atoms_pos)[np.argsort(atoms_pos)]
            struct_dict["atoms_type"] = list(dict.fromkeys(atoms_pos))
            struct_dict["atoms_num"] = [list(atoms_pos).count(i) for i in struct_dict["atoms_type"]]
            struct_dict["atoms_sequence"] = listcross(struct_dict["atoms_type"], struct_dict["atoms_num"])"""

            self.total_atom_num = int(struct_lines[atom_num_inx].split()[-1])
            self.lattice9 = np.array([list(map(float,
                                               i.split()[-3:])) \
                                               for i in struct_lines[lattice_inx+1:lattice_inx+4]])
            atoms_pos = ([i.split()[:1][0] \
                          for i in struct_lines[atom_position_bgn_inx+1:\
                                                atom_position_bgn_inx+1+total_atoms]])[:]
            self.coord_direct = np.array([list(map(float, 
                                                   i.split()[-3:])) \
                                                   for i in struct_lines[atom_position_bgn_inx+1:\
                                                                         atom_position_bgn_inx+1+total_atoms]])[np.argsort(atoms_pos)]
            atoms_pos = np.array(atoms_pos)[np.argsort(atoms_pos)]
            self.atom_type = list(dict.fromkeys(atoms_pos))
            self.atom_num = [list(atoms_pos).count(i) \
                             for i in struct_dict["atoms_type"]]
            self.atoms_sequence = listcross(self.atom_type,
                                            self.atom_num)

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
        """ """
        spg_number = [symbol_map[i] for i in self.atoms_sequence]
        spg_tuple = (np.ndarray.tolist(self.lattice9), 
                     np.ndarray.tolist(self.coord_direct), 
                     spg_number)
        lattice, positions, spgnumber = spg.refine_cell(spg_tuple, 
                                                        symprec=symprecision)

        spg_indx = np.argsort(np.array(spgnumber))
        spgnumber = list(np.array(spgnumber)[spg_indx])
        atom_type = list(set(spgnumber))
        atom_num = np.array([spgnumber.count(i) for i in atom_type])
        
        self.lattice9 = np.array(lattice)
        self.coord_direct = np.array(positions)[spg_indx, :]
        self.atom_type = [atom_data[i][1] for i in atom_type]
        self.atom_num = atom_num
    
    def one_atom_reshape(self, v):
        """
        Reshape the Vector for one-atom system 
        """
        if self.total_atom_num == 1:
            return v.reshape(1, 3)
        else:
            return v
    
    def update_struct_dict(self, new_dict):
        """ 
        Update struct dictionary with new and changed structure 

        keys
        ----
        lattice_direct
        """
        try:
            self.lattice9 = new_dict["lattice_direct"]
        except:
            pass
    
    def add_dyn(self, xmin=-1e8, xmax=1e8,
                      ymin=-1e8, ymax=1e8,
                      zmin=-1e8, zmax=1e8):
        """
        Add dyn 
        """
        self.dyn = True
        dyn_sym_w = "T"
        dyn_sym_o = "F"
        coord = self.coord_cart
        tmp = np.full((self.total_atom_num, 3), fill_value=dyn_sym_w)

        if "vasp" in self.calculator:
            if "direct" in self.calculator:
                coord = self.coord_direct
                self.calculator = "vasp_dyn_direct"
            elif "cart" in self.calculator:
                self.calculator = "vasp_dyn_cart"
        elif "qe" in self.calculator:
            dyn_sym_w = "1"
            dyn_sym_o = "0"
            tmp = np.full((self.total_atom_num, 3), fill_value=dyn_sym_w)   

        for i in range(len(coord)):
            if xmax > coord[i, 0] > xmin and \
               ymax > coord[i, 1] > ymin and \
               zmax > coord[i, 2] > zmin:
                tmp[i, 0] = dyn_sym_o
                tmp[i, 1] = dyn_sym_o
                tmp[i, 2] = dyn_sym_o
        self.atom_dyn = tmp

    def change_calculator(self, new_calculator):
        """ 
        update struct calculator option [vasp_direct, wien] 
        """
        self.calculator = new_calculator
    
    def write_struct(self, 
                     fpath:str="./", 
                     fname:str="POSCAR", 
                     calculator:str="none"):
        if calculator == "none":
            calculator = self.calculator
        else:
            if calculator not in ("vasp_direct", "vasp_cart", "qein", 
                                  "vasp_dyn_direct", "vasp_dyn_cart"):
                lexit("""
                Unsupport calculators, now struct writor only supports
                calculators of vasp_direct, vasp_cart, vasp_dyn_direct, 
                vasp_dyn_cart, vasp_dyn, wien, qein
                """)
        
        if "vasp" in calculator:
            with open(fpath+"/"+fname, "w", newline="\n") as f:
                f.write("%s\n" % self.title)
                f.write("%10.5f\n" % self.ratio)
                np.savetxt(fname=f,
                           X=self.lattice9,
                           fmt="%20.10f")
                for i in self.atom_type: 
                    f.write("{:>5s}".format(i))
                f.write("\n")
                np.savetxt(f, self.atom_num.reshape((1,-1)), fmt="%5d")

                # check dyn
                if "dyn" in calculator:
                    if self.dyn:
                        pass
                    else:
                        self.atom_dyn = [["T", "T", "T"]] * self.total_atom_num
                    f.write("Selective dynamics\n")
                    if "direct" in calculator:
                        f.write("Direct\n")
                        np.savetxt(fname=f,
                                   X=np.hstack((self.coord_direct, 
                                                self.atom_dyn)),
                                   fmt="%20s %20s %20s %6s %6s %6s")
                    elif "cart" in calculator:
                        f.write("Cartesian\n")
                        np.savetxt(fname=f,
                                   X=np.hstack((self.coord_cart, 
                                                self.atom_dyn)),
                                   fmt="%20s %20s %20s %6s %6s %6s")
                else:
                    if "direct" in calculator:
                        f.write("Direct\n")
                        np.savetxt(fname=f,
                                   X=self.coord_direct,
                                   fmt="%20s %20s %20s")
                    elif "cart" in calculator:
                        f.write("Cartesian\n")
                        np.savetxt(fname=f,
                                   X=self.coord_cart,
                                   fmt="%20s %20s %20s")

        elif calculator == "wien":
            with open(fpath+"/"+fname, "w", newline="\n") as f:
                f.write("generate by mkits\n")
                f.write("%-4sLATTICE,NONEQUIV.ATOMS%4d%11s\n" % ("P", self.struct_dict["atoms_tot"], "P1"))
                f.write("MODE OF CALC=RELA unit=bohr                                                    \n")
                f.write("%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n" % (self.struct_dict["lattice_direct"][0]/0.529177, self.struct_dict["lattice_direct"][1]/0.529177, self.struct_dict["lattice_direct"][2]/0.529177, self.struct_dict["lattice_direct"][3], self.struct_dict["lattice_direct"][4], self.struct_dict["lattice_direct"][5]))
                for _ in range(self.struct_dict["atoms_tot"]):
                    f.write("ATOM%4d: X=%-10.8f Y=%-10.8f Z=%-10.8f\n" % (_+1, self.struct_dict["pos_frac"][_,0], self.struct_dict["pos_frac"][_,1], self.struct_dict["pos_frac"][_,2]))
                    f.write("          MULT= 1          ISPLIT= 4                                           \n")
                    f.write("%-5s      NPT=  781  R0=%10.8f RMT=%8.5f     Z:%10.5f\n" % (self.struct_dict["atoms_alias"][_], self.struct_dict["atoms_r0"][_], self.struct_dict["atoms_rmt"][_], self.struct_dict["atoms_index"][_]))
                    f.write(local_rot_matrix)
                f.write("   0      NUMBER OF SYMMETRY OPERATIONS                                        \n")
        else:
            pass
                
    def gen_potcar(self, 
                   potpath: str = "./", 
                   fpath: str ="./", 
                   fname: str ="POTCAR"):
        # VASP] generate POTCAR
        atoms = self.atom_type
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
            with open(fpath+"/"+fname, "w", newline="\n") as f:
                f.writelines(potcar)
    
    def del_atom(self,
                 by="element",
                 parameters=""):
        """ """
        if by == "element":
            try:
                parameters = parameters.split()
            except:
                print("Specify the parameter like: Pt Ti O")
            
            # checking
            if isinstance(parameters, list) and all(item in self.atom_type for item in parameters):
                pass
            else:
                lexit("Some element specified cannot find in the structure.")
            
            # find the index
            del_pos_idx = []
            del_atom_idx = []
            atom_continous_num = [0] + self.atom_continous_num
            for i in parameters:
                idx = self.atom_type.index(i)
                del_atom_idx.append(idx)
                del_pos_idx += list(range(atom_continous_num[idx],
                                          atom_continous_num[idx+1]))
            
            # get complement of the del index
            # reserving index
            reserve_atom_idx = list(set(list(range(len(self.atom_type)))) - set(del_atom_idx))
            reserve_pos_idx = list(set(range(self.total_atom_num)) - set(del_pos_idx))
            reserve_atom_idx.sort()
            reserve_pos_idx.sort()  

            # delete term in struct dictionary
            self.atom_index = list_index_by_list(self.atom_index, reserve_atom_idx)
            self.atom_type = list_index_by_list(self.atom_type, reserve_atom_idx)
            self.atom_num = self.atom_num[reserve_atom_idx]
            self.atoms_sequence = listcross(self.atom_type, self.atom_num)
            self.coord_direct = list_index_by_list(self.coord_direct, reserve_pos_idx)
            self.coord_cart = list_index_by_list(self.coord_cart, reserve_pos_idx)
            if self.dyn:
                self.atom_dyn = self.atom_dyn[reserve_pos_idx]
            




class about(object):
    """
    """
    def __init__(self, inp) -> None:
        with open(inp, "r") as f:
            self.inplines = f.readlines()
        
        self.calculator = self.parse_calculator(self.inplines)

    def parse_calculator(self, inplines):
        """
        Parse the calculator for the input files

        Returns:
        --------
        [vasp_xml, vasp_h5, vasp_out]
        """
        