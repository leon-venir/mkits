# -*- coding: utf-8 -*

import numpy as np
import matplotlib.pyplot as plt
import spglib as spg
import os
import sys
import math
import logging
from mkits.database import *


"""
Class
-----

Functions
---------
parse_inputfile:

rmspace:
    Remove [space,tab] characters in strings.
vector_angle:
    Calculate the angle between two vectors.
vector_angle_cclockwise:
    Calculate the angle between two vectors conterclockwise.
listcross:
    Get the crossing product of two lists.
lexit:
    Log an error message and exit.
write2log:
    Log every action and record it.
"""

def find_lines_with_keyword(file_path, 
                            keyword, 
                            returnall=True):
    with open(file_path, 'r', encoding='utf-8') as f:
        for line_number, line in enumerate(f, 1): 
            if keyword in line:
                yield line_number, line.strip()
                if returnall:
                    continue
                else:
                    break


def list_count_elements(listinp, elements):
    _count = 0
    for i in listinp:
        if elements == i:
            _count += 1
    return _count


def extractband(recp_lattice,
                klist,
                eigenvalue,
                thresholdfactor=5):
    """
    DESCRIPTION:
    ------------
    Extract bandstructure from qeout file.

    RETURNS:
    --------
    A gnuplot data.
    """

    klist = frac2cart(
        recp_lattice,
        klist
    )
    # klist
    kvector = np.array([0])
    highpoint = np.array([0])
    threshhold = np.linalg.norm(klist[1] - klist[0])
    for i in range(1, len(klist)):
        
        _ = np.linalg.norm(klist[i] - klist[i-1])
        if _ < threshhold * thresholdfactor:
            kvector = np.append(kvector, _+kvector[-1])
        else:
            kvector = np.append(kvector, kvector[-1])
            highpoint = np.append(highpoint, kvector[-1])
    highpoint = np.append(highpoint, kvector[-1])
    kvector = kvector.reshape((len(klist), 1))
    # gen 
    lines = ["# " 
             + "".join([str( "%10.6f" % _) for _ in highpoint.tolist()])
             + " \n"]

    for i in eigenvalue.T:

        lines += convert_high2writeablelist(
            convert_array2strlist(
                np.hstack((kvector, i.reshape((len(klist), 1))))
            )
        )
        lines += ["\n", "\n"]
        
    return lines


def loadlines(lines, col, dtype=float):
    """
    DESCRIPTION:
    ------------

    PARAMETERS:
    -----------
    
    RETURN:
    -------
    """
    data = np.array(lines[0].split())[col]
    if len(lines) > 1:
        for i in range(1, len(lines)):
            data = np.vstack((data, (np.array(lines[i].split())[col])))
    return np.array(data, dtype=dtype)

def gen_centerklist():
    """
    DESCRIPTION:
    ------------
    Generate a gamma-mode klist with diff weight.

    PARAMETERS:
    -----------
    start
        3 numpy array, 
    
    RETURN:
    -------

    """


def gen_lineklist(start, end, num, weight=1.0):
    """
    DESCRIPTION:
    ------------
    Generate a line-mode klist with same weight.

    PARAMETERS:
    -----------
    start
        3 numpy array, 
    
    RETURN:
    -------
    [[0.0, 0.0, 0.0, 1.0],
     [0.0, 0.0, 0.1, 1.0],
     [0.0, 0.0, 0.2, 1.0],
    ]
    """
    start = np.append(start, weight)
    end = np.append(end, weight)
    return np.linspace(start, end, num)


def write_runsh(runsh, cmd):
    """Write shell run script."""
    with open(runsh, "w", newline="\n") as f:
        f.write(cmd)
    os.chmod(runsh, 0o775)


def dict2lines(
        dicts, 
        assign_sym="="
):
    """ 
    Convert dictionary to writeable string list.
    """
    lines = []
    for item in dicts.keys():
        if dicts[item] == "":
            lines.append(item+"\n")
        else:
            lines.append(
                item +
                assign_sym +
                str(dicts[item]) +
                "\n"
            )
    return lines


def parser_inputpara(inputstring):
    """
    PARAMETERS:
    ----------- 
    inputstring: str
        The input parameters are seperated by ",", 
        and key and attribution are seperated by ":"
    RETURN:
    -------
    dictionary
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


def write_lines(fpath, lines, mode):
    """
    DESCRIPTION:
    ------------
    Write lines to file.

    PARAMETERS:
    -----------
    fpath: str
        Absolute path the file to write.
    lines: list
        A string list.
    """
    with open(fpath, mode) as f:
        f.writelines(lines)


def round_even_odd(num, even_odd):
    """
    DESCRIPTION:
    ------------
    Round input varieties to nearest even number or odd number.
    num              3.5         2.5         2.9
    num/2            1.75        1.25        1.45
    math.modf(num/2) (0.75, 1)   (0.25, 1)   (0.45, 1)

    PARAMETERS:
    -----------
    num: float or int
        Input value.
    even_odd: int
        0 for even, 1 for odd, -1 for round only

    RETURNS:
    --------
        int             
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


def hstack_append_list(list1, list2):
    """
    DESCRIPTION:
    -----------
    Append list2 to list1 in the horizontal direction.

    PARAMETERS:
    -----------
    inp:
        A numpy array.
    fmt:
        Python format 

list1
    [["1","2"],        [["8", "9"]             [["1","2","8","9"],
     ["3","4"],    +    ["10","11"]]    ->      ["3","4","10","11"],
     ["5","7"]]                                 ["5","7"]]
    :param list1 and list2


    RETURNS:
    --------
        A string-like list terminated with \\n
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


def convert_high2writeablelist(inp, endstr="\n"):
    """
    DESCRIPTION:
    -----------
    Convert a list to writeable one-dimension string-like list.

    PARAMETERS:
    -----------
    inp:
        A numpy array.
    endstr:
        Terminated string. 
    RETURNS:
    --------
        A list suiteable for the writelines function.
    """
    if type(inp) == list:
        inp = np.array(inp)
    _list = []

    if inp.ndim == 1:
        inp = inp.reshape(1, len(inp))
    elif inp.ndim == 2:
        pass
    else:
        lexit("Error, functions.py convert_array2strlist not support array higher than 2 dim.")

    for i in range(inp.shape[0]):
        _ = ""
        for j in range(inp.shape[1]):
            _ += inp[i, j]
        _list.append(_+endstr)
    return(_list)


def convert_array2strlist(inp, 
                          fmt="{:20.10f}"):
    """
    DESCRIPTION:
    -----------
    Convert an array or a list to writeable string list.

    PARAMETERS:
    -----------
    inp:
        A numpy array.
    fmt:
        Python format 
    RETURNS:
    --------
        A string-like list terminated with \\n
    """
    if type(inp) == list:
        inp = np.array(inp)
    _list = []

    if inp.ndim == 1:
        inp = inp.reshape(1, len(inp))
    elif inp.ndim == 2:
        pass
    else:
        lexit("Error, functions.py convert_array2strlist not support array higher than 2 dim.")
    for i in range(inp.shape[0]):
        _ = []
        for j in range(inp.shape[1]):
            _.append(fmt.format(inp[i, j]))
        _list.append(_)
    return(_list)


def parse_inputfile(lines, 
                    comment_sym="#", 
                    assign_sym="=", 
                    seprate_sym=","):
    """
    DESCRIPTION:
    -----------
    Parse 

    PARAMETERS:
    -----------
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
    # delete comments and \n
    for i in range(len(lines)):
        lines[i] = lines[i].strip()
        try:
            harshsym_index = lines[i].index(comment_sym)
            lines[i] = lines[i][:harshsym_index]
        except:
            continue
    # expand several keys on the same line
    newlines = []
    for i in range(len(lines)):
        if seprate_sym in lines[i]:
            newlines += lines[i].split(seprate_sym)
        else:
            newlines.append(lines[i])
    # delete blank line and find the keys - values
    for i in range(len(newlines)):
        if newlines[i] == "\n":
            continue
        else:
            try:
                assign_sym_index = newlines[i].index(assign_sym)
                key = newlines[i][:assign_sym_index].strip()
                value = newlines[i][assign_sym_index+1:].strip()
                keywords[key] = value
            except:
                continue
    return keywords


def frac2cart(cart_lattice, fraction_pos):
    """
    DESCRIPTION
    -----------
    Convert fractional coordinates to cartisian coordinates

    PARAMETERS
    ----------
    cartesian lattice
        [[4.45, 0, 0], [-2.22, 3.86, 0], [0, 0, 17.50]]
    fractional coordinates: [[], []]
    
    RETURNS
    -------
    """

    cart_ = np.zeros(3)

    if fraction_pos.shape == (3,):
        fraction_pos = np.vstack((fraction_pos,
                                  np.array([[0,0,0]])))
        for _ in range(len(fraction_pos)):
            cart_ = np.vstack((cart_, 
                               np.sum((cart_lattice.T*fraction_pos[_]).T,
                                      axis=0)))
        return cart_[1:-1][0]
    else:
        for _ in range(len(fraction_pos)):
            cart_ = np.vstack((cart_, 
                               np.sum((cart_lattice.T*fraction_pos[_]).T, 
                                      axis=0)))
        return cart_[1:]


def cart2frac_single(lattice, cart_pos):
    """
    DESCRIPTION
    -----------
    Fractional coordinates to cartesian coordinates
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


def rmspace(string):
    ns=""
    for i in string:
        if(not i == " "):
            ns+=i
    return ns


def rmenter(string):
    ns=""
    for i in string:
        if(not i == "\n"):
            ns+=i
    return ns
def rmtab(string):
    ns=""
    for i in string:
        if(not i == "\t"):
            ns+=i
    return ns
def rmisspace(string):
    ns=""
    for i in string:
        if(not i.isspace):
            ns+=i
    return ns


def vector_angle(vector1, vector2, unit="deg"):
    """
    Calculate the angle between two vectors.
    Parameters
    ----------
    vector1: array or list
        The first vector
    vector2: array or list
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
    return np.rad2deg(angle) if unit == "deg" else angle


def vector_angle_cclockwise(vector1, vector2, unit="deg"):
    """
    Calculate the angle between two vectors conterclockwise.
    Parameters
    ----------
    vector1: array or list
        The first vector
    vector2: array or list
        The second vector
    unit:str
        units, option [rad, deg]
    Returns
    -------
    """
    if isinstance(vector1, list) or isinstance(vector2, list):
        pass
    else:
        vector1 = vector1.tolist()
        vector2 = vector2.tolist()
    if len(vector1) != 2 or len(vector2) != 2:
        lexit("The lenght of the vector has to be 2!")
    vector1 = np.arctan2(*vector1[::-1])
    vector2 = np.arctan2(*vector2[::-1])
    if unit == "deg":
        return np.rad2deg((vector1 - vector2) % (2 * np.pi))


def lattice_conversion(give_lattice):
    """
    Convert lattice parameters 3x3 <-> 6.
    1. give cartesian lattice array; return vectors with angle:
    array([a_x, a_y, a_z], [b_x, b_y, b_z], [c_x, c_y, c_z]) 
    -> array([a, b, c, alpha, beta, gamma])
    2. give vectors with angle; return cartisian lattice array:
    array[a, b, c, alpha, beta, gamma] 
    -> array([a_x, a_y, a_z], [b_x, b_y, b_z], [c_x, c_y, c_z])
    """

    # cartesian -> vector
    if len(give_lattice) == 3:
        lattice_a = np.sqrt(np.dot(give_lattice[0, :], 
                                   give_lattice[0, :]))
        lattice_b = np.sqrt(np.dot(give_lattice[1, :], 
                                   give_lattice[1, :]))
        lattice_c = np.sqrt(np.dot(give_lattice[2, :], 
                                   give_lattice[2, :]))
        angle_alpha = vector_angle(give_lattice[0, :], 
                                   give_lattice[2, :])
        angle_beta  = vector_angle(give_lattice[2, :], 
                                   give_lattice[1, :])
        angle_gamma = vector_angle(give_lattice[1, :], 
                                   give_lattice[0, :])
        return(np.array([lattice_a, lattice_b, lattice_c, 
                         angle_alpha, angle_beta, angle_gamma]))

    # vector -> cartesian
    elif len(give_lattice) == 6:
        # convert arc to degree
        new_lattice = np.hstack((give_lattice[:3], 
                                 give_lattice[3:]/180*np.pi))
        # a in x-axis; b in xy-plane
        vector_x = np.array([new_lattice[0], 0, 0])
        vector_y = np.array([np.cos(new_lattice[5]), 
                             np.sin(new_lattice[5]), 0])*new_lattice[1]
        x = np.cos(new_lattice[3])
        y = (np.linalg.norm(vector_y)*np.cos(new_lattice[4])-
             vector_y[0]*np.cos(new_lattice[3]))/vector_y[1]
        z = np.sqrt(1-x**2-y**2)
        vector_z_direction = np.array([x,y,z])
        vector_z = vector_z_direction / np.linalg.norm(vector_z_direction) * new_lattice[2]
        return np.array([vector_x, vector_y, vector_z])
        
    else:
        lexit("Error, give the lattice with following format: array([a_x, a_y, a_z], [b_x, b_y, b_z], [c_x, c_y, c_z]) or array([a, b, c, alpha, beta, gamma])")


def listcross(list1, list2):
    """
    Get the crossing product of two lists. List2 must be a numeric type list 
    and len(list1) = len(list2)

    Example:
        listcross(['Ge', 'Te'], np.array([2.0, 5.0]))
        listcross(['Ge', 'Te'], [2.0, 5.0])
        output: ['Ge', 'Ge', 'Te', 'Te', 'Te', 'Te', 'Te']
    """
    if (type(list1) != list or type(list2) != list) and type(list2) != np.ndarray:
        lexit("The input lists are not list type.")
    try:
        num_list = [int(i) for i in list2]
    except:
        lexit("Wrong numeric list")
    crosslist = []
    for i in range(len(num_list)):
        crosslist += [list1[i]] * num_list[i]
    return crosslist


# logging
logging.basicConfig(level=logging.DEBUG,
                    filename='output.log',
                    datefmt='%Y/%m/%d %H:%M:%S',
                    format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)
def lexit(message, code=1):
    """Log an error message and exit."""
    logger.error(message)
    sys.exit(code)
def write2log(message):
    """Log every action and record it."""
    logger.info(message)
