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
    Convert fractional coordinates to cartisian coordinates.
    Performs the linear transformation: r_cart = frac_pos * lattice_matrix
    r = u*a + v*b + w*c

    PARAMETERS
    ----------
    cart_lattice : 3x3 numpy array
        The lattice vectors stored as rows:
        [[ax, ay, az], 
         [bx, by, bz], 
         [cx, cy, cz]]
    fraction_pos : (3,) or (n, 3) numpy array/list
        Fractional coordinates.
    
    RETURNS
    -------
    numpy array
        Cartesian coordinates with the same shape as input fraction_pos.
    """
    lattice = np.array(cart_lattice, dtype=float)
    pos = np.array(fraction_pos, dtype=float)

    return pos @ lattice


def cart2frac(lattice, cart_pos):
    """
    DESCRIPTION
    -----------
    Convert Cartesian coordinates to fractional coordinates.
    Works for both single atom and multiple atoms.
    
    Mathematical formula: r_frac = r_cart * inverse(Lattice_Matrix)

    PARAMETERS
    ----------
    lattice : 3x3 numpy array or list of 6 params
        The lattice matrix (row vectors) or [a, b, c, alpha, beta, gamma].
    cart_pos : (3,) or (N, 3) numpy array
        Cartesian coordinates.
    
    RETURNS
    -------
    numpy array
        Fractional coordinates with the same shape as input cart_pos.
    """
    # 1. 处理晶格矩阵
    lattice = np.array(lattice, dtype=float)
    
    # 如果传入的是 6 个参数，先转换为 3x3 矩阵
    if lattice.shape == (6,) or lattice.size == 6:
        lattice = lattice_conversion(lattice)
    
    # 2. 计算逆矩阵
    # 这是转换的核心：分数坐标 = 笛卡尔坐标 @ 晶格矩阵的逆
    # 这一步只需要做一次，不需要在循环中重复做
    try:
        inv_lattice = np.linalg.inv(lattice)
    except np.linalg.LinAlgError:
        raise ValueError("Lattice matrix is singular (Volume is zero). Check your lattice parameters.")

    # 3. 坐标转换 (矩阵乘法)
    # np.asarray 确保输入是 array，@ 运算符自动处理 (3,) 和 (N, 3) 的情况
    return np.asarray(cart_pos) @ inv_lattice

def cart2frac_single(lattice, cart_pos):
    """
    Just for legacy
    
    :param lattice: Description
    :param cart_pos: Description
    """
    return cart2frac(lattice, cart_pos)


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
    Convert lattice parameters between 3x3 matrix and 6 parameters.
    
    Mode 1: 3x3 Matrix -> [a, b, c, alpha, beta, gamma]
    Mode 2: [a, b, c, alpha, beta, gamma] -> 3x3 Matrix
            (Standard convention: a // x, b in xy plane)
    
    Parameters
    ----------
    give_lattice : np.ndarray
        Either (3, 3) matrix or (6,) vector.
    
    Returns
    -------
    np.ndarray
    """
    give_lattice = np.array(give_lattice, dtype=float)

    # Mode 1: Cartesian Matrix (3x3) -> Parameters (6)
    if give_lattice.shape == (3, 3):
        # 1. Calculate lengths (norms)
        # axis=1 means norm along the row vector
        lengths = np.linalg.norm(give_lattice, axis=1)
        a, b, c = lengths[0], lengths[1], lengths[2]
        
        # 2. Calculate angles
        # alpha: angle between b (row 1) and c (row 2)
        # beta:  angle between a (row 0) and c (row 2)
        # gamma: angle between a (row 0) and b (row 1)
        
        # Use dot product: u . v = |u||v|cos(theta)
        # Clip values to [-1, 1] to avoid numerical errors (e.g. 1.000000002) causing NaN in arccos
        
        cos_alpha = np.dot(give_lattice[1], give_lattice[2]) / (b * c)
        cos_beta  = np.dot(give_lattice[0], give_lattice[2]) / (a * c)
        cos_gamma = np.dot(give_lattice[0], give_lattice[1]) / (a * b)
        
        # Clip to ensure valid arccos range
        cos_alpha = np.clip(cos_alpha, -1.0, 1.0)
        cos_beta  = np.clip(cos_beta,  -1.0, 1.0)
        cos_gamma = np.clip(cos_gamma, -1.0, 1.0)

        alpha = np.degrees(np.arccos(cos_alpha))
        beta  = np.degrees(np.arccos(cos_beta))
        gamma = np.degrees(np.arccos(cos_gamma))
        
        return np.array([a, b, c, alpha, beta, gamma])

    # Mode 2: Parameters (6) -> Cartesian Matrix (3x3)
    elif give_lattice.size == 6:
        # Flatten input in case shape is (6, 1) or similar
        params = give_lattice.flatten()
        a, b, c = params[0], params[1], params[2]
        alpha, beta, gamma = np.radians(params[3:6])
        
        # Standard Convention:
        # a is along x-axis
        # b is in xy-plane
        
        # 1. Vector a
        # a = (a, 0, 0)
        
        # 2. Vector b
        # b_x = b * cos(gamma)
        # b_y = b * sin(gamma)
        
        # 3. Vector c
        # c_x = c * cos(beta)
        # c_y derived from: b.c = b*c*cos(alpha)
        # c_z derived from: c^2 = cx^2 + cy^2 + cz^2
        
        val = (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
        
        # Construct matrix rows
        # Row 1: a vector
        v_a = [a, 0.0, 0.0]
        
        # Row 2: b vector
        v_b = [b * np.cos(gamma), b * np.sin(gamma), 0.0]
        
        # Row 3: c vector
        cz = np.sqrt(1 - np.cos(beta)**2 - val**2)
        v_c = [c * np.cos(beta), 
               c * val, 
               c * cz]
        
        return np.array([v_a, v_b, v_c])

    else:
        # 简单的错误处理，为了兼容你原来的代码风格，可以保留 lexit 或 raise
        raise ValueError("Lattice must be 3x3 matrix or 6 parameters array.")

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
