from mkits.globle import *
import pandas as pd


"""
:func get_data_line_range: get data range from critic2
"""

def get_data_line_range(output, characters):
    """get data range from critic2, start with characters and end with %%"""
    start = 0
    for i in range(len(output)):
        if characters in output[i]: 
            start = i
            break
    for i in range(start+1, len(output)):
        if '%%' in output[i]:
            end = i
            break
        else:
            end = i
    return start, end


def get_auto_cp(output):
    """"""
    start, end = get_data_line_range(output, '%% auto')

    data = [['cp_index', 'coordinates_x', 'coordinates_y', 'coordinates_z', 'Type',
            'f', '|grad_f|', 'del2_f', 'Hessian_x', 'Hessian_y', 'Hessian_z', 'Ellipticity']]
    j = 0
    for i in range(start, end):
        if '+ Critical point' in output[i]:
            data.append([output[i][:-1].split(' ')[-1]])
            j +=1
        if '  Crystallographic coordinates:' in output[i]:
            data[j].append(output[i][:-1].split()[-3])
            data[j].append(output[i][:-1].split()[-2])
            data[j].append(output[i][:-1].split()[-1])
        if '  Type :' in output[i]:
            data[j].append(output[i][:-1].split(' ')[-1])
        if '  Field value (f):' in output[i]:
            data[j].append(output[i][:-1].split(' ')[-1])
        if '  Gradient norm (|grad f|):' in output[i]:
            data[j].append(output[i][:-1].split(' ')[-1])
        if '  Laplacian (del2 f):' in output[i]:
            data[j].append(output[i][:-1].split(' ')[-1])
        if '  Hessian eigenvalues:' in output[i]:
            data[j].append(output[i][:-1].split()[-3])
            data[j].append(output[i][:-1].split()[-2])
            data[j].append(output[i][:-1].split()[-1])
        if '  Ellipticity (l_1/l_2 - 1):' in output[i]:
            data[j].append(output[i][:-1].split(' ')[-1])
    for i in data:
        if len(i) != 12: i.append("nan")

    return data


def get_cpreport(output):
    """get [gtf_kir, vtf_kir, htf_kir]"""
    start, end = get_data_line_range(output, '%% cpreport verylong')
    
    data = [['cp_index', 'gtf_kir', 'vtf_kir', 'htf_kir']]
    j = 0  # set the line_index where write data
    for i in range(start, end):
        if '+ Critical point' in output[i]:
            data.append([output[i][:-1].split(' ')[-1]])
            j += 1
        if 'gtf_kir' in output[i]:
            data[j].append(output[i][:-1].split(' ')[-1])
        if 'vtf_kir' in output[i]:
            data[j].append(output[i][:-1].split(' ')[-1])
        if 'htf_kir' in output[i]:
            data[j].append(output[i][:-1].split(' ')[-1])
    return data


def critc2_extract_data(inputfile:str="more.cro", outfile:str="more.crx"):
    """
    :param inputfile
    :param outfile: string, the name of output file
    :filegen outfile:  
    """

    with open(inputfile, "r") as f:
        inputlines = f.readlines()

    cp_data = get_auto_cp(inputlines)
    qt_data = get_cpreport(inputlines)
    datatowrite = hstack_append_list(qt_data, cp_data)

    # add slash n 
    datatowrite = hstack_append_list(datatowrite, [["\n"]]*len(datatowrite))
    # join list and replace (3,-1) with "bond"
    #datatowrite = ", ".join(item for innerlist in datatowrite for item in innerlist)
    for _ in range(len(datatowrite)):
        datatowrite[_] = ", ".join(item for item in datatowrite[_])
    datatowrite = " ".join(item for item in datatowrite)

    datatowrite = datatowrite.replace("(3,-3)", "artinuclei")
    datatowrite = datatowrite.replace("(3,-1)", "bond")
    datatowrite = datatowrite.replace("(3,1)", "ring")
    datatowrite = datatowrite.replace("(3,3)", "cage")

    with open(outfile, "w") as f:
        f.write(datatowrite)


def critic2_read_crx(inpfile:str="more.crx"):
    """read the output files from function: critc2_extract_data"""
    pd.read_csv(inpfile, sep=",", header=None)