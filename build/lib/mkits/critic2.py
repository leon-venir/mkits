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

    # get data from "Analysis of system bonds"
    data2 = [["bcp_ncp", "bcp_end", "bcp_r1_bohr", "bcp_r2_bohr", "bcp_r1or2", "bcp_degree"]]
    bcp_analysis = 0
    bcp_analysis_beg = 0
    bcp_analysis_end = 0
    for i in range(start, end):
        if "* Analysis of system bonds" in output[i]:
            bcp_analysis = i
    for i in range(500):
        if output[i+bcp_analysis][0] != "#" and output[i+bcp_analysis][0] != "*":
            bcp_analysis_beg = i+bcp_analysis
            break
    for i in range(500):
        if not output[i+bcp_analysis].split():
            bcp_analysis_end = i+bcp_analysis
            break

    total_cp = int(data[-1][0])
    for i in range(bcp_analysis_beg, bcp_analysis_end):
        #data2.append([str(int(output[i][0:4])), output[i][4:15].replace(" ", "")[:2] + "-" + output[i][13:22].replace(" ", "")[:2], str(float(output[i][22:31])), str(float(output[i][31:40])), str(float(output[i][40:49])), str(float(output[i][49:57]))])
        dataline = output[i].split()
        data2.append([dataline[0], dataline[1]+dataline[3], dataline[5], dataline[6], dataline[7], dataline[8]])
    ncp_bcp_beg = int(data2[1][0])
    ncp_bcp_end = int(data2[-1][0])
    nan1 = [["nan"] * 6] * (ncp_bcp_beg - 1)
    nan2 = [["nan"] * 6] * (total_cp - ncp_bcp_end)
    data_bond_analysis = [data2[0]]
    data_bond_analysis += nan1
    data_bond_analysis += data2[1:]
    data_bond_analysis += nan2
    data = hstack_append_list(data, data_bond_analysis)

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

    # delete spurious bond critical points, ie bond angle < 160
    for _ in range(len(datatowrite)-1, 0, -1):
        if "(3,-1)" in datatowrite[_]:
            if float(datatowrite[_][-1]) < 160:
                del(datatowrite[_])

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