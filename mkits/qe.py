# -*- coding: utf-8 -*

import numpy as np
import matplotlib.pyplot as plt
import spglib as spg
import os
import sys
import logging
from mkits.functions import *
from mkits.database import *
from mkits.structure import *
from mkits import structure
from mkits import functions
from mkits import database


"""
Class
-----


Functions
---------

parse_qein_block:
    Paser 

"""


def parse_qechg(lines):
    """
    DESCRIPTION:
    ------------
    Parse charge density from the qe hdf5 file.

    PARAMETERS:
    -----------
    lines
        Lines containing all info.

    RETURNS:
    --------
    Return a dictionary.
    """


def qe_extractband(inp, 
                   out="./band.dat", 
                   thresholdfactor=5):
    """
    DESCRIPTION:
    ------------
    Extract bandstructure from qeout file.

    RETURNS:
    --------
    A gnuplot data.
    """
    with open(inp, "r") as f:
        lines = f.readlines()
    qeout = parse_qeout(lines)

    recp_lattice = qeout["reciprocal_parameter"]
    klist = qeout["klist"][-1]
    eigenvalue = qeout["eigenvalue"][-1]

    lines = functions.extractband(
        recp_lattice,
        klist,
        eigenvalue,
        thresholdfactor
    )
        
    if out == "none":
        return lines
    else:
        with open(out, "w") as f:
            f.writelines(lines)


def parse_qeout(lines):
    """
    DESCRIPTION:
    ------------
    Parse info of the qe output file.

    atomic_position: 
        [ relax step
            [ 
                [cartesian x, y, z of atom1]
                [cartesian x, y, z of atom2] 
            ]

            [
                [fractional x, y, z]
                []
            ]
        ]

    PARAMETERS:
    -----------
    lines
        Lines containing all info.

    RETURNS:
    --------
    Return a dictionary.
    """
    results = [
        # value with in keys
        "ntyp", "nat", "ecutwfc", "ecutrho"
        # useful value
        "volume", "total_energy", "vbm", "cbm", "electrons",
        "latticeunit", "klist_num", "eigen_num"
        # value block
        "klist", "eigenvalue", "occupation", "atomic_index"
        # multi 
        "cell_paramater", "atomic_position", "reciprocal_parameter"
    ]
    results = dict.fromkeys(results)

    # start index
    crystal_axes = 0
    reciprocal_axes = 0
    cartesian_axes = 0
    kpoints_index = []
    atomic_position_index = []
    cell_paramater_index = []
    scf_num = 0


    for i in range(len(lines)):
        # value 
        if "lattice parameter" in lines[i]:
            results["latticeunit"] = float(lines[i][:-1].split()[-2])
        elif "unit-cell volume          =" in lines[i]:
            results["volume"] = float(lines[i][:-1].split()[-2])
        elif "number of atoms/cell" in lines[i]:
            results["nat"] = int(lines[i][:-1].split()[-1])
        elif "number of atomic types" in lines[i]:
            results["ntyp"] = int(lines[i][:-1].split()[-1])
        elif "number of electrons" in lines[i]:
            results["electrons"] = float(lines[i][:-1].split()[-1])
        elif "kinetic-energy cutoff" in lines[i]:
            results["ecutwfc"] = float(lines[i][:-1].split()[-2])
        elif "charge density cutoff" in lines[i]:
            results["ecutrho"] = float(lines[i][:-1].split()[-2])
        elif "!    total energy" in lines[i]:
            results["total_energy"] = float(lines[i][:-1].split()[-2])
        elif "number of k points" in lines[i]:
            results["klist_num"] = int(lines[i][:-1].split()[4])
        elif "End of self-consistent calculation" in lines[i]:
            scf_num += 1
        # get block start index
        elif "crystal axes" in lines[i]:
            crystal_axes = i
        elif "reciprocal axes" in lines[i]:
            reciprocal_axes = i
        elif "Cartesian axes" in lines[i]:
            cartesian_axes = i        
        elif "         k =" in lines[i]:
            kpoints_index.append(i)
        elif "ATOMIC_POSITIONS" in lines[i]:
            atomic_position_index.append(i)
        elif "CELL_PARAMETERS" in lines[i]:
            cell_paramater_index.append(i)  
        #elif "End final coordinates" in lines[i]:
        #    break   
        
    # read crystal_axes
    cell_paramater = functions.loadlines(
        lines[crystal_axes+1: crystal_axes+4], [3,4,5]
    ) * results["latticeunit"] * database.uc_bohr2ang
    
    # read reciprocal_axes
    reciprocal_parameter = functions.loadlines(
        lines[reciprocal_axes+1: reciprocal_axes+4], [3,4,5]
    ) * 2 * database.cons_pi / results["latticeunit"] / database.uc_bohr2ang
    results["reciprocal_parameter"] = reciprocal_parameter

    # read first atomic position cartesian_axes
    atomic_position = functions.loadlines(
        lines[cartesian_axes+3: cartesian_axes+results["nat"]+3], 
        [6,7,8]
    )
    _ = functions.loadlines(
        lines[cartesian_axes+3: cartesian_axes+results["nat"]+3], 
        [1],
        str
    )
    _ = [database.symbol_map[_[0]] for _ in _]
    results["atomic_index"] = _

    # read kpoints_index
    klist = np.array(
        [float(lines[kpoints_index[0]][13:20]),
         float(lines[kpoints_index[0]][20:27]),
         float(lines[kpoints_index[0]][27:34])]
    )
    for i in kpoints_index[1:]:
        klist = np.vstack(
            (
                klist,
                np.array(
                    [float(lines[i][13:20]),
                     float(lines[i][20:27]),
                     float(lines[i][27:34])]
                )
            )
        )
    results["klist"] = klist.reshape(-1, results["klist_num"], 3)
    
    # read atomic_position_index
    # [atomic_index, x, y, z]
    atomic_position = np.ones([1, 3])
    if atomic_position_index:
        for i in atomic_position_index:
            _ = functions.loadlines(
                lines[i+1: i+results["nat"]+1], [1,2,3]
            )
            atomic_position = np.vstack((
                atomic_position, _
            ))
    results["atomic_position"] = atomic_position[1:, :].reshape(
        -1, results["nat"], 3
    )
    
    
    # read cell_paramater_index
    if cell_paramater_index:
        for i in cell_paramater_index:
            _ = functions.loadlines(
                lines[i+1: i+4], [0,1,2]
            )
            cell_paramater = np.vstack((
                cell_paramater, _
            ))
    results["cell_paramater"] = cell_paramater.reshape(-1, 3, 3)

    # eigenvalue
    eigennum = " ".join(lines[kpoints_index[0]+1: kpoints_index[1]])
    eigennum = len(eigennum.split())
    results["eigen_num"] = eigennum
    if eigennum%8 == 0:
        eigennum = eigennum//8
    else:
        eigennum = eigennum//8 +1 
    
    eigen = np.array([
        float(_) for _ in " ".join(
            lines[kpoints_index[0]+2: kpoints_index[0]+2+eigennum]
        ).split()
    ])
    for i in kpoints_index[1:]:
        eigen = np.vstack((
            eigen,
            np.array([
                float(_) for _ in " ".join(
                    lines[i+2: i+2+eigennum]
                ).split()
            ])
        ))
    results["eigenvalue"] = eigen.reshape(
        -1, results["klist_num"], results["eigen_num"]
    )

    return results


def qe_conv_kspacing():
    """

    """
    """

    
    elif calculation == "conv_ecutwfc":
        cal_detail = "caltype=qe_conv_ecutwfc"
        write_runsh(wkdir+"/caldetails", cal_detail)

        system_block["ecutwfc"] = "xxxxx"
        system_block["ecutrho"] = "yyyyy"
        qe_scf()
        cmd = "# conv test calculation \n"
        cmd += "for ecutwfc in 35 40 45 50 55 60 65 70; do\n"
        cmd += '        sed -e "s/xxxxx/$ecutwfc/g" -e "s/yyyyy/$((ecutwfc*%s))/g" scf.in > inp\n' % "10" if is_uspp else "4"
        cmd += "        " + execcode % ("inp", "ecutwfc_$ecutwfc.out") + "\n"
        cmd += "        rm -rf ./outdir\n"
        cmd += "done\n"
        write_runsh(wkdir+"/run.sh", cmd)
        os.chmod(wkdir+"/run.sh", 0o775)
    
    elif calculation == "conv_degauss" and metal:
        cal_detail = "caltype=qe_conv_degauss"
        write_runsh(wkdir+"/caldetails", cal_detail)
        
        qe_scf()
        cmd = "# conv test calculation \n"
        cmd += "for degauss in 0.01 0.03 0.05 0.1 0.15 0.2 0.25 0.3; do\n"
        cmd += "        sed s/xxxxx/$degauss/g scf.in > inp\n"
        cmd += "        " + execcode % ("inp", "degauss_$degauss.out") + "\n"
        cmd += "        rm -rf ./outdir\n"
        cmd += "        sleep 5s\n"
        cmd += "done\n"
        write_runsh(wkdir+"/run.sh", cmd)
        os.chmod(wkdir+"/run.sh", 0o775)
    
    elif calculation == "conv_kmesh":
        cal_detail = "caltype=qe_conv_kmesh"
        write_runsh(wkdir+"/caldetails", cal_detail)

        for kmesh in [0.15, 0.2, 0.25, 0.3, 0.4, 0.5]:
            kpoints_block = qe_kgen(struct_inp=struct_inp, kspacing=kmesh)
            qe_scf(fpath=wkdir, fname="kspacing_%4.2f.in" % kmesh, kpoints_block=kpoints_block)
        
        cmd = "# conv test calculation \n"
        cmd += "for kmesh in 0.15 0.20 0.25 0.30 0.40 0.50; do\n"
        cmd += "        cp kspacing_$kmesh.in inp\n"
        cmd += "        " + execcode % ("inp", "kspacing_$kmesh.out") + "\n"
        cmd += "        rm -rf ./outdir\n"
        cmd += "done\n"
        write_runsh(wkdir+"/run.sh", cmd)
        os.chmod(wkdir+"/run.sh", 0o775)
    
    elif calculation == "conv_vacuum":
 
    
    
    """


def qe_autokgen(
        lattice, 
        kspacing=0.15, 
        oddeven="odd"
):
    """
    DESCRIPTION:
    ------------
    Generate auto k-mesh based on the given k-spacing.
    ---> need add more info

    PARAMETERS:
    -----------
    struct_class: class
        Structural object in structure.py.
    kspacing: float
        Distance between each k-point.
    oddeven: str
        Odd- or even mesh of kpoints.
        none, both odd and even mesh is allowed.
        odd, only odd mesh is allowed.
        even, only even mesh is allowed.

    RETURN:
    -------
    "K_POINTS automatic\n", "7 7 7 1 1 1\n", "\n"]
    """

    reciprocal = 2*np.pi/lattice

    if oddeven=="none":
        oddeven = -1
    elif oddeven=="odd" or oddeven=="ODD":
        oddeven = 1
    elif oddeven=="even" or oddeven=="EVEN":
        oddeven = 0
    n1 = int(
        np.max(
            np.array([1, 
                      round_even_odd(reciprocal[0]/(kspacing), 
                                     oddeven)])
        )
    )
    n2 = int(
        np.max(
            np.array([1, 
                      round_even_odd(reciprocal[1]/(kspacing), 
                                     oddeven)])
        )
    )
    n3 = int(
        np.max(
            np.array(
                [1,
                 round_even_odd(reciprocal[2]/(kspacing), 
                                oddeven)]
            )
        )
    )
    return ["K_POINTS automatic\n", "%d %d %d  1 1 1\n" % (n1, n2, n3), "\n"]


def qe_upf_valence_e(lines):
    """
    DESCRIPTION:
    Get the number of valence electrons in UPF.

    PARAMETERS:
    -----------
    lines: list
        A string list from readlines.
    
    RETURN:
    -------
    int, number of valence electron
    """
    z_val = 0

    for line in lines:
        if "Z valence" in line:
            line = line.replace("=", " ")
            line = line.replace(",", " ")
            line = line.replace('"', " ")
            z_val = float(line.split()[0])
        elif "z_valence" in line:
            line = line.replace("=", " ")
            line = line.replace(",", " ")
            line = line.replace('"', " ")
            z_val = float(line.split()[-1])
    return z_val


def qeblock2lines(lines, block):
    """
    DESCRIPTION:
    ------------
    Append QE input block to lines.

    PARAMETERS:
    -----------

    RETURNS:
    --------
    Return a writeable lines containing all blocks.
    """
    if block["qeblock"] == "system":
        del(block["qeblock"])
        lines.append("&SYSTEM\n")
        for _ in block.keys():
            if _ in qe_system_key:
                lines.append(_+"="+str(block[_])+"\n")
        lines.append("/\n\n")
    
    return lines


def addkey2qeinblock(indict, block):
    """
    DESCRIPTION:
    ------------
    Add keys to lines or update keys.

    PARAMETERS:
    -----------
    indict:
        Dictionary contains 
    block:

    RETURNS:
    --------
    Return a writeable blocks containing all keys.
    """


def parse_qein_block(lines):
    """
    DESCRIPTION:
    ------------
    Parse all blocks of the qe input file.

    PARAMETERS:
    -----------
    lines
        Lines containing all 

    RETURNS:
    --------
    Return a dictionary containing all blocks.
    """
    qeinblock = {}
    blockname = []
    start_idx = []
    end_idx = []
    for i in range(len(lines)):
        if rmspace(lines[i])[0] == "&":
            blockname.append(rmspace(lines[i])[1:])
            start_idx.append(i)
    for i in start_idx:
        for j in range(i, len(lines)):
            if rmspace(lines[j])[0] == "/":
                end_idx.append(j)
                break
    for i in range(len(blockname)):
        qeinblock[blockname[i][:-1]] = lines[start_idx[i]+1: end_idx[i]]
    
    # search ATOMIC_POSITIONS, KPOINTS, CELL
    system_block = parse_inputfile(lines=qeinblock["SYSTEM"],
                                   comment_sym="!",
                                   assign_sym="=",
                                   seprate_sym=",")
    try:
        tot_atom = int(system_block["nat"])
    except:
        lexit("The QE input file doesn't contain 'nat'!")
    for i in range(len(lines)):
        if "ATOMIC_POSITIONS" in lines[i].upper():
            qeinblock["ATOMIC_POSITIONS"] = lines[i:i+tot_atom+1]
        elif "K_POINTS" in lines[i].upper():
            qeinblock["K_POINTS"] = lines[i: i+2]
        elif "CELL_PARAMETERS" in lines[i].upper():
            qeinblock["CELL_PARAMETERS"] = lines[i: i+4]
            
    return qeinblock


def qe_geninput(
        calculation="scf", 
        wpath="./",
        wname="pwscf",
        struct_inp="pwscf.in", 
        dryrun=False, 
        upfpath="./", 
        metal=True, 
        mag=False, 
        soc=False,
        execcode='mpirun --bind-to core -np $NTASKS', 
        params="gga=pbe",
        **kwargs
):
    """
    DESCRIPTION:
    ------------
    Generate QE inputs.

    PARAMETERS:
    -----------
    inp, string
        The name of the input structure 

    Attributs:
    params: dict
        Dictionary of 
    is_uspp:
        If USPP upf is used.
    ----------
    """

    func_help = """
    KWARGS:
    -------
    ggau: 
        "Ti=4.2,Cu=6.2"
    dynrange:
        XXX

    PARAMETERS:
        gga=[pbe,lda,pbesol]
        kspacing=0.3  #
        oddeven: [odd, even, none] 
        nonlinear: [F, T]
    
    """
    if dryrun:
        print(func_help)
        lexit("Show help with dryrun.")

    # =================================================================================
    # global setting
    # =================================================================================

    params_default = {
        "kspacing": "0.3",
        "oddeven": "none",
        "nonlinear": "F",
        "nbndfactor": "1.5"
    }

    # import default blocks
    params_default.update(qe_control_block)
    params_default.update(qe_system_block)
    params_default.update(qe_electrons_block)
    params_default.update(qe_ions_block)
    params_default.update(qe_cell_block)

    # set default parameters
    if calculation == "scf":
        params_default["kspacing"] = "0.2"
    elif calculation in ["vcrelax", "relax"]:
        params_default["kspacing"] = "0.35"
    elif calculation in ["dos", "bands"]:
        params_default["kspacing"] = "0.1"
    else:
        lexit("Error: need gga.")
    # update parameters
    params_default.update(
        functions.parser_inputpara(inputstring=params)
    )
    
    # 
    electron_in_band = 2
    if params_default["nonlinear"] == "T":
        electron_in_band = 1

    # generate the working directory
    if wpath == "./": 
        wpath = os.path.abspath("./")
    wkdir = wpath
    if wname != "none":
        wkdir = wpath+"/"+wname+"/"
    if not os.path.exists(wpath):
        os.mkdir(wpath)
    if not os.path.exists(wkdir):
        os.mkdir(wkdir)
    if not os.path.exists(wkdir+"/pseudo"):
        os.mkdir(wkdir+"/pseudo")
    
    # =================================================================================
    # internal functions
    # =================================================================================
    def keysindicts(dicts, key, value):
        # check if a key is in dicts, if not, set to a default value
        if key not in dicts.keys():
            dicts[key] = value
        return dicts
    
    def params2lines(params, block):
        # block: control, system, electrons, ions
        dicts = {}
        if block == "control":
            dicts = {"&CONTROL": ""}
            for _ in params.keys():
                if _ in database.qe_control_key:
                    dicts[_] = params[_]
        elif block == "system":
            dicts = {"&SYSTEM": ""}
            for _ in params.keys():
                if _ in database.qe_system_key:
                    dicts[_] = params[_]
        elif block == "electrons":
            dicts = {"&ELECTRONS": ""}
            for _ in params.keys():
                if _ in database.qe_electrons_key:
                    dicts[_] = params[_]
        elif block == "ions":
            dicts = {"&IONS": ""}
            for _ in params.keys():
                if _ in database.qe_ions_key:
                    dicts[_] = params[_]
        elif block == "cell":
            dicts = {"&CELL": ""}
            for _ in params.keys():
                if _ in database.qe_cell_key:
                    dicts[_] = params[_]
        elif block == "rism":
            dicts = {"&RISM": ""}
            for _ in params.keys():
                if _ in database.qe_system_key:
                    dicts[_] = params[_]
        dicts["/"] = ""
        return functions.dict2lines(dicts)

    # =================================================================================
    # Lattice parameters
    # =================================================================================
    
    # read structure
    qestruct = structure.struct(struct_inp)
    params_default.update(
        qestruct.write_struct(
            calculator="qein",
            dyn=True,
            write2file=False
        )[0]
    )

    dyn = False
    # dynrange or not
    if "dynrange" in kwargs.keys():
        dyn = True
        selec_dyn = kwargs["dynrange"]
        selec_dyn = parser_inputpara(selec_dyn)
        dynrange = {"xmin": -1e8, "ymin": -1e8, "zmin": -1e8,
                    "xmax": 1e8,  "ymax": 1e8,  "zmax": 1e8}
        for key in selec_dyn.keys():
            dynrange[key] = float(selec_dyn[key])

    lattice_lines = qestruct.write_struct(
        calculator="qein",
        dyn=dyn,
        write2file=False
    )[1] 

    # =========================================================================
    # ATOMIC_SPECIES lines 
    # =========================================================================
    # write upf and get electron of valence
    atomic_species_lines = ["ATOMIC_SPECIES\n"]
    val_electron = 0
    is_uspp = False
    atomic_indexs = [
        int(_) for _ in qestruct.position[1:,0]
    ]
    # atomic_index [16, 83]
    atomic_index = list(set(atomic_indexs))
    atomic_num = [
        atomic_indexs.count(_) for _ in atomic_index
    ]
    # UPF
    upflist = os.listdir(upfpath)
    upflist_format = [_.capitalize()[:2] for _ in upflist]
    upflist_format = [_.replace("_", "") for _ in upflist_format]
    upflist_format = [_.replace(".", "") for _ in upflist_format]
    #
    for i in atomic_index:
        atom = atom_data[i][1]
        atomic_mass = atom_data[i][3]
        upf = upflist[upflist_format.index(atom)]
        with open(upfpath+"/"+upf, "r") as f:
            upflines = f.readlines()
        for j in range(len(upflines)):
            if "USPP" in upflines[j] or "ltrasoft" in upflines[j]:
                is_uspp = True
                break
        # write upf
        write_lines(wkdir+"/pseudo/"+upf, upflines, "w")
        # write block
        atomic_species_lines.append(
            "%s%15.5f   %s\n" % (atom, atomic_mass, upf)
        )
        # valence number
        val_electron += int(
            qe_upf_valence_e(upflines) * atomic_num[atomic_index.index(i)]
        )
    
    # uspp or not
    if "ecutrho" not in params_default.keys():
        if is_uspp:
            params_default["ecutrho"] = "%.2f" % (float(params_default["ecutwfc"])*10)
        else:
            params_default["ecutrho"] = "%.2f" % (float(params_default["ecutwfc"])*4)
   
    # metal or not
    if metal:
        params_default["nbnd"] = str(
            int(val_electron/electron_in_band*float(params_default["nbndfactor"]))
        )
        params_default["occupations"] = "'smearing'"
        params_default["smearing"] = "'gaussian'"
        params_default["degauss"] = "0.01"
        params_default["nbndfactor"] = "1.2"
    else:
        if calculation in ["relax", "vcrelax", "scf"]:
            params_default["nbndfactor"] = "1.0"
        else:
            params_default["nbndfactor"] = "1.5"
        params_default["nbnd"] = str(
            int(val_electron/electron_in_band*float(params_default["nbndfactor"]))
        )
        params_default["occupations"] = "'fixed'"
    
    # =========================================================================
    # KPOINTS lines 
    # =========================================================================
    oddeven = "none"
    if "oddeven" in params_default:
        oddeven = params_default["oddeven"]
    kpoints_lines = qe_autokgen(
        lattice=qestruct.lattice6[:3],
        kspacing=float(params_default["kspacing"]),
        oddeven=oddeven
    )
    
    if calculation == "bands":
        kpoints_lines = ["K_POINTS\n"]
        if "klist" in kwargs.keys():
            klist = kwargs["klist"].split(";")
            for i in klist:
                ends = [float(j) for j in i.split(" ")]
                k = functions.gen_lineklist(
                    ends[:3], ends[3:6], int(ends[6])
                )
                k = functions.convert_array2strlist(k)
                k = functions.convert_high2writeablelist(k)               
                kpoints_lines += k
            kpoints_lines.insert(1, str(len(kpoints_lines)-1)+"\n")
        else:
            # need to add 
            functions.lexit("Error, need specify klist\n" + 
                            "klist=0 0 0 0 0 0.5 10;0 0 0.5 0.0 0.5 0.5 10\n")

    # =================================================================================
    # write blocks: 
    # =================================================================================
    
    filein = "/pwscf.in"
    blocklines = []

    # update final setting
    if calculation == "scf":
        params_default["calculation"] = "'scf'"
        for _ in ["control", "system", "electrons"]:
            blocklines += params2lines(params_default, _)
        # write sh
        cmd = "# scf calculation \n"
        cmd += execcode + " %s -nimage 1 -npool 1 -ntg 1 -inp %s > %s" % ("pw.x", "scf.in", "scf.out") + "\n"
        functions.write_runsh(wkdir+"/run_scf.sh", cmd)

    elif calculation == "relax":
        filein = "/relax.in"
        params_default["calculation"] = "'relax'"
        params_default = keysindicts(params_default, "forc_conv_thr", "1.0D-4")
        for _ in ["control", "system", "electrons", 
                  "ions"]:
            blocklines += params2lines(params_default, _)
        # write sh
        cmd = "# scf calculation \n"
        cmd += execcode + " %s -nimage 1 -npool 1 -ntg 1 -inp %s > %s"  % ("pw.x", "relax.in", "relax.out") + "\n"
        functions.write_runsh(wkdir+"/run_relax.sh", cmd)

    elif calculation == "vcrelax":
        filein = "/vcrelax.in"
        params_default["calculation"] = "'vc-relax'"
        params_default = keysindicts(params_default, "forc_conv_thr", "1.0D-4")
        for _ in ["control", "system", "electrons",
                  "ions", "cell"]:
            blocklines += params2lines(params_default, _)
        # write sh
        cmd = "# scf calculation \n"
        cmd += execcode + " %s -nimage 1 -npool 1 -ntg 1 -inp %s > %s"  % ("pw.x", "vcrelax.in", "vcrelax.out") + "\n"
        functions.write_runsh(wkdir+"/run_vcrelax.sh", cmd)
        
    elif calculation == "nscf":
        filein = "/nscf.in"
        params_default["calculation"] = "'nscf'"
        for _ in ["control", "system", "electrons"]:
            blocklines += params2lines(params_default, _)
        # write sh
        cmd = "# nscf calculation \n"
        cmd += execcode + " %s -nimage 1 -npool 1 -ntg 1 -inp %s > %s"  % ("pw.x", "nscf.in", "nscf.out") + "\n"
        functions.write_runsh(wkdir+"/run_nscf.sh", cmd)
    elif calculation == "bands" \
        or calculation == "band":
        filein = "/bands.in"
        del(params_default["restart_mode"])
        params_default["calculation"] = "'bands'"
        for _ in ["control", "system", "electrons"]:
            blocklines += params2lines(params_default, _)
        # write sh
        cmd = "# band calculation \n"
        cmd += execcode + " %s -nimage 1 -npool 1 -ntg 1 -inp %s > %s"  % ("pw.x", "bands.in", "bands.out") + "\n"
        functions.write_runsh(wkdir+"/run_bands.sh", cmd)
    
    # =================================================================================
    # write lines: 
    # =================================================================================
    functions.write_lines(
        fpath=wkdir+filein,
        lines=blocklines,
        mode="w"
    )    
    functions.write_lines(
        fpath=wkdir+filein,
        lines=atomic_species_lines+lattice_lines+kpoints_lines,
        mode="a"
    )

    # gga+u or not
    if "ggau" in kwargs.keys():
        ujlines = ["HUBBARD {ortho-atomic}\n"]
        
        functions.write_lines(
            fpath=wkdir+filein,
            lines=ujlines,
            mode="a"
        )
    
    


    
