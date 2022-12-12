from mkits.globle import *
from mkits.database import *
import ase.io
import numpy as np
import os
import subprocess


"""
:func qe_add_ifpos: add
:func qe_input_parser
:func qe_conv_extractor: extract data from convergence test calculation, [input files are generated from func qe_gen_input, ie conv_encut, conv_kmesh
:func qe_input_lines_parser: get 
:func qe_get_struct: return a standard output for QE input file
"""


def qe_add_ifpos():
    """
    :param wkdir: str
    :param fname: str
    :param qe_inp: str
    :param coords: str : x_1=0,x_2=1,y_1=0,y_2=1,z_1=0,z_2=1
    :param ifpos: str : "000"
    :write 
    write if_pos p
    """



def qe_conv_extractor(wkdir:str="./"):
    """
    extract data from convergence test calculation, [input files are generated from func vasp_gen_input, ie encut, kmesh
    :param wkdir: working directory
    """
    caltype = parser_inputlines(wkdir+"/caldetails")["caltype"]

    conv_name_list = []
    for _ in os.listdir(wkdir):
        if ".out" in _:
            conv_name_list.append(_)
        
    if caltype == "qe_conv_kmesh":
        kpoints = np.array([])
        final_ene = np.array([])
        for conv_f in conv_name_list:
            try:
                qe_out = qe_output_extractor(wkdir+"/"+conv_f)
                kpoints = np.append(kpoints, qe_out["kpoints"])
                final_ene = np.append(final_ene, qe_out["final_ene"])
            except:
                continue
        datatowrite = np.vstack((kpoints, final_ene)).T
        datatowrite = datatowrite[np.argsort(datatowrite[:,0])]
        np.savetxt(fname=wkdir+"/kmesh_ene.dat", X=datatowrite, header="kpoints final_ene(eV)")
        with open(wkdir+"/kmesh_ene.gnu", "w") as f:
            f.write(gnu2dline % ("kmesh_ene.png", "Convergence test of kpoints", "Kpoints", "Total energy(eV)", "kmesh_ene.dat"))
    elif caltype == "qe_conv_ecutwfc":
        ecutwfc = np.array([])
        final_ene = np.array([])
        for conv_f in conv_name_list:
            try:
                qe_out = qe_output_extractor(wkdir+"/"+conv_f)
                ecutwfc = np.append(ecutwfc, qe_out["ecutwfc"])
                final_ene = np.append(final_ene, qe_out["final_ene"])
            except:
                continue
        datatowrite = np.vstack((ecutwfc, final_ene)).T
        datatowrite = datatowrite[np.argsort(datatowrite[:,0])]
        np.savetxt(fname=wkdir+"/ecutwfc_ene.dat", X=datatowrite, header="Ecutwfc(Ry) final_ene(eV)")
        with open(wkdir+"/ecutwfc_ene.gnu", "w") as f:
            f.write(gnu2dline % ("ecutwfc_ene.png", "Convergence test of Kinetic energy cutoff", "Kinetic energy cutoff(Ry)", "Total energy(eV)", "ecutwfc_ene.dat"))
    elif caltype == "qe_conv_degauss":
        degauss = np.array([])
        final_ene = np.array([])
        for conv_f in conv_name_list:
            try:
                if "degauss_" in conv_f:
                    degauss = np.append(degauss, float(conv_f[8:-4]))
                qe_out = qe_output_extractor(wkdir+"/"+conv_f)
                final_ene = np.append(final_ene, qe_out["final_ene"])
            except:
                continue
        datatowrite = np.vstack((degauss, final_ene)).T
        datatowrite = datatowrite[np.argsort(datatowrite[:,0])]
        np.savetxt(fname=wkdir+"/degauss_ene.dat", X=datatowrite, header="Degauss")
        with open(wkdir+"/degauss_ene.gnu", "w") as f:
            f.write(gnu2dline % ("degauss_ene.png", "Convergence test of Kinetic energy cutoff", "Kinetic energy cutoff(Ry)", "Total energy(eV)", "degauss_ene.dat"))


def qe_output_extractor(output_file:str="tmp.out"):
    """
    extract data from QE output file
    :return dictionary: final_ene: final total energy in eV
                      : fermi_ene: fermi energy in eV
                      : ecutwfc: kinetic-energy cutoff in Ry
    """
    with open(output_file, "r") as f:
        output_lines = f.readlines()
    
    # input parameters
    kpoints = int(subprocess.getoutput('grep "number of k points" %s' % output_file).split()[4])
    ecutwfc = float(subprocess.getoutput('grep "kinetic-energy cutoff" %s' % output_file).split()[-2])
    
    # total energy, fermi energy
    final_ene = float(subprocess.getoutput('grep "!    total energy " %s' % output_file).split()[-2]) * uc_ry2ev
    fermi_ene = float(subprocess.getoutput('grep "the Fermi energy is " %s' % output_file).split()[-2])

    # return
    return {
        "kpoints": kpoints,
        "ecutwfc": ecutwfc,
        "final_ene": final_ene,
        "fermi_ene": fermi_ene
    }


def qe_input_parser(qe_inp_file:str, block:str="control"):
    """
    :param block: optional [control, system, electrons, ions, cell]
    """
    with open(qe_inp_file, "r") as f:
        lines = f.readlines()
    
    qe_inp_control_indx = qe_input_lines_parser(lines=lines)

    if "control" in block or "CONTROL" in block:
        return parser_inputlines(lines[qe_inp_control_indx["CONTROL_beg"]:qe_inp_control_indx["CONTROL_end"]], comment_letter="!")
    elif "system" in block or "SYSTEM" in block:
        return parser_inputlines(lines[qe_inp_control_indx["SYSTEM_beg"]:qe_inp_control_indx["SYSTEM_end"]], comment_letter="!")
    elif "electrons" in block or "ELECTRONS" in block:
        return parser_inputlines(lines[qe_inp_control_indx["ELECTRONS_beg"]:qe_inp_control_indx["ELECTRONS_end"]], comment_letter="!")
    elif "ions" in block or "IONS" in block:
        return parser_inputlines(lines[qe_inp_control_indx["IONS_beg"]:qe_inp_control_indx["IONS_end"]], comment_letter="!")
    elif "cell" in block or "CELL" in block:
        return parser_inputlines(lines[qe_inp_control_indx["CELL_beg"]:qe_inp_control_indx["CELL_end"]], comment_letter="!")
    elif "k_points" in block or "K_POINTS" in block:
        return parser_inputlines(lines[qe_inp_control_indx["K_POINTS_beg"]:qe_inp_control_indx["K_POINTS_end"]], comment_letter="!")
    elif "cell_parameters" in block or "CELL_PARAMETERS" in block:
        return parser_inputlines(lines[qe_inp_control_indx["CELL_PARAMETERS_beg"]:qe_inp_control_indx["CELL_PARAMETERS_end"]], comment_letter="!")
    elif "atomic_positions" in block or "ATOMIC_POSITIONS" in block:
        return parser_inputlines(lines[qe_inp_control_indx["ATOMIC_POSITIONS_beg"]:qe_inp_control_indx["ATOMIC_POSITIONS_end"]], comment_letter="!")
    elif "temp" in block or "temp" in block:
        return parser_inputlines(lines[qe_inp_control_indx["_beg"]:qe_inp_control_indx["_end"]], comment_letter="!")
    else:
        lexit("Wrong block name.")


def qe_input_lines_parser(lines:list):
    """
    &CONTROL
    &SYSTEM
    &ELECTRONS
    &IONS
    &CELL
    ATOMIC_SPECIES
    K_POINTS
    CELL_PARAMETERS
    ATOMIC_POSITIONS
    """

    lines_control = {
        "CONTROL_beg": -1,
        "CONTROL_end": -1,
        "SYSTEM_beg": -1,
        "SYSTEM_end": -1,
        "ELECTRONS_beg": -1,
        "ELECTRONS_end": -1,
        "IONS_beg": -1,
        "IONS_end": -1,
        "CELL_beg": -1,
        "CELL_end": -1,
        "ATOMIC_SPECIES_beg": -1,
        "ATOMIC_SPECIES_end": -1,
        "K_POINTS_beg": -1,
        "K_POINTS_end": -1,
        "CELL_PARAMETERS_beg": -1,
        "CELL_PARAMETERS_end": -1,
        "ATOMIC_POSITIONS_beg": -1,
        "ATOMIC_POSITIONS_end": -1
    }

    for i in range(len(lines)):
        if "&CONTROL" in lines[i]:
            lines_control["CONTROL_beg"] = i
            for j in range(i, len(lines)):
                if "/" in lines[j]:
                    lines_control["CONTROL_end"] = j
                    break
        elif "&SYSTEM" in lines[i]:
            lines_control["SYSTEM_beg"] = i
            for j in range(i, len(lines)):
                if "/" in lines[j]:
                    lines_control["SYSTEM_end"] = j
                    break
        elif "&ELECTRONS" in lines[i]:
            lines_control["ELECTRONS_beg"] = i
            for j in range(i, len(lines)):
                if "/" in lines[j]:
                    lines_control["ELECTRONS_end"] = j
                    break
        elif "&IONS" in lines[i]:
            lines_control["IONS_beg"] = i
            for j in range(i, len(lines)):
                if "/" in lines[j]:
                    lines_control["IONS_end"] = j
                    break
        elif "&CELL" in lines[i]:
            lines_control["CELL_beg"] = i
            for j in range(i, len(lines)):
                if "/" in lines[j]:
                    lines_control["CELL_end"] = j
                    break
        elif "ATOMIC_SPECIES" in lines[i]:
            lines_control["ATOMIC_SPECIES_beg"] = i
            for j in range(i, len(lines)):
                if lines[j].split():
                    pass
                else:
                    lines_control["ATOMIC_SPECIES_end"] = j
                    break
        elif "K_POINTS" in lines[i]:
            lines_control["K_POINTS_beg"] = i
            for j in range(i, len(lines)):
                if lines[j].split():
                    pass
                else:
                    lines_control["K_POINTS_end"] = j
                    break
        elif "CELL_PARAMETERS" in lines[i]:
            lines_control["CELL_PARAMETERS_beg"] = i
            for j in range(i, len(lines)):
                if lines[j].split():
                    pass
                else:
                    lines_control["CELL_PARAMETERS_end"] = j
                    break
        elif "ATOMIC_POSITIONS" in lines[i]:
            lines_control["ATOMIC_POSITIONS_beg"] = i
            for j in range(i, len(lines)):
                if lines[j].split():
                    pass
                else:
                    lines_control["ATOMIC_POSITIONS_end"] = j
                    break
    return lines_control


def qe_get_struct(file_inp:str, fractional_position:bool=True):
    crystals = ase.io.read(filename=file_inp)
    ase.io.write(filename="./tmp", images=crystals, format="espresso-in")
    with open("./tmp", "r") as f:
        crystals_lines = f.readlines()
    lines_control = qe_input_lines_parser(crystals_lines)

    qe_struct = {
        "&SYSTEM": crystals_lines[lines_control["SYSTEM_beg"]:lines_control["SYSTEM_end"]+1],
        "ATOMIC_SPECIES": crystals_lines[lines_control["ATOMIC_SPECIES_beg"]:lines_control["ATOMIC_SPECIES_end"]+1],
        "K_POINTS": crystals_lines[lines_control["K_POINTS_beg"]:lines_control["K_POINTS_end"]+1],
        "CELL_PARAMETERS": crystals_lines[lines_control["CELL_PARAMETERS_beg"]:lines_control["CELL_PARAMETERS_end"]+1],
        "ATOMIC_POSITIONS": crystals_lines[lines_control["ATOMIC_POSITIONS_beg"]:lines_control["ATOMIC_POSITIONS_end"]+1],
        "CELL_UNITS": crystals_lines[lines_control["CELL_PARAMETERS_beg"]].split()[1],
        "atoms_num": crystals.get_atomic_numbers()
    }

    # cartesian to fractional
    if fractional_position:
        #cell = crystals.get_cell().round(16)
        cell = crystals.cell.cellpar().round(16)
        cart_pos = crystals.get_positions()
        atom_num = crystals.get_atomic_numbers()
        frac_pos = cart2frac(cart_lattice=cell, cart_pos=cart_pos)

        atomic_position_block = ["ATOMIC_POSITIONS crystal\n"]
        for i in range(len(atom_num)):
            atomic_position_block.append("%6s%20.16f%20.16f%20.16f\n" % (atom_data[atom_num[i]][1], frac_pos[i, 0], frac_pos[i, 1], frac_pos[i, 2]))
    qe_struct["ATOMIC_POSITIONS"] = atomic_position_block

    return qe_struct


def qe_upf_valence_e(upf_file:str):
    """
    :param upf_file: upf file name
    :return int, number of valence electron
    """
    z_val = 0
    with open(upf_file, "r") as f:
        lines = f.readlines()
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


def qe_get_upf(atom_num, atomic_species_block:list, upfpath:str="", wkdir:str="./"):
    """
    :param atom_num: array with atomic number, make sure all upf named with following format: Atom_*
    :param atomic_species_block: str
    :param upfpath: path to 
    :return : atomic_species_block, valence_num, is_uspp
    :write  : write UPF files to working directory
    """
    z_valence = {}
    valence_num = 0
    is_uspp = False

    upflist = os.listdir(upfpath)

    for i in range(len(atomic_species_block)):
        if "None" in atomic_species_block[i]:
            atom = atomic_species_block[i].split()[0]
            for upf in upflist:
                if atom+"_" in upf[:3]: # or str.lower(atom) in upf[:3] or str.upper(atom) in upf[:3]:
                    # read UPF from upf_path
                    with open(upfpath+"/"+upf, "r") as f:
                        upflines = f.readlines()
                    # check if the UPF is USPP
                    for j in range(len(upflines)):
                        if "USPP" in upflines[j] or "ltrasoft" in upflines[j]:
                            is_uspp = True
                            break
                    # wirte UPF file to the working directory
                    with open(wkdir+"/"+upf, "w") as f:
                        f.writelines(upflines) 
                    z_valence[symbol_map[atom]] = qe_upf_valence_e(upfpath+"/"+upf)
                    atomic_species_block[i] = atomic_species_block[i].replace("None", upf)
    for num in atom_num:
        valence_num += z_valence[int(num)]
    return atomic_species_block, valence_num, is_uspp


def qe_in_write(fpath:str="./", fname:str="tmp.in", control_block:dict={}, system_block:dict={}, electrons_block:dict={}, ions_block:dict={}, cell_block:dict={}, atomic_species_block:list=[], kpoints_block:list=[], cell_parameters_block:list=[], atomic_positions_block:list=[], dos_block={}, fcp_block={}, rism_block={}):
    """"""
    if control_block:
        with open(fpath+"/"+fname, "w") as f:
            f.write("&CONTROL\n")
            for item in control_block.keys():
                f.write(item+"="+control_block[item]+"\n")
            f.write("/\n\n")
    else:
        lexit("&CONTROL block should not be blank.")
    if system_block:
        with open(fpath+"/"+fname, "a") as f:
            f.write("&SYSTEM\n")
            for item in system_block.keys():
                f.write(item+"="+system_block[item]+"\n")
            f.write("/\n\n")
    if electrons_block:
        with open(fpath+"/"+fname, "a") as f:
            f.write("&ELECTRONS\n")
            for item in electrons_block.keys():
                f.write(item+"="+electrons_block[item]+"\n")
            f.write("/\n\n")
    if ions_block:
        with open(fpath+"/"+fname, "a") as f:
            f.write("&IONS\n")
            for item in ions_block.keys():
                f.write(item+"="+ions_block[item]+"\n")
            f.write("/\n\n")
    if cell_block:
        with open(fpath+"/"+fname, "a") as f:
            f.write("&CELL\n")
            for item in cell_block.keys():
                f.write(item+"="+cell_block[item]+"\n")
            f.write("/\n\n")
    if atomic_species_block:
        with open(fpath+"/"+fname, "a") as f:
            f.writelines(atomic_species_block)
            f.write("\n")
    if kpoints_block:
        with open(fpath+"/"+fname, "a") as f:
            f.writelines(kpoints_block)
            f.write("\n")
    if cell_parameters_block:
        with open(fpath+"/"+fname, "a") as f:
            f.writelines(cell_parameters_block)
            f.write("\n")
    if atomic_positions_block:
        with open(fpath+"/"+fname, "a") as f:
            f.writelines(atomic_positions_block)
            f.write("\n")
    if dos_block:
        with open(fpath+"/"+fname, "a") as f:
            f.writelines(dos_block)
            f.write("\n")
    if fcp_block:
        with open(fpath+"/"+fname, "a") as f:
            f.writelines(dos_block)
            f.write("\n")
    if rism_block:
        with open(fpath+"/"+fname, "a") as f:
            f.writelines(rism_block)
            f.write("\n")


def qe_kgen(struct_inp:str="cu.cif", kspacing=0.15, oddeven="odd"):
    """
    :return ["K_POINTS automatic\n", "7 7 7 1 1 1\n", "\n"]
    """

    kspacing = float(kspacing)
    crystals = ase.io.read(filename=struct_inp)
    lattice = crystals.cell.cellpar()[:3]
    reciprocal = 2*np.pi/lattice

    if oddeven=="none":
        oddeven = -1
    elif oddeven=="odd" or oddeven=="ODD":
        oddeven = 1
    elif oddeven=="even" or oddeven=="EVEN":
        oddeven = 0
    n1 = int(np.max(np.array([1, round_even_odd(reciprocal[0]/(kspacing), oddeven)])))
    n2 = int(np.max(np.array([1, round_even_odd(reciprocal[1]/(kspacing), oddeven)])))
    n3 = int(np.max(np.array([1, round_even_odd(reciprocal[2]/(kspacing), oddeven)])))
    return ["K_POINTS automatic\n", "%d %d %d  1 1 1\n" % (n1, n2, n3), "\n"]


def write_runsh(runsh, cmd):
        with open(runsh, "a") as f:
            f.write(cmd)
            f.write("\n")


def qe_geninput(calculation:str="scf", wpath:str="./", struct_inp:str="cu.cif", dryrun:bool=False, upf_path:str="./", metal:bool=True, functional:str="pbe", execcode:str='srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS pw.x -nimage 1 -npool 1 -ntg 1 -inp %s > %s', params_file:str="caldetail"):
    """
    Generate input file for QE calculations
    :param       calculation: [scf, nscf, dos, relax, vc-relax, conv_kmesh, conv_ecutwfc]
    :                       : scf
    :                       : nscf
    :                       : dos
    :                       : 
    :                       : relax
    :                       : vc-relax
    :                       : 
    :                       : 
    :param             wpath: 
    :param        struct_inp: 
    :param            dryrun: 
    :param          upf_path: 
    :param             metal: 
    :param        functional: 
    :param          execcode: 
    :param       params_file: 
    :                       : 
    """
    func_help = """
    Generate input file for QE calculations
    --caltype       : optional      opt  ->  
                    :               scf  ->  
                    :              band  ->   
                    :               dos  ->   
    """
    if dryrun:
        print(func_help)
        lexit("Show help with dryrun.")
    
    # 
    control_block = qe_control_block
    system_block = qe_system_block
    electrons_block = qe_electrons_block
    ions_block = qe_ions_block
    cell_block = qe_cell_block
    dos_block = qe_dos_block

    #
    if wpath == "./": wpath = os.path.abspath("./")
    wkdir = wpath+"/"+calculation+"/" 
    if not os.path.exists(wpath):
        os.mkdir(wpath)
    if not os.path.exists(wkdir):
        os.mkdir(wkdir)
    
    # default parameters
    params_default = {
        "kspacing": "0.2",
        "oddeven": "odd"
    }
    try:
        params = parser_inputlines(params_file)
        params_default.update(params)
        for key in params_default.keys():
            if key in qe_control_key:
                control_block[key] = params_default[key]
            elif key in qe_system_key:
                system_block[key] = params_default[key]
            elif key in qe_electrons_key:
                electrons_block[key] = params_default[key]
            elif key in qe_ions_key:
                ions_block[key] = params_default[key]
            elif key in qe_cell_key:
                cell_block[key] = params_default[key]
    except:
        pass

    # structure parameters
    qe_struct = qe_get_struct(file_inp=struct_inp)
    atomic_species_block, valence_num, is_uspp = qe_get_upf(atom_num=qe_struct["atoms_num"], atomic_species_block=qe_struct["ATOMIC_SPECIES"], upfpath=upf_path, wkdir=wkdir)
    qe_struct["ATOMIC_SPECIES"] = atomic_species_block
    system_block.update(parser_inputlines(qe_struct["&SYSTEM"]))
    kpoints_block = qe_kgen(struct_inp=struct_inp, kspacing=params_default["kspacing"])

    # USPP or not, ecutrho = 4*ecutwfc for paw , 10*ecutwfc for uspp
    if is_uspp:
        system_block["ecutrho"] = str(float(system_block["ecutwfc"])*10)
    else:
        system_block["ecutrho"] = str(float(system_block["ecutwfc"])*4)

    # metal or not
    if metal:
        system_block["nbnd"] = str(int(valence_num/2*1.3))
        system_block["occupations"] = '"smearing"'
        system_block["smearing"] = '"gaussian"'
        system_block["degauss"] = "0.01"
    else:
        system_block["nbnd"] = str(int(valence_num/2))
        system_block["occupations"] = '"fixed"'

    # ========================================== #
    # scf
    # ========================================== #
    def qe_scf(fpath:str=wkdir, fname:str="scf.in", control_block:dict=control_block, system_block:dict=system_block, electrons_block:dict=electrons_block, atomic_species_block:list=qe_struct["ATOMIC_SPECIES"], kpoints_block:list=kpoints_block, cell_parameters_block:list=qe_struct["CELL_PARAMETERS"], atomic_positions_block:list=qe_struct["ATOMIC_POSITIONS"]):
        """  """
        control_block["calculation"] = '"scf"'
        qe_in_write(fpath=wkdir, fname=fname, control_block=control_block, system_block=system_block, electrons_block=electrons_block, atomic_species_block=atomic_species_block, kpoints_block=kpoints_block, cell_parameters_block=cell_parameters_block, atomic_positions_block=atomic_positions_block)
    
    def qe_nscf(fpath:str=wkdir, fname:str="nscf.in", control_block:dict=control_block, system_block:dict=system_block, electrons_block:dict=electrons_block, atomic_species_block:list=qe_struct["ATOMIC_SPECIES"], kpoints_block:list=kpoints_block, atomic_positions_block:list=qe_struct["ATOMIC_POSITIONS"]):
        """  """
        control_block["calculation"] = '"nscf"'
        control_block["verbosity"] = '"high"'
        #control_block.pop("restart_mode", None)
        system_block["nosym"] = '.TRUE.'
        electrons_block["conv_thr"] = '1.D-8'

        # you need denser k-mesh
        kpoints_block = qe_kgen(struct_inp=struct_inp, kspacing=float(params_default["kspacing"])*0.4)
        # 
        print("Move all the generated file to the folder of scf calculation.")

        qe_in_write(fpath=wkdir, fname=fname, control_block=control_block, system_block=system_block, electrons_block=electrons_block, atomic_species_block=atomic_species_block, kpoints_block=kpoints_block, atomic_positions_block=atomic_positions_block)
    
    def qe_dos(fpath:str=wkdir, fname:str="dos.in", dos_block:dict=dos_block):
        """  """
        qe_in_write(fpath=wkdir, fname=fname, dos_block=dos_block)


    def qe_vcrelax(control_block:dict=control_block, system_block:dict=system_block, electrons_block:dict=electrons_block, ions_block:dict=ions_block, cell_block:dict=cell_block, atomic_species_block:list=qe_struct["ATOMIC_SPECIES"], kpoints_block:list=kpoints_block, cell_parameters_block:list=qe_struct["CELL_PARAMETERS"], atomic_positions_block:list=qe_struct["ATOMIC_POSITIONS"]):
        """ """
        #
        control_block["calculation"] = '"vc-relax"'
        control_block["forc_conv_thr"] = "1.0d-4"

        qe_in_write(fpath=wkdir, fname="vcrelax.in", control_block=control_block, system_block=system_block, electrons_block=electrons_block, ions_block=ions_block, cell_block=cell_block, atomic_species_block=atomic_species_block, kpoints_block=kpoints_block, cell_parameters_block=cell_parameters_block, atomic_positions_block=atomic_positions_block)
    

    def qe_relax(control_block:dict=control_block, system_block:dict=system_block, electrons_block:dict=electrons_block, ions_block:dict=ions_block, atomic_species_block:list=qe_struct["ATOMIC_SPECIES"], kpoints_block:list=kpoints_block, cell_parameters_block:list=qe_struct["CELL_PARAMETERS"], atomic_positions_block:list=qe_struct["ATOMIC_POSITIONS"]):
        """ """
        #
        control_block["calculation"] = '"relax"'
        control_block["forc_conv_thr"] = "1.0d-4"

        qe_in_write(fpath=wkdir, fname="relax.in", control_block=control_block, system_block=system_block, electrons_block=electrons_block, ions_block=ions_block, atomic_species_block=atomic_species_block, kpoints_block=kpoints_block, cell_parameters_block=cell_parameters_block, atomic_positions_block=atomic_positions_block)


    if calculation == "scf":
        qe_scf()
        cmd = "# scf calculation \n"
        cmd += execcode % ("scf.in", "scf.out") + "\n"
        write_runsh(wkdir+"/run.sh", cmd)
        os.chmod(wkdir+"/run.sh", 0o775)

        cal_detail = "caltype=qe_scf"
        write_runsh(wkdir+"/caldetails", cal_detail)
    
    elif calculation == "nscf":
        qe_nscf()
        cmd = "# nscf calculation \n"
        cmd += execcode % ("nscf.in", "nscf.out") + "\n"
        write_runsh(wkdir+"/run.sh", cmd)
        os.chmod(wkdir+"/run.sh", 0o775)

        cal_detail = "caltype=qe_nscf"
        write_runsh(wkdir+"/caldetails", cal_detail)

    elif calculation == "vc-relax":
        cal_detail = "caltype=qe_vc-relax"
        write_runsh(wkdir+"/caldetails", cal_detail)

        qe_vcrelax()
        cmd = "# opt vc-relax calculation \n"
        cmd += execcode % ("vcrelax.in", "vcrelax.out") + "\n"
        write_runsh(wkdir+"/run.sh", cmd)
        os.chmod(wkdir+"/run.sh", 0o775)

    
    elif calculation == "relax":
        qe_relax()
        cmd = "# opt relax calculation \n"
        cmd += execcode % ("relax.in", "relax.out") + "\n"
        write_runsh(wkdir+"/run.sh", cmd)
        os.chmod(wkdir+"/run.sh", 0o775)

        cal_detail = "caltype=qe_relax"
        write_runsh(wkdir+"/caldetails", cal_detail)
    
    elif calculation == "conv_ecutwfc":
        cal_detail = "caltype=qe_conv_ecutwfc"
        write_runsh(wkdir+"/caldetails", cal_detail)

        system_block["ecutwfc"] = "xxxxx"
        qe_scf()
        cmd = "# conv test calculation \n"
        cmd += "for ecutwfc in 35 40 45 50 55 60 65 70; do\n"
        cmd += "        sed s/xxxxx/$ecutwfc/g scf.in > inp\n"
        cmd += "        " + execcode % ("inp", "ecutwfc_$ecutwfc.out") + "\n"
        cmd += "        rm -rf ./outdir\n"
        cmd += "done\n"
        write_runsh(wkdir+"/run.sh", cmd)
        os.chmod(wkdir+"/run.sh", 0o775)
    
    elif calculation == "conv_degauss" and metal:
        cal_detail = "caltype=qe_conv_degauss"
        write_runsh(wkdir+"/caldetails", cal_detail)

        system_block["degauss"] = "xxxxx"
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
        """"""

