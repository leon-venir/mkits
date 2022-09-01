import numpy as np
import xml.etree.ElementTree as ET
import os
import re
import ase.io
import ase.geometry
import subprocess
from mkits.globle import *
from mkits.sysfc import *
from mkits.database import *


"""
:func vaspxml_parser            : parse the vasprun.xml file
:func vasp_dos_extractor        : extract data from vasprun.xml and create a gnuplot format dos_num.dat
:func vasp_band_extractor       : extract data from EIGNVAL and create a gnuplot format BANDCAR
:func parse_vasp_incar          : parse the INCAR file
:func parse_vasp_incarstring    : parse string-like INCAR
:func parse_vasp_outcar         : parse the OUTCAR file 
:func vasp_potcar_gen           : generate the POTCAR file
:func vasp_kpoints_gen          : generate the KPOINTS file
:func vasp_build_low_dimension  : generate 2-d structure from vasp-format crystal structures
:func vasp_gen_input            : generate input files for VASP valculations
:func vasp_struct_diff          : compare the absolute value of the difference between 2 structure in angstrom)
:func vasp_build_low_dimension  : build 2-dimensinal structure
:func xdatcar_parser            : path to XDATCAR, parse XDATCAR and return a series of struct dict
:func mse_xdatcar               : Mean squared error of structure compared with crystalline structure during md process
:func extract_conv_test         : extract data from convergence test calculation, [input files are generated from
                                : func vasp_gen_input, ie encut, kmesh, ]
:func getvbmcbm                 : extract the valence maximum band and conduction minimum band
:func arbitrary_klist           : Generate arbitrary klist for effective mass calculation
"""


def vaspxml_parser(select_attrib:str, xmlfile:str="vasprun.xml"):
    """
    
    :param select_arrib : optional [tdos, pdos, eigen, final_ene, paramters, incar]
    :return final_ene: the 
    :return klist: array like kpoints list and the the accumulative k length [frac_reciprocal_x,y,z, kpath_length]
    :return highsympoint: 
    :return eigen: return two eigenval arrays with a shape of (kpoint, eigenvalues, occupation)
    :return paramters: parameters_dict key  : NELECT
                                            : NBANDS
    :return incar: incar_dict: 
    """ 
    try:
        vaspxml = ET.parse(xmlfile)
        vasproot = vaspxml.getroot()
        calculation = vasproot.find("calculation")
        atominfo = vasproot.find("atominfo")
    except:
        lexit("Cannot find "+xmlfile+".")

    # parser
    if select_attrib == "dos":
        calculation = vasproot.find("calculation")
        atominfo = vasproot.find("atominfo")

        dos = calculation.find("dos")
        # fermi : dos[0]
        # total dos: dos[1] -> data: dos[1][0][5][0]
        # partial dos: dos[2] -> data: dos[2][0][5][0]
        fermi = dos[0].text
        dos_lines_tot = "#" + "%s" % fermi + "".join(dos[1][0][5][0].itertext())
        # column field: energy s py pz ...
        col_field = ''
        for _ in dos[2][0][3:-2]: col_field += _.text
        # ion num
        ion_num = int(atominfo[0].text)
        # partial dos
        dos_lines_partial = [ "#" + col_field + "%s" % fermi + "".join(dos[2][0][-1][_][0].itertext()) for _ in range(ion_num) ]
        return dos_lines_tot, dos_lines_partial, ion_num
    elif select_attrib == "final_ene":
        calculation = vasproot.find("calculation")
        return float(calculation.find("energy")[-1].text)
    elif select_attrib == "eigen":
        calculation = vasproot.find("calculation")
        eigenvalues = calculation.find("eigenvalues")

        eigen_spin1 = []
        eigen_spin2 = []
        for k in range(len(eigenvalues[0][5][0])):
            for e in range(len(eigenvalues[0][5][0][k])):
                eigen_spin1.append([float(i) for i in eigenvalues[0][5][0][k][e].text.split()])
        return np.array(eigen_spin1).reshape(-1, len(eigenvalues[0][5][0][0]), 2), np.array(eigen_spin2(-1, len(eigenvalues[0][5][0][0]), 2)) if eigen_spin2 else eigen_spin2
    elif select_attrib == "klist":
        kpoints = vasproot.find("kpoints")
        klist = []
        kpoints = kpoints.findall("varray")
        attributes = [_.attrib["name"] for _ in kpoints]
        kpoints = kpoints[attributes.index("kpointlist")]
        for _ in range(len(kpoints)):
            klist.append([float(i) for i in kpoints[_].text.split()])
        klist = np.array(klist)
        reciprocal = vaspxml_parser(select_attrib="cell", xmlfile=xmlfile)["rec_basis_end"]
        kpath_abs = [0]
        cartesian_reciprocal = frac2cart(reciprocal, klist)
        for _ in range(1, len(klist)):
            kpath_abs.append(kpath_abs[_-1]+np.sqrt(np.dot(cartesian_reciprocal[_]-cartesian_reciprocal[_-1], cartesian_reciprocal[_]-cartesian_reciprocal[_-1])))
        return np.hstack((klist, cartesian_reciprocal, np.array(kpath_abs[:]).reshape(-1,1)))
    elif select_attrib == "highsympoints":
        kpoints = vasproot.find("kpoints")
        highsympoints_frac = []
        highsympoints = kpoints.findall("generation")
        attributes = [_.attrib["param"] for _ in highsympoints]
        highsympoints = kpoints[attributes.index("listgenerated")]
        for _ in range(1, len(highsympoints)):
            highsympoints_frac.append([float(i) for i in highsympoints[_].text.split()])
        reciprocal = vaspxml_parser(select_attrib="cell", xmlfile=xmlfile)["rec_basis_end"]
        cartesian_reciprocal = frac2cart(reciprocal, np.array(highsympoints_frac))
        highsympoints_abs = [0]
        for _ in range(1, len(highsympoints_frac)):
            highsympoints_abs.append(highsympoints_abs[_-1]+np.square(np.dot(cartesian_reciprocal[_]-cartesian_reciprocal[_-1], cartesian_reciprocal[_]-cartesian_reciprocal[_-1])))
        return np.array(highsympoints_abs)
    elif select_attrib == "parameters":
        parameters = vasproot.find("parameters")
        parameters = parameters.findall("separator")
        attributes = [_.attrib["name"] for _ in parameters]
        parameters_electronic = parameters[attributes.index("electronic")]
        attributes_electronic = [_.attrib["name"] for _ in parameters_electronic]

        parameters_dict = {}
        parameters_dict["NELECT"] = float(parameters_electronic[attributes_electronic.index("NELECT")].text)
        parameters_dict["NBANDS"] = float(parameters_electronic[attributes_electronic.index("NBANDS")].text)
        return parameters_dict
    elif select_attrib == "incar":
        incars = vasproot.find("incar")
        attributes = [_.attrib["name"] for _ in incars]
        incars_dict = {}
        for _ in attributes:
            incars_dict[_] = incars[attributes.index(_)].text
        return incars_dict
    elif select_attrib == "cell":
        structure = vasproot.findall("structure")
        structure_beg = structure[0]
        structure_end = structure[-1]
        basis_beg = []
        basis_end = []
        for _ in range(len(structure_beg[0][0])):
            basis_beg.append([float(i) for i in structure_beg[0][0][_].text.split()])
        for _ in range(len(structure_end[0][0])):
            basis_end.append([float(i) for i in structure_end[0][0][_].text.split()])

        rec_basis_beg = []
        rec_basis_end = []
        for _ in range(len(structure_beg[0][2])):
            rec_basis_beg.append([float(i) for i in structure_beg[0][2][_].text.split()])
        for _ in range(len(structure_end[0][2])):
            rec_basis_end.append([float(i) for i in structure_end[0][2][_].text.split()])

        position_beg = []
        position_end = []
        for _ in range(len(structure_beg[1])):
            position_beg.append([float(i) for i in structure_beg[1][_].text.split()])
        for _ in range(len(structure_end[0][2])):
            position_end.append([float(i) for i in structure_end[1][_].text.split()])

        volum_beg = float(structure_beg[0][1].text)
        volum_end = float(structure_end[0][1].text)

        return {
            "basis_beg": np.array(basis_beg),
            "basis_end": np.array(basis_end),
            "volum_beg": np.array(volum_beg),
            "volum_end": np.array(volum_end),
            "rec_basis_beg": np.array(rec_basis_beg),
            "rec_basis_end": np.array(rec_basis_end),
            "position_beg": np.array(position_beg),
            "position_end": np.array(position_end)
        }
    else:
        lexit("Cannot find selected attribution of %s.")


def vasp_dos_extractor(xmlfile="vasprun.xml"):
    """  """
    dos_tot, dos_partial, ion_num = vaspxml_parser("dos", xmlfile)
    with open("dos_tot.dat", "w") as f:
        f.write(dos_tot)
    for _ in range(ion_num):
        with open ("dos_ion%d.dat" % _, "w") as f:
            f.write(dos_partial[_])


def vasp_band_extractor_previous(eign='EIGENVAL', poscar="POSCAR", kpoint:str="KPOINTS"):
    """
    get eigen value from EIGNVAL
    electron num; kpoints; states
    :param eign: string, input file EIGENVAL
    :param poscar: string, input file POSCAR
    :param kpoint: string, the input file of KPOINTS for band structure calculation, "line mode" need to be used in KPOINTS
    """
    try:
        kpoints_states_num = np.loadtxt(eign, skiprows=5, max_rows=1, usecols=[0,1,2])
        kpoints_states = np.loadtxt(eign, skiprows=7, usecols=[0,1,2])
        kpoints_states = np.reshape(kpoints_states, (-1, int(kpoints_states_num[2]+1), 3))
    except:
        lexit("Cannot find eigen value: ", eign)
    
    try:
        kpoints_high_sym = np.loadtxt(kpoint, skiprows=4, usecols=[0,1,2])
        kpoints_line = np.loadtxt(kpoint, skiprows=1, max_rows=1)
        #print(kpoints_line)
    except:
        lexit("Cannot find KPOINTS")

    # get reciprocal lattice
    try:
        ucell = ase.io.read(poscar)
        recip = ase.geometry.Cell.reciprocal(ucell.cell)
    except:
        lexit("Cannot find ", poscar)

    # get kpath 
    kpath_abs = np.array([0])
    kpath_high_point = np.array([0])
    for i in range(int(len(kpoints_high_sym)/2)):
        n2 = kpoints_high_sym[2*i+1,0]*recip[0,:] + kpoints_high_sym[2*i+1,1]*recip[1,:] + kpoints_high_sym[2*i+1,2]*recip[2,:]
        n1 = kpoints_high_sym[2*i,0]*recip[0,:] + kpoints_high_sym[2*i,1]*recip[1,:] + kpoints_high_sym[2*i,2]*recip[2,:]
        n1n2 = np.linalg.norm(n2-n1)
        kpath_abs = np.hstack((kpath_abs, np.linspace(kpath_abs[-1], kpath_abs[-1]+n1n2, int(kpoints_line))))
        kpath_high_point = np.append(kpath_high_point, kpath_high_point[-1]+n1n2)
    kpath_abs = kpath_abs[1:]
    kpath_high_point = np.array([kpath_high_point[1:]])


    # write file
    with open('BANDCAR', 'w') as f:
        f.write('# high symmetry points: ')
        np.savetxt(f, kpath_high_point)
        kpath_abs_T = np.reshape(kpath_abs, (-1,1))
        for i in range(len(kpoints_states[0,1:,0])):
            f.write('# %s \n' % str(i+1))
            states_2_write = np.hstack((kpath_abs_T, kpoints_states[:, i+1, 1:]))
            np.savetxt(f, states_2_write)
            f.write('\n\n')


def vasp_band_extractor(fpath:str="./", fname:str="BANDCAR", xmlfile:str="vasprun.xml"):
    """
    get eigen value from EIGNVAL
    electron num; kpoints; states
    :param fpath: string, input file EIGENVAL
    :param fname: string, input file POSCAR
    :param xmlfile: string, the input file of KPOINTS for band structure calculation, "line mode" need to be used in KPOINTS
    """
    kpoint_high_sym = vaspxml_parser(select_attrib="highsympoints", xmlfile=xmlfile)
    kpath = vaspxml_parser(select_attrib="klist", xmlfile=xmlfile)[:, -1]
    eigen1, eigen2 = vaspxml_parser(select_attrib="eigen", xmlfile=xmlfile)

    # write file
    with open(fpath+"/"+fname, "w") as f:
        f.write('# high symmetry points: \n')
        f.write('# %s\n' % ''.join("{:15.8f}".format(i) for i in kpoint_high_sym))
        for _ in range(len(eigen1[0,:,0])):
            f.write('# %s \n' % str(_+1))
            np.savetxt(f, np.hstack((kpath.reshape(-1, 1), eigen1[:,_,:])), fmt="%12.8f" "%10.4f" "%10.4f")
            f.write('\n\n')  


def parse_vasp_incar(fpath="./", fname="INCAR"):
    """ parse the INCAR file and return a dictionary """
    incar_dict = {}
    try:
        linesincar = open(fpath+fname, "r").readlines()
    except:
        lexit("Cannot find INCAR")
    # delete comments
    for i in range(len(linesincar)):
        if "#" in linesincar[i]:
            linesincar[i] = linesincar[i][:linesincar[i].index("#")]
        else:
            linesincar[i] = linesincar[i][:-1]
    linesincar = "\n".join(linesincar)
    linesincar = linesincar.replace(" ", "")
    linesincar = re.split(",|\n", linesincar)
    linesincar = [i for i in linesincar if i]
    for line in linesincar:
        para = re.split("=", line)
        incar_dict[para[0].upper()] = para[-1]
    return incar_dict


def parse_vasp_incarstring(fstring):
    incar_dict = {}
    linesincar = (fstring.replace(" ", "")).split("\n")
    for i in range(len(linesincar)):
        if "#" in linesincar[i]:
            linesincar[i] = linesincar[i][:linesincar[i].index("#")]
        else:
            linesincar[i] = linesincar[i][:]
    linesincar = [i for i in linesincar if i]
    for line in linesincar:
        para = line.split("=")
        incar_dict[para[0]] = para[-1]
    return incar_dict


def parse_vasp_outcar(fpath="./", fname="OUTCAR"):
    outcar_dict = {}
    try:
        linesoutcar = open(fpath+"/"+fname, "r").readlines()
    except:
        lexit("Cannot find OUTCAR.")
    # get optmization information
    outcar_dict["opt_state"] = "none"
    for i in range(1, 50):
        if "reached required accuracy - stopping structural energy minimisation" in linesoutcar[-i]:
            outcar_dict['opt_state'] = "opted"

    # get last total energy
    cmd = 'grep "free  energy   TOTEN  =" %sOUTCAR | tail -1' %fpath
    if subprocess.getoutput(cmd):
        outcar_dict["energy"] = subprocess.getoutput(cmd).split()[-2]
    return outcar_dict


def vasp_potcar_gen(poscar_dict, potpath, fpath="./"):
    """ 
    generate POTCAR
    :param poscar_dict  : dictionary from Class struct
    :param potpath      : path to POTCAR database
    """
    atoms = poscar_dict["atoms_type"]
    potcar = []
    for atom in atoms:
        try:
            potcar += open(str(potpath)+"/"+atom+"/POTCAR", "r").readlines()
        except:
            lexit("The POTCAR of " + atom + " doesn't exit.")
    if os.path.exists(fpath+"/POTCAR"):
        lexit("POTCAT exits, rename it and re-excute the code.")
    else:
        with open(fpath+"/"+"POTCAR", "w") as f:
            f.writelines(potcar)


def vasp_kpoints_gen(poscar_dict, kspacing=0.3, kmesh="none", fpath="./", fname="KPOINTS"):
    """
    Generate kpoints
    :param poscar_dict  : dictionary, structure dict
    :param kspacing     : float
    :param kmesh        : string: odd, even, none(default); generate even or odd kmesh
    """
    
    a1 = np.sqrt(np.dot(poscar_dict["lattice"][0,:], poscar_dict["lattice"][0,:]))
    a2 = np.sqrt(np.dot(poscar_dict["lattice"][1,:], poscar_dict["lattice"][1,:]))
    a3 = np.sqrt(np.dot(poscar_dict["lattice"][2,:], poscar_dict["lattice"][2,:]))
    b1 = 2*np.pi/a1
    b2 = 2*np.pi/a2
    b3 = 2*np.pi/a3
    if kmesh=="none":
        oddeven = -1
    elif kmesh=="odd" or kmesh=="ODD":
        oddeven = 1
    elif kmesh=="even" or kmesh=="EVEN":
        oddeven = 0

    n1 = int(np.max(np.array([1, round_even_odd(b1/(kspacing), oddeven)])))
    n2 = int(np.max(np.array([1, round_even_odd(b2/(kspacing), oddeven)])))
    n3 = int(np.max(np.array([1, round_even_odd(b3/(kspacing), oddeven)])))
    
    with open(fpath+fname, "w") as f:
        f.write("mesh auto\n0\nG\n%s %s %s\n0 0 0" %(str(n1), str(n2), str(n3)))


def vasp_build_low_dimension(poscar="POSCAR", direction="z", vacuum=20, type="add"):
    """
    :param vacuum   :float, vacuum
    :param type     :string add, add vacuum to structure; fit, fit vacuum to vacuum, here vacuum is the lattice parameter""" 
    poscar = struct(poscar)
    struct_dict = poscar.return_dict()
    lattice_para = struct_dict["lattice"]
    fraction = struct_dict["pos_frac"]
    abs_fraction = fraction*np.array([1,1,lattice_para[2,2]])
    if type == "add": pass
    elif type == "fit": vacuum = vacuum - lattice_para[2,2]
    lattice_para[2,2] = lattice_para[2,2]+vacuum
    new_fraction = np.vstack((fraction[:,0], fraction[:,1], (abs_fraction[:,2]+vacuum/2)/lattice_para[2,2])).T

    struct_dict["lattice"] = lattice_para
    struct_dict["pos_frac"] = new_fraction
    poscar.update_struct_dict(struct_dict)
    poscar.write_struct("./", "POSCAR_vacuum.vasp")


def write_incar(incar_dict, fpath="./", fname="INCAR"):
    """ write INCAR """
    with open(fpath+fname, "w") as f:
        for item in incar_dict.keys():
            f.write(item+"="+incar_dict[item]+"\n")


def update_incar(incar_dict, new_dict, new_key):
    """ add new_dict to incar_dict, if the key is duplicated then updates it """
    for key in new_key:
        try:
            if key in incar_dict.keys() and key in new_dict.keys():
                incar_dict[key] = new_dict[key]
            elif key not in incar_dict.keys() and key in new_dict.keys():
                incar_dict.update({key: new_dict[key]})
            else:
                pass
        except:
            pass
    return incar_dict


def vasp_gen_input(dft="scf", potpath="./", poscar="POSCAR", dryrun=False, prec="Normal", wpath="./", execode="mpirun -np $SLURM_NTASKS vasp_std", params="gga=pbelda"):
    func_help = """
    generate inputs for vasp
    :param dft        : string, optional [opt, scf, scfband, conv_encut, conv_kmesh, conv_sigma, xanes, ifc3-phono3py, born]
    :param potpath    : string, the path containing potpaw_LDA.54, potpaw_PBE.54, put all required POTCAR in one directory, and renames it as POTCAR_Ti_sv or POTCAR_Bi
    :param poscar     : 
    :param dryrun     :
    :param prec       : string, optional [low, normal, accurate, l, n, a],
    :param wpath      : working directory
    :param params     : addtitional parameters for calculation, example: gga=ps,gga=rev-vdw-DF2,...
                      : [dft=any]           ->    gga = [pbelda, pbesol, hse06, hsesol, rev-vdW-DF2, optB88, optPBE] default PBE
                      :                     ->    oddeven = [odd, even]
                      :                     ->    kmesh = [0.1, 0.15, ...]
                      : [dft=opt]           ->    mulisif [263], can using with mulprec [Low-Normal-Normal]
                      : [dft=scf]           ->    
                      : [dft=conv_encut]    ->    encut [200-300-400-500], kmesh [0.05-0.1-0.2]
                      : [dft=conv_kmesh]    ->    kspacing [0.05-0.07-0.1-0.15-0.2-0.3]
                      : [dft=xanes]         ->    hole [1s,2s,2p,...]
                      : [dft=born]          ->    
                      : [dft=if3-phono3py]  ->    dim = ["2 2 2", ...]; 
    """
    if dryrun:
        print(func_help)
        lexit(func_help)
    
    poscar = struct(poscar)
    kspacing = 0.07
    if wpath == "./": wpath = os.path.abspath("./")
    wdir = wpath+"/"+dft+"/"
    val_electron = 0
    incar = {}
    
    params = parser_inputpara(params) 
    gga = params["gga"]

    # =================================================================================
    # global setting
    # =================================================================================
    # gen directory
    if not os.path.exists(wpath):
        os.mkdir(wpath)
    if not os.path.exists(wdir):
        os.mkdir(wdir)

    # gen POTCAR and get electron of valence
    if os.path.exists(wdir+"/POTCAR"):
        os.remove(wdir+"/POTCAR")
    potcar_lines = []
    for _ in range(len(poscar.return_dict()["atoms_type"])):
        atom = poscar.return_dict()["atoms_type"][_]
        with open("%s/POTCAR_%s" %(potpath, atom)) as f:
            lines = f.readlines()
        val_electron += float(lines[1][:-1])*poscar.return_dict()["atoms_num"][_]
        potcar_lines += lines
    with open(wdir+"/POTCAR", "w") as f:
        f.writelines(potcar_lines)
    #for _ in poscar.return_dict()["atoms_type"]:
    #    subprocess.run("cat %s/POTCAR_%s* >> %s/POTCAR" %(potpath, _, wdir), shell=True)
    
    def write_runsh(runsh, cmd):
        with open(runsh, "a") as f:
            f.write(cmd)

    def ch_functional(incar):
        """change functional from parameters"""
        update_incar(incar, incar_functionals[gga], list(incar_functionals[gga].keys()))
        return incar
    
    # =================================================================================
    # scf
    # =================================================================================
    if dft == "scf":
        """  """
        incar.update(incar_glob)
        incar.update(incar_scf)
        update_incar(incar, params, ["ENCUT", "PREC"])
        incar = ch_functional(incar)

        vasp_kpoints_gen(poscar.return_dict(), kspacing=float(params["kmesh"]) if "kmesh" in params else 0.15, kmesh=params["oddeven"] if "oddeven" in params else "odd", fpath=wdir, fname="KPOINTS_scf")

        poscar.write_struct(fpath=wdir, fname="POSCAR_init")

        dftgga = dft+"_"+gga
        write_incar(incar, fpath=wdir, fname="INCAR_%s" % dftgga)
        cmd = "# scf calculation\n"
        cmd += "cp INCAR_%s INCAR\n" % dftgga
        cmd += "cp POSCAR_init POSCAR\n"
        cmd += "cp KPOINTS_scf KPOINTS\n"
        cmd += execode + "\n"
        cmd += "cp OUTCAR OUTCAR_%s\n" % dftgga
        cmd += "cp CHGCAR CHGCAR_%s\n" % dftgga
        cmd += "cp WAVECAR WAVECAR_%s\n" % dftgga
        cmd += "cp vasprun.xml vasprun_%s.xml\n" % dftgga
        write_runsh(wdir+"/run.sh", cmd)
        os.chmod(wdir+"/run.sh", 0o775)

    # =================================================================================
    # opt
    # =================================================================================
    if dft == "opt":
        """  """
        incar.update(incar_glob)
        incar.update(incar_opt)
        update_incar(incar, params, ["ENCUT", "PREC", "POTIM"])
        incar = ch_functional(incar)
        # for multi-isif
        if "mulisif" in params:
            if "mulprec" in params:
                optprec = params["mulprec"].split("-")
            else:
                optprec = ["Normal"]*len(params["mulisif"])
            for _ in range(len(params["mulisif"])):
                incar["ISIF"] = params["mulisif"][_]
                incar["PREC"] = optprec[_]

                dftgga = dft+str(_)+"_"+gga+"_isif"+params["mulisif"][_]
                write_incar(incar, fpath=wdir, fname="INCAR_%s_%s_%s_isif%s" %(dftgga, str(_), gga, params["mulisif"][_]))
                
                cmd = "# opt %s isif=%s calculation\n" % (str(_), params["mulisif"][_])
                cmd += "cp INCAR_%s_%s_%s_isif%s INCAR" % (dftgga, str(_), gga, params["mulisif"][_])
                cmd += execode + "\n"
                cmd += "cp OUTCAR OUTCAR_%s_%s_%s_isif%s\n" % (dftgga, str(_), gga, params["mulisif"][_])
                cmd += "cp vasprun.xml vasprun_%s_%s_%s_isif%s.xml\n" % (dftgga, str(_), gga, params["mulisif"][_])
                cmd += "cp CONTCAR POSCAR\n"
                write_runsh(wdir+"/run.sh", cmd)
                os.chmod(wdir+"/run.sh", 0o775)
        else:
            write_incar(incar, fpath=wdir, fname="INCAR_%s_%s_isif%s" %(dft, gga, incar["ISIF"]))
            dftgga = dft+"_"+gga+"_isif"+incar["ISIF"]
            cmd = "# opt isif=%s calculation\n"
            cmd += "cp INCAR_%s INCAR\n" % dftgga
            cmd += execode + "\n"
            cmd += "cp OUTCAR OUTCAR_%s\n" % dftgga
            cmd += "cp vasprun.xml vasprun_%s.xml\n" % dftgga
            cmd += "cp CONTCAR POSCAR\n" 
            write_runsh(wdir+"/run.sh", cmd)
            os.chmod(wdir+"/run.sh", 0o775)

    # =================================================================================
    # convergence test: conv_encut
    # =================================================================================
    if dft == "conv_encut":
        """ generate input files for encut convergence test """
        incar.update(incar_glob)
        incar.update(incar_scf)
        update_incar(incar, params, ["PREC"])
        incar = ch_functional(incar)

        if "encut" in params:
            encut = params["encut"].split("-")
        else:
            encut = ["250", "300", "350", "400", "450", "500", "550", "600"]
        incar["ENCUT"] = "xxx"

        vasp_kpoints_gen(poscar.return_dict(), kspacing=float(params["kmesh"]) if "kmesh" in params else 0.15, kmesh=params["oddeven"] if "oddeven" in params else "odd", fpath=wdir, fname="KPOINTS_scf")


        write_incar(incar, fpath=wdir, fname="INCAR_%s_%s" %(dft, gga))
        poscar.write_struct(fpath=wdir, fname="POSCAR_init")
        dftgga = dft+"_"+gga
        cmd = "# encut convergence test\n"
        cmd += "cp POSCAR_init POSCAR\n"
        cmd += "cp KPOINTS_scf KPOINTS\n"
        cmd += "for encut in %s; do\n" % " ".join(encut)
        cmd += "sed -e s/xxx/$encut/g INCAR_%s > INCAR\n" % dftgga
        cmd += execode + "\n"
        cmd += "cp vasprun.xml vasprun_%s_encut$encut.xml\n" % (dftgga)
        cmd += "done\n"
        write_runsh(wdir+"/run.sh", cmd)
        os.chmod(wdir+"/run.sh", 0o775)

    if dft == "conv_kmesh":
        """ generate input files for kmesh convergence test """
        incar.update(incar_glob)
        incar.update(incar_scf)
        update_incar(incar, params, ["ENCUT", "PREC"])
        incar = ch_functional(incar)

        if "encut" in params:
            kspacing = params["kspacing"].split("-")
        else:
            kspacing = [0.1, 0.12, 0.14, 0.15, 0.2, 0.25, 0.3]
        
        for _ in kspacing:
            vasp_kpoints_gen(poscar.return_dict(), _, kmesh="odd", fpath=wdir, fname="KPOINTS_kconv_k%s" % _)
        
        write_incar(incar, fpath=wdir, fname="INCAR_%s_%s" %(dft, gga))
        poscar.write_struct(fpath=wdir, fname="POSCAR_init")
        dftgga = dft+"_"+gga
        cmd = "# encut convergence test\n"
        cmd += "cp POSCAR_init POSCAR\n"
        cmd += "cp INCAR_%s_%s INCAR\n" %(dft, gga)
        cmd += "for kspacing in %s; do\n" % " ".join(str(_) for _ in kspacing)
        cmd += "cp KPOINTS_kconv_k$kspacing KPOINTS\n"
        cmd += execode + "\n"
        cmd += "cp vasprun.xml vasprun_%s_kmesh$kspacing.xml\n" % (dftgga)
        cmd += "done\n"
        write_runsh(wdir+"/run.sh", cmd)
        os.chmod(wdir+"/run.sh", 0o775)
    
    if dft == "conv_vacuum":
        """ generate input files for vacuum convergence test (2d matertials)"""
    
    if dft == "conv_layer":
        """ generate input files for atomic layers convergence test (2d matertials)"""

    # =================================================================================
    # dos, band, effective mass
    # =================================================================================

    # update setting: "d" for d-elements, "f" for f-elements; set the LMAXMIX=2/4/6
    if max(poscar.return_dict()["atoms_index"]) > 57:
        dfstates_="f"
    elif max(poscar.return_dict()["atoms_index"]) > 20:
        dfstates_="d"
    else:
        dfstates_="n"
    
    if dft == "dos":
        """ """
        print("Make sure WAVECAR_scf and CHGCAR_scf are in the same directory, otherwise, the calculation will be perfomed with dense k-mesh and comsume more time.")


    # =================================================================================
    # XANES
    # =================================================================================
    if dft == "xanes":
        """ generate input files for XANES simulation """
        incar.update(incar_glob)
        incar.update(incar_prop["xanes"])
        update_incar(incar, params, ["ENCUT", "PREC", "SIGMA"])
        incar = ch_functional(incar)
        dftgga = dft+"_"+gga

        # [Ti S]*[2 1]->[Ti Ti S] -> 
        poscar_dict = poscar.return_dict()
        atom_list = listcross(poscar_dict["atoms_type"], poscar_dict["atoms_num"])
        
        # the first atom
        new_atom_type = [atom_list[0]] + list(dict.fromkeys(atom_list[1:]))
        new_atom_num = [1] + [atom_list[1:].count(elem) for elem in list(dict.fromkeys(atom_list[1:]))]
        poscar_dict["atoms_type"] = new_atom_type
        poscar_dict["atoms_num"] = np.array(new_atom_num)
        poscar.update_struct_dict(poscar_dict)
        poscar.write_struct(fpath=wdir, fname="POSCAR_%s_atom%d" % (dftgga, 1))
        potcar_lines = []
        for atom in poscar_dict["atoms_type"]:
            with open("%s/POTCAR_%s" %(potpath, atom)) as f:
                lines = f.readlines()
            val_electron += float(lines[1][:-1])
            potcar_lines += lines
        with open(wdir+"/POTCAR_%s_atom%d" % (dftgga, 1), "w") as f:
            f.writelines(potcar_lines)

        # the last atom
        new_atom_type = list(dict.fromkeys(atom_list[:-1])) + [atom_list[-1]]
        new_atom_num = [atom_list[:-1].count(elem) for elem in list(dict.fromkeys(atom_list[:-1]))] + [1]
        poscar_dict["atoms_type"] = new_atom_type
        poscar_dict["atoms_num"] = np.array(new_atom_num)
        poscar.update_struct_dict(poscar_dict)
        poscar.write_struct(fpath=wdir, fname="POSCAR_%s_atom%d" % (dftgga, len(atom_list)))
        potcar_lines = []
        for atom in poscar_dict["atoms_type"]:
            with open("%s/POTCAR_%s" %(potpath, atom)) as f:
                lines = f.readlines()
            val_electron += float(lines[1][:-1])
            potcar_lines += lines
        with open(wdir+"/POTCAR_%s_atom%d" % (dftgga, len(atom_list)), "w") as f:
            f.writelines(potcar_lines)

        # the others
        for _ in range(1, len(atom_list)-1):
            new_atom_type = list(dict.fromkeys(atom_list[:_])) + [atom_list[_]] + list(dict.fromkeys(atom_list[_+1:]))
            new_atom_num  = [atom_list[:_].count(elem) for elem in list(dict.fromkeys(atom_list[:_]))] + [1] + [atom_list[_+1:].count(elem) for elem in list(dict.fromkeys(atom_list[_+1:]))]
            poscar_dict["atoms_type"] = new_atom_type
            poscar_dict["atoms_num"] = np.array(new_atom_num)
            poscar.update_struct_dict(poscar_dict)
            poscar.write_struct(fpath=wdir, fname="POSCAR_%s_atom%d" % (dftgga, _+1))

            # gen POTCAR and get electron of valence
            potcar_lines = []
            for atom in poscar_dict["atoms_type"]:
                with open("%s/POTCAR_%s" %(potpath, atom)) as f:
                    lines = f.readlines()
                val_electron += float(lines[1][:-1])
                potcar_lines += lines
            with open(wdir+"/POTCAR_%s_atom%d" % (dftgga, _+1), "w") as f:
                f.writelines(potcar_lines)

        vasp_kpoints_gen(poscar.return_dict(), kspacing=float(params["kmesh"]) if "kmesh" in params else 0.15, kmesh=params["oddeven"] if "oddeven" in params else "odd", fpath=wdir, fname="KPOINTS_scf")

        with open(wdir+"/extractor.sh", "w") as f:
            f.write(extract_xanes_vasp)
        os.chmod(wdir+"/extractor.sh", 0o775)
        
        write_incar(incar, fpath=wdir, fname="INCAR_%s" % dftgga)
        cmd = "# xanes for all atoms\n"
        cmd += "cp KPOINTS_scf KPOINTS\n"
        cmd += "for atom in {1..%d}; do\n" % int(np.sum(poscar.return_dict()["atoms_num"]))
        cmd += "cp POSCAR_%s_atom${atom} POSCAR\n" % dftgga
        cmd += "cp POTCAR_%s_atom${atom} POTCAR\n" % dftgga
        cmd += "sed s/xxx/${atom}/g INCAR_%s > INCAR\n" % dftgga
        cmd += execode + "\n"
        cmd += "./extractor.sh\n"
        cmd += "cp OUTCAR OUTCAR_%s_atom${atom}\n" % dftgga
        cmd += "cp vasprun.xml vasprun_%s_atom${atom}.xml\n" % dftgga
        cmd += "cp xanes.dat xanes_%s_atom${atom}.dat\n" % dftgga
        cmd += "done\n"
        write_runsh(wdir+"/run.sh", cmd)
        os.chmod(wdir+"/run.sh", 0o775)
    
    # =================================================================================
    # BORN effective charge
    # =================================================================================
    if dft == "born":
        """calculate BORN effective charge with DFPT method"""
        incar.update(incar_glob)
        incar.update(incar_prop["born"])
        update_incar(incar, params, ["ENCUT", "PREC", "EDIFF"])
        incar = ch_functional(incar)

        vasp_kpoints_gen(poscar.return_dict(), kspacing=float(params["kmesh"]) if "kmesh" in params else 0.15, kmesh=params["oddeven"] if "oddeven" in params else "odd", fpath=wdir, fname="KPOINTS_scf")

        poscar.write_struct(fpath=wdir, fname="POSCAR_init")

        dftgga = dft+"_"+gga
        write_incar(incar, fpath=wdir, fname="INCAR_%s" % dftgga)
        cmd = "# born effective charge calculation\n"
        cmd += "cp INCAR_%s INCAR\n" % dftgga
        cmd += "cp POSCAR_init POSCAR\n"
        cmd += "cp KPOINTS_scf KPOINTS\n"
        cmd += execode + "\n"
        cmd += "cp OUTCAR OUTCAR_%s\n" % dftgga
        cmd += "cp vasprun.xml vasprun_%s.xml\n" % dftgga
        write_runsh(wdir+"/run.sh", cmd)
        os.chmod(wdir+"/run.sh", 0o775)

    # =================================================================================
    # third order force constants with phono3py
    # =================================================================================
    if dft =="ifc3-phono3py":
        """  """
        incar.update(incar_glob)
        incar.update(incar_prop["dfpt"])
        update_incar(incar, params, ["ENCUT", "PREC", "EDIFF"])
        incar["LWAVE"] = ".FALSE."
        incar["LCHAGE"] = ".FALSE."
        incar = ch_functional(incar)

        # check phonon3py installation 
        if "not find" in subprocess.getoutput("which phono3py"):
            lexit("Cannot find phono3py in the PATH, may activate conda before mkits.")
        
        # default setting
        dim = params["dim"] if "dim" in params else "2 2 2"

        vasp_kpoints_gen(poscar.return_dict(), kspacing=float(params["kmesh"]) if "kmesh" in params else 0.15, kmesh=params["oddeven"] if "oddeven" in params else "odd", fpath=wdir, fname="KPOINTS_scf")

        poscar.write_struct(fpath=wdir, fname="POSCAR_init")

        dftgga = dft+"_"+gga
        write_incar(incar, fpath=wdir, fname="INCAR_%s" % dftgga)

        # mkdir inp and out directories
        inpdir = wdir + "/inpdir"
        if not os.path.exists(inpdir):
            os.mkdir(inpdir)
        cmd = "phono3py -c POSCAR_init -d --dim=\"%s\";" % dim
        cmd += "mv POSCAR-* %s" % inpdir.replace(" ", "\ ")
        subprocess.run(cmd, shell=True, cwd=wdir)
        
            
    # =================================================================================
    # finalizing
    # =================================================================================


def struct_diff(structname1, structname2):
    """ 
    compare the absolute value of the difference between structure 1 and 2 
    :param struct1, struct2: string, name of structures; or dictionary """
    if type(structname1) == str:
        struct1 = struct(structname1).return_dict()
    elif type(structname1) == dict:
        struct1 = structname1
    else:
        lexit("Input parameter errors, make sure the type of input is string or dictionary exclusively.")
    if type(structname2) == str:
        struct2 = struct(structname2).return_dict()
    elif type(structname2) == dict:
        struct2 = structname2
    else:
        lexit("Input parameter errors, make sure the type of input is string or dictionary exclusively.")
    
    return np.sum([np.sqrt(np.dot(vec, vec)) for vec in (frac2cart(struct1["lattice"], struct1["pos_frac"]) - frac2cart(struct2["lattice"], struct2["pos_frac"]))])


def vasp_build_low_dimension(poscar="POSCAR", direc="z", vacuum="20", type="add"):
    """
    :param poscar: string, input file
    :param direc: string, optional [z, ]
    :param vacuum: float, vacuum
    :param type: string, optinal [add, fit]
    """
    struct_dict = struct(poscar).return_dict()
    lattice_para = struct_dict["lattice"]
    fraction = struct_dict["pos_frac"]
    if direc == "z":
        abs_fraction = fraction*np.array([1,1,lattice_para[2,2]])
        if type == "add": vacuum = float(vacuum)
        elif type == "fit": vacuum = float(vacuum) - lattice_para[2,2]
        lattice_para[2,2] = lattice_para[2,2]+vacuum
        new_fraction = np.vstack((fraction[:,0], fraction[:,1], (abs_fraction[:,2]+vacuum/2)/lattice_para[2,2])).T

    struct_dict["lattice"] = lattice_para
    struct_dict["pos_frac"] = new_fraction
    return struct(poscar).update_struct_dict(struct_dict)


def xdatcar_parser(xdatcar):
    """ 
    parse XDATCAR and return a series of struct dict 
    :param xdatcar: xdatcar file
    :return : 3-dimensional array, [struct_index, [2-dimensional array of cartesian coordinates]]
    """ 
    with open(xdatcar, "r") as f:
        xdatcar_lines = f.readlines()
        
    coe = float(xdatcar_lines[1].split()[0])
    tot_atoms = int(np.sum(np.loadtxt(xdatcar, skiprows=6, max_rows=1)))
    lattice = np.loadtxt(xdatcar, skiprows=2, max_rows=3) * coe
    tot_struct = 0
    for _ in range(len(xdatcar_lines)-1, 1, -1):
        if "Direct" in xdatcar_lines[_]:
            tot_struct = int(xdatcar_lines[_].split()[-1])
            break

    skiprow_ = 8
    frac_pos_ = np.loadtxt(xdatcar, skiprows=skiprow_, max_rows=tot_atoms)
    xdatcar_cart = frac2cart(lattice, np.where(frac_pos_<0, 1+frac_pos_, frac_pos_)) 
    skiprow_ += tot_atoms+1

    for _ in range(tot_struct-1):
        frac_pos_ = np.loadtxt(xdatcar, skiprows=skiprow_, max_rows=tot_atoms)
        xdatcar_cart = np.dstack((xdatcar_cart, frac2cart(lattice, np.where(frac_pos_<0, 1+frac_pos_, frac_pos_))))
        skiprow_ += tot_atoms+1

    return xdatcar_cart


def extract_xdatcar(xdatcar:str, idx:list, fpath:str="./"):
    """
    parse XDATCAR and write the idx_th structure to "XDATCAR_idx.vasp"
    :param xdatcar: string, xdatcar file name
    :param idx: integer list, the index of 
    """
    with open(xdatcar, "r") as f:
        xdatcar_lines = f.readlines()
    
    coe = float(xdatcar_lines[1].split()[0])
    tot_atoms = int(np.sum(np.loadtxt(xdatcar, skiprows=6, max_rows=1)))
    skiprow = 8
    tot_struct = 0
    for _ in range(len(xdatcar_lines)-1, 1, -1):
        if "Direct" in xdatcar_lines[_]:
            tot_struct = int(xdatcar_lines[_].split()[-1])
            break
    
    frac_pos = np.loadtxt(xdatcar, skiprows=skiprow, max_rows=tot_atoms)
    frac_pos = np.where(frac_pos<0, 1+frac_pos, frac_pos)
    skiprow += tot_atoms+1

    for _ in range(tot_struct-1):
        frac_pos_new = np.loadtxt(xdatcar, skiprows=skiprow, max_rows=tot_atoms)
        frac_pos_new = np.where(frac_pos_new<0, 1+frac_pos_new, frac_pos_new)
        frac_pos = np.dstack((frac_pos, frac_pos_new))
        skiprow += tot_atoms+1
    
    print(frac_pos.shape)

    for _ in idx:
        with open(fpath+"/XDATCAR_"+str(_)+".vasp", "w") as f:
            f.writelines(xdatcar_lines[:7])
            f.write("Direct configuration=  %4d\n" % _)
            np.savetxt(f, frac_pos[:, :, _-1], fmt="%20.16f")
    

def mse_xdatcar(crys_struct:str, xdatcar:str, fpath:str="./"):
    """
    :param crys_struct: the original structure file
    :param xdatcar: XDATCAR file
    """

    crystal_ = struct(crys_struct).return_dict()
    crys_cart = frac2cart(crystal_["lattice"], crystal_["pos_frac"])
    xdatcar = xdatcar_parser(xdatcar)

    mse = np.array([])
    for _ in range(len(xdatcar[0,0,:])):
        mse = np.append(mse, np.sum(np.square(xdatcar[:,:,_]-crys_cart), axis=1).mean())
    
    np.savetxt("%s/mse.dat" % fpath, np.vstack((np.linspace(1, len(mse), len(mse)), mse)).T, fmt="%20.10f")


def extract_conv_test(wdir="./"):
    """
    extract data from convergence test calculation, [input files are generated from func vasp_gen_input, ie encut, kmesh
    :param wdir: working directory
    """
    conv_name_list = []
    for _ in os.listdir(wdir):
        if "vasprun_conv_" in _:
            conv_name_list.append(_)
    # check the type of convergence: vasprun_conv_[kmesh, encut]_[gga]
    if "vasprun_conv_encut_" in conv_name_list[0]:
        ene = [vaspxml_parser("final_ene", wdir+"/"+_) for _ in conv_name_list]
        encut = [float((_.replace(".xml", "")).split("_encut")[-1]) for _ in conv_name_list]
        np.savetxt(wdir+"/encut-ene.dat", (np.vstack((np.array(encut), np.array(ene))).T)[np.argsort(encut)[::-1], :], header="encut(eV) ene(eV)")
        with open(wdir+"/encut-ene.gnu", "w") as f:
            f.write(gnu2dline % ("encut-ene.png", "Convergence test of kinetic energy cutoff", "Kinetic energy cutoff(eV)", "Total energy(eV)", "encut-ene.dat"))
    
    if "vasprun_conv_kmesh_" in conv_name_list[0]:
        ene = [vaspxml_parser("final_ene", wdir+"/"+_) for _ in conv_name_list]
        kmesh = [float((_.replace(".xml", "")).split("_kmesh")[-1]) for _ in conv_name_list]
        np.savetxt(wdir+"/kmesh-ene.dat", (np.vstack((np.array(kmesh), np.array(ene))).T)[np.argsort(kmesh)[::-1], :], header="k-mesh ene(eV)")
        with open(wdir+"/kmesh-ene.gnu", "w") as f:
            f.write(gnu2dline % ("kmesh-ene.png", "Convergence test of k-mesh", "K-Spacing", "Total energy(eV)", "kmesh-ene.dat"))
    

def getvbmcbm(xmlfile:str="vasprun.xml"):
    """
    extract the valence maximum band and conduction minimum band
    :param xmlfile: string, input file
    :return 
    """
    klist = vaspxml_parser(select_attrib="klist", xmlfile=xmlfile)
    eigen1, eigen2 = vaspxml_parser(select_attrib="eigen", xmlfile=xmlfile)
    electron = vaspxml_parser(select_attrib="parameters", xmlfile=xmlfile)["NELECT"]
    incars = vaspxml_parser(select_attrib="incar", xmlfile=xmlfile)

    # check occupation: SO, 
    if "LSORBIT" in incars or "lsorbit" in incars:
        if "T" in incars["LSORBIT"] or "t" in incars["LSORBIT"] or "T" in incars["lsorbit"] or "t" in incars["lsorbit"]:
            vbm_index = int(electron)
            cbm_index = int(vbm_index + 1)
    else:
        vbm_index = int(electron/2)
        cbm_index = int(vbm_index + 1)
    
    vbm_band = eigen1[:,int(vbm_index-1),0]
    cbm_band = eigen1[:,int(cbm_index-1),0]
    vbm_kindex = np.argsort(vbm_band)[-1:-6:-1]
    cbm_kindex = np.argsort(cbm_band)[:5]

    return {
        "vbm_kindex": vbm_kindex,
        "vbm_kpoint": klist[vbm_kindex, :3],
        "vbm_ene": eigen1[vbm_kindex, vbm_index-1, 0],
        "cbm_kindex": cbm_kindex,
        "cbm_kpoint": klist[cbm_kindex, :3],
        "cbm_ene": eigen1[cbm_kindex, cbm_index-1, 0]
    }

    
def arbitrary_klist(kend:str, kmesh:str="7-7-7", fpath:str="./", fname:str="KPOINTS"):
    """
    Generate arbitrary klist for effective mass calculation
    :param kend: string,    "0.5,0.5,0.5,0.01,0.01,0.01" --> the center of the kmesh cuboid with the length of the sides
    :                       "0.0,0.0,0.0,0.0,0.0,0.5  --> two ends of the klist line
    :param kmesh: int, 
    :param fpath: str, 
    :param fname: str, 
    """
    kend = [float(_) for _ in kend.split(",")]
    kmesh = [int(_) for _ in kmesh.split("-")]
    khead = """Automatic generation
%d
Reciprocal lattice
"""

    # cube mesh
    if len(kmesh) == 3:
        x_ = np.linspace(kend[0]-kend[3]/2, kend[0]+kend[3]/2, num=kmesh[0])
        y_ = np.linspace(kend[1]-kend[4]/2, kend[1]+kend[4]/2, num=kmesh[1])
        z_ = np.linspace(kend[2]-kend[5]/2, kend[2]+kend[5]/2, num=kmesh[2])
        xx, yy, zz = np.meshgrid(x_, y_, z_, indexing='ij')
        x, y, z = xx.flatten(), yy.flatten(), zz.flatten()
        
        with open(fpath+"/"+fname, "w") as f:
            f.write(khead % int(kmesh[0]*kmesh[1]*kmesh[2]))
            for _ in range(int(kmesh[0]*kmesh[1]*kmesh[2])):
                f.write("%15.9f%15.9f%15.9f%10.6f\n" % (x[_], y[_], z[_], 1))

    # line mesh
    elif len(kmesh) == 1:
        if kmesh[0] == 7:
            kmesh[0] = 51
        x_ = np.linspace(kend[0], kend[3], num=kmesh[0])
        y_ = np.linspace(kend[1], kend[4], num=kmesh[0])
        z_ = np.linspace(kend[2], kend[5], num=kmesh[0])

        with open(fpath+"/"+fname, "w") as f:
            f.write(khead % int(kmesh[0]))
            for _ in range(int(kmesh[0])):
                f.write("%15.9f%15.9f%15.9f%10.6f\n" % (x_[_], y_[_], z_[_], 1))

    else:
        lexit("Make sure the parameter of kend with correct format: 0.5,0.5,0.5 for cube mesh or 0.0,0.0,0.0,0.0,0.0,0.5 for line mesh")
    
    