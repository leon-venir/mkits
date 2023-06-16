import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import os
import re
import ase.io
import ase.geometry
import ase.dft.kpoints
import subprocess
from mkits.globle import *
from mkits.sysfc import *
from mkits.database import *


"""
:func vasp_show_incar           : show several INCARs 
:func gapshifter                : shift gap
:func vasp_opt_2d               : 
:func vasp_gen_IBZKPT           : generate k-points list
:func xml_block                 : search the block with specific block_name and return the block or the index
:func vasp_split_IBZKPT         : split IBZKPT for huge kpoints calculations
:func vasp_dp_effect_mass       : 
:func vasp_defomation_potential : generate structures for deformation potential calculation
:func gen_gnu_bandcar           : generate gnuplot scripts 
:func add_selective_dyn         : add the 4-, 5-, 6th columns for selective dynamic POSCAR
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


def vasp_show_incar(*arguments:str):
    """
    """
    maxrow = 0
    maxchr = 0

    for arg in arguments:
        with open(arg, "r") as f:
            lines = f.readlines()
        if len(lines) > maxrow:
            maxrow = len(lines)
        for line in lines:
            if len(line) >  maxchr:
                maxchr = len(line)
    maxchr += 2

    totlines = []
    for arg in arguments:
        with open(arg, "r") as f:
            lines = f.readlines()
        if len(lines) < maxrow:
            totlines += lines + ["\n"]*int(maxrow - len(lines))
        else:
            totlines += lines
    
    for i in range(len(totlines)):
        totlines[i] = totlines[i].replace("\n", " "*(maxchr - len(totlines[i])))
    
    #
    for i in range(maxrow):
        print("|".join(totlines[i::maxrow]))


def gapshifter(code:str, shift:float, inp:str, out:str="shifted.o", absolute:bool=True):
    """
    shift to/or shift specific value
    :param code: str [vasp, wien2k]
    :param absolute: bool, default True shift gap to specific value; otherwise shift specific value
    :param shift: float, the value of gap 
    :param inp: string, the input file
    :param out: string, the output file
    """
    with open(inp, "r") as f:
        inplines = f.readlines()

    if code == "vasp" and ".xml" in inp:
        # get the number of valence bands
        valenceband = 0
        incar = vaspxml_parser(select_attrib="incar", xmlfile=inp)
        nelect = vaspxml_parser(select_attrib="parameters", xmlfile=inp)["NELECT"]
        
        if "LSORBIT" in incar:
            if "T" in incar["LSORBIT"] or "t" in incar["LSORBIT"]:
                valenceband = int(nelect)
            else:
                valenceband = int(nelect/2)
        else:
            valenceband = int(nelect/2)
        
        if absolute:
            energygap = getvbmcbm(inp)["ene_gap"]
            realshift = shift - energygap

        else:
            realshift = shift
        
        # get eigenvalue block and 
        eigen_beg = 0
        eigen_end = 0
        for i in range(len(inplines)):
            if "<eigenvalues>" in inplines[i]:
                eigen_beg = i
            elif "</eigenvalues>" in inplines[i]:
                eigen_end = i

        with open(out, "w") as f:
            # the part before eigenvalues
            f.writelines(inplines[:eigen_beg])

            # eigenvalue part
            eigen_indx = 1
            for i in range(eigen_beg, eigen_end+1):
                if "kpoint" in inplines[i]:
                    eigen_indx = 1
                    f.write(inplines[i])
                elif "<r>" in inplines[i] and "<r>" in inplines[i]:
                    if eigen_indx <= valenceband:
                        f.write(inplines[i])
                        eigen_indx += 1
                    else:
                        f.write("       <r>%10.4f%s" % (float(inplines[i][11:20])+realshift, inplines[i][20:]))
                        eigen_indx += 1
                else:
                    f.write(inplines[i])

            # the part after eigenvalues
            f.writelines(inplines[eigen_end+1:])
        
        
    elif code == "vasp" and ".hdf5" in inp:
        pass
    elif code == "wien2k":
        pass


def vasp_opt_2d():
    """
    optimization -> position(ISIF=2) -> ab/ac -> aob/aoc
    :param CONTCAR_to_opt: 
    :param opt_choice:
    :param step: the step of shape changes in percentage
    :param runsh: 
    """


def vasp_gen_IBZKPT(wkdir:str="./", fname:str="IBZKPT_666", poscar:str="POSCAR", kmesh:str="6,6,6", write_file:bool=True):
    """
    generate irreducible k-list and the weights
    :param wkdir:
    :param fname: 
    :param poscar:
    :param kmesh:
    :write
    """
    ase.dft.kpoints.ibz_points()



def xml_block(xmlfile, block_name:str, block_indx:int=-1, block:bool=True):
    """
    :param xmlfile: str: absolute path to xml file or a list of xml lines
    :param block_name: the target block name to search { <eigenvalues> }, the ending gives as { </eigenvalues> }
    :param block_indx: default return the last block with an index of -1
    :param block: return the block or just return the index of the target block
    :return list: line list 
    search the block with specific block_name and return the block or the index
    """
    if isinstance(xmlfile, str):
        with open(xmlfile, "r") as f:
            lines = f.readlines()
    elif isinstance(xmlfile, list):
        lines = xmlfile

    beg_indx = []
    end_indx = []
    for i in range(len(lines)):
        if "<%s" % block_name in lines[i]:
            beg_indx.append(i)
        if "</%s" % block_name in lines[i]:
            end_indx.append(i)
    
    if block:
        return lines[beg_indx[block_indx]:end_indx[block_indx]+1]
    else:
        return beg_indx[block_indx], end_indx[block_indx]


def vasp_split_IBZKPT(wkdir:str="./", fname:str="vasprun_merged.xml", ibzkpt:str="IBZKPT", kperibzkpt:int=30, split_merge:str="split"):
    """
    This function can split a huge k-points calculation to several smaller ones and merge 
    the results into an unitary vasprun_merged.xml file.

    :param wkdir            : absolute path to the working folder
    :param fname            : the name of the merged xml file
    :param ibzkpt           : the IBZKPT file
    :param kperibzkpt       : the number of k-points in each sub-calculation
    :param split_merge      : optional [split,merge] 
                            : "split" the IBZKPT to several files, make sure you have all 
                            :         the INPUTs, including INCAR with ICHARG=11, POSCAR, 
                            :         POTCAR, CHGCAR_scf_pbelda and WAVECAR_scf_pbelda (optional)  
                            :         from previous SCF calculation. Please modify the  
                            :         run_IBZKPT_split.sh file with your own settings. 
                            : "merge" the splitted vasprun.xml into an unitary one.
    """
    if split_merge == "split" or split_merge == "s":
        print("Make sure you have all the required INPUTs, including INCAR with ICHARG=11, POSCAR, POTCAR, CHGCAR_scf_pbelda and WAVECAR_scf_pbelda (optional) from previous SCF calculation. Please modify the run_IBZKPT_split.sh file with your own settings.")

        with open(ibzkpt, "r") as f:
            lines = f.readlines()
        total_k = int(lines[1].split()[0])
        start_line_indx = 3
        
        group_num = total_k//kperibzkpt

        for i in range(1, group_num+1):
            with open(wkdir+"/IBZKPT_%d" % i, "w") as f:
                f.write("Automatically generated mesh\n")
                f.write("%8s\n" % kperibzkpt)
                f.write("Reciprocal lattice\n")
                f.writelines(lines[start_line_indx:start_line_indx+kperibzkpt])
            start_line_indx += kperibzkpt
        
        if total_k%kperibzkpt != 0:
            with open("%s/IBZKPT_%d" % (wkdir, int(group_num+1)), "w") as f:
                f.write("Automatically generated mesh\n")
                f.write("%8d\n" % int(total_k%kperibzkpt))
                f.write("Reciprocal lattice\n")
                f.writelines(lines[start_line_indx:])
            group_num += 1
        
        cmd = "# split dos\n"
        cmd += "for ibz in {1..%d}; do \n" % group_num
        cmd += "        cp WAVECAR_scf_pbelda WAVECAR\n"
        cmd += "        cp CHGCAR_scf_pbelda CHGCAR\n"
        cmd += "        cp IBZKPT_$ibz KPOINTS \n"
        cmd += "        mpirun -np $SLURM_NTASKS vasp_std\n"
        cmd += "        mv vasprun.xml vasprun_IBZKPT_$ibz.xml\n"
        cmd += "done\n"
        with open("%s/run_IBZKPT_split.sh" % wkdir, "w") as f:
            f.write(cmd)
        os.chmod("%s/run_IBZKPT_split.sh" % wkdir, 0o775)
        
    elif split_merge == "merge" or split_merge == "m":
        xmlfile_list = subprocess.getoutput("ls %s/vasprun_IBZKPT_*.xml" % wkdir).split()
        kpoints_beg, kpoints_end = xml_block(xmlfile=xmlfile_list[0], block_name="kpoints", block=False)
        eigen_beg, eigen_end = xml_block(xmlfile=xmlfile_list[0], block_name="eigenvalues", block=False)
        #print(kpoints_beg, kpoints_end)
        #print(eigen_beg, eigen_end)
        #print(xmlfile_list)
        with open(xmlfile_list[0], "r") as f:
            lines = f.readlines()
        block_before_kpoints = lines[:kpoints_beg]
        block_kpoints2eigen = lines[kpoints_end+1:eigen_beg]
        block_after_eigen = lines[eigen_end+1:]

        block_kpoints = []
        block_weight = []
        block_eigen = []

        for i in range(len(xmlfile_list)):
            xmlfile = wkdir + "/vasprun_IBZKPT_%d.xml" % (i+1)
            kpoints_lines = xml_block(xmlfile=xmlfile, block_name="kpoints", block=True)
            block_kpoints += xml_block(xmlfile=kpoints_lines, block_name="varray", block_indx=0, block=True)[1:-1]
            block_weight += xml_block(xmlfile=kpoints_lines, block_name="varray", block_indx=1, block=True)[1:-1]

            # eigenvalue
            eigenvalues_lines = xml_block(xmlfile=xmlfile, block_name="eigenvalues", block=True)
            #print(eigenvalues_lines)
            block_eigen += eigenvalues_lines[9:-4]
        
        # add index in block_eigen
        kpoint_index = 1
        for i in range(len(block_eigen)):
            if '<set comment="kpoint' in block_eigen[i]:
                block_eigen[i] = '      <set comment="kpoint %d">\n' % kpoint_index
                kpoint_index += 1

        total_block = []
        total_block += block_before_kpoints
        total_block += [' <kpoints>\n', '  <varray name="kpointlist" >\n']
        total_block += block_kpoints
        total_block += ['  </varray>\n']
        total_block += ['  <varray name="weights" >\n']
        total_block += block_weight
        total_block += ['  </varray>\n', ' </kpoints>\n']
        total_block += block_kpoints2eigen
        total_block += ['  <eigenvalues>\n', '   <array>\n', '    <dimension dim="1">band</dimension>\n', '    <dimension dim="2">kpoint</dimension>\n', '    <dimension dim="3">spin</dimension>\n', '    <field>eigene</field>\n', '    <field>occ</field>\n', '    <set>\n', '     <set comment="spin 1">\n']
        total_block += block_eigen
        total_block += ['     </set>\n', '    </set>\n', '   </array>\n', '  </eigenvalues>\n']
        total_block += block_after_eigen


        with open(wkdir+"/"+fname, "w") as f:
            f.writelines(total_block)        

    else:
        lexit("Only ")


def vasp_dp_effect_mass(fpath:str="./", xmlfile:str="vasprun.xml", carrier:str="electron", specifyrange=False):
    """
    1, Read a xml file with a line-like kpoints, and calculate the effective mass with parabolic approximation.
    :param fpath:
    :param xmlfile:
    :param carrier: str, optional [electron, hole]
    """
    vbm_index, cbm_index = vasp_cbmvbm_index(xmlfile="%s/%s" % (fpath, xmlfile))
    print("The number of valence bands top is %d." % vbm_index)
    print("The number of conduction bands bottom is %d." % cbm_index)
    if carrier == "electron":
        # get eigen states, kpath
        eigen, eigen2 = vaspxml_parser(select_attrib="eigen", xmlfile="%s/%s" % (fpath, xmlfile))
        kpath = vaspxml_parser(select_attrib="klist", xmlfile="%s/%s" % (fpath, xmlfile))[:, -1]
        cbm_band = eigen[:, cbm_index, 0]
        # get best-fitting range
        if specifyrange:
            bestfitting = np.array([[specifyrange[0], specifyrange[1], 0]])
        else:
            bestfitting = best_polyfit_range(kpath, cbm_band, 2, 5, 20)
        
        # data to write and units conversion
        print(bestfitting)
        dat_to_write = np.vstack((kpath[int(bestfitting[0,0]):int(bestfitting[0,1])]/uc_ang2m, cbm_band[int(bestfitting[0,0]):int(bestfitting[0,1])]*uc_ev2j))
        with open("%s/%s" % (fpath, "effect_electron.dat"), "w") as f:
            f.write("# kpath(angstrom^-1)    eigen_states(J)\n")
            f.write("# rsqures: %15.9f\n" % bestfitting[0, 2])
            np.savetxt(f, dat_to_write.reshape(-1, 2))
        # calculate effective mass
        print(kpath[int(bestfitting[0,0]):int(bestfitting[0,1])])
        print(cbm_band[int(bestfitting[0,0]):int(bestfitting[0,1])])
        effective_mass_calculator(kpath=kpath[int(bestfitting[0,0]):int(bestfitting[0,1])]/uc_ang2m, eigen=cbm_band[int(bestfitting[0,0]):int(bestfitting[0,1])]*uc_ev2j, fpath=fpath, fname="eff_mass_electron.png", plot=True)

    elif carrier == "hole":
        pass
    else:
        lexit("Carrier type: electron or hole")


def vasp_cbmvbm_index(xmlfile:str="vasprun.xml"):
    """
    1, get the index number of conduction bottom band and valence top band
    :param xmlfile: the input file
    :return two integers vbm_index, cbm_index
    """
    nelect = vaspxml_parser(select_attrib="parameters", xmlfile=xmlfile)["NELECT"]
    incar = vaspxml_parser(select_attrib="incar", xmlfile=xmlfile)
    if nelect % 2 == 1:
        nelect += 1
    if "LSORBIT" in incar:
        cbm_index = int(nelect)
        vbm_index = int(nelect-1)
    else:
        cbm_index = int(np.ceil(nelect/2.))
        vbm_index = int(np.ceil(nelect/2.-1))
    return vbm_index, cbm_index


def vasp_defomation_potential(fpath:str="./", poscar_file:str="POSCAR", direction:str="ac", step:float=0.005, num:int=5, init_analy:str="init"):
    """
    :param poscar:
    :param direction:
    :param step:
    :param num:
    :param init_analy: optional [init, analysis]
    """
    if init_analy == "init": 
        deformation = center_array(center=1, step=step, total_num=num)
        if "a" in direction:
            for _ in deformation:
                poscar = struct(fpath + "/" + poscar_file)
                poscar_dict = poscar.return_dict()
                poscar_dict["lattice"] = poscar_dict["lattice"]*np.array([[_], [1], [1]])
                poscar.update_struct_dict(poscar_dict)
                poscar.write_struct(fpath=fpath, fname="POSCAR_a_%5.3f" % _)
        if "c" in direction:
            for _ in deformation:
                poscar = struct(fpath + "/" + poscar_file)
                poscar_dict = poscar.return_dict()
                poscar_dict["lattice"] = poscar_dict["lattice"]*np.array([[1], [1], [_]])
                poscar.update_struct_dict(poscar_dict)
                poscar.write_struct(fpath=fpath, fname="POSCAR_c_%5.3f" % _)

    elif init_analy == "analysis":
        xml_file_list = subprocess.getoutput("ls %s/*.xml" % fpath).split()
        if "a" in direction:
            with open("%s/dp_a.dat" % fpath, "w") as f:
                f.write("%5s%15s%15s\n" % ("#dp(%)", "VBM", "CBM"))
                for _ in xml_file_list:
                    if "_a_" in _:
                        vbm_index, cbm_index = vasp_cbmvbm_index(xmlfile=_)
                        eigen, eigen2 = vaspxml_parser(select_attrib="eigen", xmlfile=_)
                        cbm = np.min(eigen[:, cbm_index, 0])
                        vbm = np.max(eigen[:, vbm_index, 0])
                        f.write("%5s%15.5f%15.5f\n" % (_[-9:-4], vbm, cbm))
            dp_dat_a = np.loadtxt("%s/dp_a.dat" % fpath)
            dp_vbm = np.polyfit(dp_dat_a[:, 0], dp_dat_a[:, 1], deg=1)
            dp_cbm = np.polyfit(dp_dat_a[:, 0], dp_dat_a[:, 2], deg=1)
            dp_vbm_x = np.linspace(dp_dat_a[0, 0], dp_dat_a[-1, 0], num=100)
            dp_cbm_x = np.linspace(dp_dat_a[0, 0], dp_dat_a[-1, 0], num=100)
            dp_vbm_y = dp_vbm[0]*dp_vbm_x + dp_vbm[1]
            dp_cbm_y = dp_cbm[0]*dp_cbm_x + dp_cbm[1]
            plt.plot(dp_vbm_x, dp_vbm_y, c="red", label="VBM, dp=%7.3f" % dp_vbm[0])
            plt.plot(dp_cbm_x, dp_cbm_y, c="blue", label="CBM, dp=%7.3f" % dp_cbm[0])
            plt.scatter(dp_dat_a[:, 0], dp_dat_a[:, 1], c="red", marker="x")
            plt.scatter(dp_dat_a[:, 0], dp_dat_a[:, 2], c="blue", marker="x")
            plt.legend(frameon=False)
            plt.xlabel("Strain in a-axis (%)")
            plt.ylabel("Energy (eV)")
            plt.savefig("%s/dp_a.png" % fpath)
            plt.close()
        if "b" in direction:
            with open("%s/dp_b.dat" % fpath, "w") as f:
                f.write("%5s%15s%15s\n" % ("#dp(%)", "VBM", "CBM"))
                for _ in xml_file_list:
                    if "_b_" in _:
                        vbm_index, cbm_index = vasp_cbmvbm_index(xmlfile=_)
                        eigen, eigen2 = vaspxml_parser(select_attrib="eigen", xmlfile=_)
                        cbm = np.min(eigen[:, cbm_index, 0])
                        vbm = np.max(eigen[:, vbm_index, 0])
                        f.write("%5s%15.5f%15.5f\n" % (_[-9:-4], vbm, cbm))
            dp_dat_b = np.loadtxt("%s/dp_b.dat" % fpath)
            dp_vbm = np.polyfit(dp_dat_b[:, 0], dp_dat_b[:, 1], deg=1)
            dp_cbm = np.polyfit(dp_dat_b[:, 0], dp_dat_b[:, 2], deg=1)
            dp_vbm_x = np.linspace(dp_dat_b[0, 0], dp_dat_b[-1, 0], num=100)
            dp_cbm_x = np.linspace(dp_dat_b[0, 0], dp_dat_b[-1, 0], num=100)
            dp_vbm_y = dp_vbm[0]*dp_vbm_x + dp_vbm[1]
            dp_cbm_y = dp_cbm[0]*dp_cbm_x + dp_cbm[1]
            plt.plot(dp_vbm_x, dp_vbm_y, c="red", label="VBM, dp=%7.3f" % dp_vbm[0])
            plt.plot(dp_cbm_x, dp_cbm_y, c="blue", label="CBM, dp=%7.3f" % dp_cbm[0])
            plt.scatter(dp_dat_b[:, 0], dp_dat_b[:, 1], c="red", marker="x")
            plt.scatter(dp_dat_b[:, 0], dp_dat_b[:, 2], c="blue", marker="x")
            plt.legend(frameon=False)
            plt.xlabel("Strain in b-axis (%)")
            plt.ylabel("Energy (eV)")
            plt.savefig("%s/dp_b.png" % fpath)
            plt.close()
        if "c" in direction:
            with open("%s/dp_c.dat" % fpath, "w") as f:
                f.write("%5s%15s%15s\n" % ("#dp(%)", "VBM", "CBM"))
                for _ in xml_file_list:
                    if "_c_" in _:
                        vbm_index, cbm_index = vasp_cbmvbm_index(xmlfile=_)
                        eigen, eigen2 = vaspxml_parser(select_attrib="eigen", xmlfile=_)
                        cbm = np.min(eigen[:, cbm_index, 0])
                        vbm = np.max(eigen[:, vbm_index, 0])
                        f.write("%5s%15.5f%15.5f\n" % (_[-9:-4], vbm, cbm))
            dp_dat_c = np.loadtxt("%s/dp_c.dat" % fpath)
            dp_vbm = np.polyfit(dp_dat_c[:, 0], dp_dat_c[:, 1], deg=1)
            dp_cbm = np.polyfit(dp_dat_c[:, 0], dp_dat_c[:, 2], deg=1)
            dp_vbm_x = np.linspace(dp_dat_c[0, 0], dp_dat_c[-1, 0], num=100)
            dp_cbm_x = np.linspace(dp_dat_c[0, 0], dp_dat_c[-1, 0], num=100)
            dp_vbm_y = dp_vbm[0]*dp_vbm_x + dp_vbm[1]
            dp_cbm_y = dp_cbm[0]*dp_cbm_x + dp_cbm[1]
            plt.plot(dp_vbm_x, dp_vbm_y, c="red", label="VBM, dp=%7.3f" % dp_vbm[0])
            plt.plot(dp_cbm_x, dp_cbm_y, c="blue", label="CBM, dp=%7.3f" % dp_cbm[0])
            plt.scatter(dp_dat_c[:, 0], dp_dat_c[:, 1], c="red", marker="x")
            plt.scatter(dp_dat_c[:, 0], dp_dat_c[:, 2], c="blue", marker="x")
            plt.legend(frameon=False)
            plt.xlabel("Strain in c-axis (%)")
            plt.ylabel("Energy (eV)")
            plt.savefig("%s/dp_c.png" % fpath)
            plt.close()

    else:
        lexit("Available choice of init_analy: init, analysis")



def add_selective_dyn(inp:str="POSCAR", out:str="POSCAR_init" ):
    """"""
    pass


def vaspxml_parser(select_attrib:str, xmlfile:str="vasprun.xml"):
    """
    
    :param select_arrib : optional [tdos, pdos, eigen, final_ene, paramters, incar]
    :return    final_ene: the total energy
    :               kist: 
    :               tdos:
    :               pdos:
    :              eigen:
    :         parameters: return a dictionary parameters_dict, available key: NELECT, NBAND, efermi
    :              incar:
    :            fermi_e: 
    :return final_ene         : the 
    :return klist             : array like kpoints list and the the accumulative k length [frac_reciprocal_x,y,z, abs_reciprocal_x,y,z, kpath_length]
    :return highsympoint      : 
    :return eigen             : return two eigenval arrays with a shape of (kpoint, [eigenvalues, occupation])
    :return paramters         : parameters_dict key  : NELECT
                                                     : NBANDS
    :return incar: incar_dict : 
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
        dos = calculation.find("dos")
        # fermi : dos[0]
        # total dos: dos[1] -> data: dos[1][0][5][0]
        # partial dos: dos[2] -> data: dos[2][0][5][0]
        fermi = dos[0].text
        dos_lines_tot = "#" + "%s" % fermi + "".join(dos[1][0][5][0].itertext())
        # column field: energy s py pz ...
        col_field = ''
        for _ in dos[2][0][3:-2]: 
            col_field += _.text
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

        is_ispin = False
        incar = vaspxml_parser(select_attrib="incar", xmlfile=xmlfile)
        if "ISPIN" in incar:
            if incar["ISPIN"] == 2:
                is_ispin = True

        eigen_spin1 = []
        eigen_spin2 = []
        for k in range(len(eigenvalues[0][5][0])):
            for e in range(len(eigenvalues[0][5][0][k])):
                eigen_spin1.append([float(i) for i in eigenvalues[0][5][0][k][e].text.split()])

        if is_ispin:
            return np.array(eigen_spin1).reshape(-1, len(eigenvalues[0][5][0][0]), 2), np.array(eigen_spin2(-1, len(eigenvalues[0][5][0][0]), 2))
        else:
            return np.array(eigen_spin1).reshape(-1, len(eigenvalues[0][5][0][0]), 2), eigen_spin2

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
            highsympoints_abs.append(highsympoints_abs[_-1]+np.sqrt(np.dot(cartesian_reciprocal[_]-cartesian_reciprocal[_-1], cartesian_reciprocal[_]-cartesian_reciprocal[_-1])))
        return np.array(highsympoints_abs)
    elif select_attrib == "parameters":
        parameters = vasproot.find("parameters")
        parameters = parameters.findall("separator")
        attributes = [_.attrib["name"] for _ in parameters]
        parameters_electronic = parameters[attributes.index("electronic")]
        attributes_electronic = [_.attrib["name"] for _ in parameters_electronic]

        calculation = vasproot.find("calculation")
        dos = calculation.find("dos")
        efermi = float(dos[0].text)

        parameters_dict = {}
        parameters_dict["NELECT"] = float(parameters_electronic[attributes_electronic.index("NELECT")].text)
        parameters_dict["NBANDS"] = float(parameters_electronic[attributes_electronic.index("NBANDS")].text)
        parameters_dict["efermi"] = efermi

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
        for _ in range(len(structure_end[0][1])):
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


def vasp_band_extractor_previous(eign:str='EIGENVAL', poscar:str="POSCAR", kpoint:str="KPOINTS"):
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
    efermi = vaspxml_parser(select_attrib="parameters", xmlfile=xmlfile)["efermi"]

    # write file
    with open(fpath+"/"+fname, "w") as f:
        f.write('# high symmetry points:  ')
        f.write(' %s\n' % ''.join("{:15.8f}".format(i) for i in kpoint_high_sym))
        f.write('# efermi: %.8f\n' % efermi)
        for _ in range(len(eigen1[0,:,0])):
            f.write('# %s \n' % str(_+1))
            np.savetxt(f, np.hstack((kpath.reshape(-1, 1), eigen1[:,_,:])), fmt="%12.8f" "%10.4f" "%10.4f")
            f.write('\n\n')
    
    # write the gnuplot script
    x_max = kpoint_high_sym[-1] * 1.001
    xtics = ''
    K_label = 1
    for i in kpoint_high_sym:
        xtics += '"K%d" %.8f, ' % (K_label, i)
        K_label += 1
    high_points = ''.join("{:15.8f}".format(i) for i in kpoint_high_sym[1:-1])

    with open(fpath+"/"+fname+".gnu", "w") as f:
        f.write(gnubandcar % (fname, x_max, xtics[:-2], high_points, efermi))


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
            potcar += open(str(potpath)+"/POTCAR_"+atom, "r").readlines()
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
    poscar.write_struct(fpath="./", fname="POSCAR_vacuum.vasp")


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


def vasp_gen_input(dft="scf", potpath="./", poscar="POSCAR", dryrun=False, wpath="./", wname:str="none", execode="srun --mpi=pmi2 vasp_std", params="gga=pbelda", **kwargs):
    func_help = """
    generate inputs for vasp
    --dft       : optional   opt  ->  
                :            scf  ->  
                :           band  ->  
                :            dos  ->  
                :     conv_encut  ->  
                :     conv_kmesh  ->  
                :     conv_sigma  ->  
                :          xanes  ->  
                :  ifc3-phono3py  ->  
                :           born  ->  
                :          elast  ->  
                :             dp  -> 
                :       aimd-nvt  -> 
                :       aimd-npt  ->
                :       aimd-nve  ->
                :
                :
                :
    --potpath   : absolute path to the PP library and renames them as POTCAR_Ti
    --poscar    : poscar
    --dryrun    : show this help
    --wpath     : the parent path containing the working diectory 
    --wname     : The name of working folder. If specified, generate working
                  direction with given name instead of default one.
    --param     : dft=any         ->  gga      = [pbelda, pbesol, hse06, hsesol, rev-vdW-DF2, 
                :                 ->             optB88, optPBE]
                :                 ->  oddeven  = [odd, even]
                :                 ->  kspacing = 0.15
                :                 ->  charged  = -1
                :                 ->  any tag in INCAR
                : dft=opt         ->  mulisif  = 2/6/3
                :                 ->  mulprec  = Low/Normal/Normal
                :                 ->  
                :                 ->  
                :                 ->  



    write       : create a folder containing all the necessary input files for VASP, a file 
                : of cal_details, 

                      : [dft=opt]           ->    
                      : [dft=scf]           ->    
                      : [dft=band]          ->    NBAND=200, 
                      : [dft=conv_encut]    ->    ENCUT=200-300-400-500, kmesh [0.05-0.1-0.2]
                      : [dft=conv_kmesh]    ->    kspacing [0.05-0.07-0.1-0.15-0.2-0.3]
                      : [dft=xanes]         ->    hole [1s,2s,2p,...]
                      : [dft=born]          ->    
                      : [dft=if3-phono3py]  ->    dim = ["2 2 2", ...]
                      : [dft=dp]            ->    direct=[a b c], step=0.05, range=0.1, e.g. --param direct=ac,step=0.05,range=0.1
                      : [dft=ml_heat]       -> 
    """
    if dryrun:
        print(func_help)
        exit()
    
    poscar = struct(poscar)
    val_electron = 0
    incar = {}
    
    params = parser_inputpara(params)
    gga = params["gga"]

    # =================================================================================
    # global setting
    # =================================================================================
    # gen directory
    if wpath == "./": 
        wpath = os.path.abspath("./")
    wdir = wpath+"/"+dft+"/"
    if wname != "none":
        wdir = wpath+"/"+wname+"/"
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

    # global incar
    incar.update(incar_glob)
    # update setting: "d" for d-elements, "f" for f-elements; set the LMAXMIX=2/4/6
    if max(poscar.return_dict()["atoms_index"]) > 57:
        incar["LMAXMIX"] = "6"
    elif max(poscar.return_dict()["atoms_index"]) > 20:
        incar["LMAXMIX"] = "4"
    else:
        incar["LMAXMIX"] = "2"

    # charged system with charged=-1 in parameters
    if "charged" in params:
        incar["NELECT"] = str(val_electron - int(params["charged"]))
    
    # write run.sh file and add x permission
    def write_runsh(runsh, cmd):
        with open(runsh, "a") as f:
            f.write(cmd)
        os.chmod(runsh, 0o775)

    # change gga functionals
    def ch_functional(incar, wdir=wdir):
        """change functional from parameters"""
        update_incar(incar, incar_functionals[gga], list(incar_functionals[gga].keys()))
        return incar
    
    # set GGA+U, need inputs from kwargs ["ggau"]
    if "ggau" in params:
        val_u = {}
        val_j = {}
        incar["# gga+u setting "] = " "
        incar["LDAU"] = ".TRUE."
        incar["LDAUTYPE"] = "2"
        incar["LDAUL"] = ""
        incar["LDAUU"] = ""
        incar["LDAUJ"] = ""
        try:
            val_u = parser_inputpara(kwargs["ggau"])
        except:
            pass
        try:
            val_j = parser_inputpara(kwargs["ggaj"])
        except:
            pass
        for atom in poscar.return_dict()["atoms_type"]:
            if atom in val_u:
                incar["LDAUL"] = incar["LDAUL"] + "  2"
                incar["LDAUU"] = incar["LDAUU"] + "  " + val_u[atom]
            else:
                incar["LDAUL"] = incar["LDAUL"] + " -1"
                incar["LDAUU"] = incar["LDAUU"] + "  0.00"
            if atom in val_j:
                incar["LDAUJ"] = incar["LDAUJ"] + "  " + val_j[atom]
            else:
                incar["LDAUJ"] = incar["LDAUJ"] + "  0.00"

        
    
    # =================================================================================
    # scf
    # =================================================================================
    def dft_scf(incar=incar):
        incar.update(incar_scf)
        update_incar(incar, params, incar_tag)
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
        cmd += "mv OUTCAR OUTCAR_%s\n" % dftgga
        cmd += "mv CHGCAR CHGCAR_%s\n" % dftgga
        cmd += "mv WAVECAR WAVECAR_%s\n" % dftgga
        cmd += "mv vasprun.xml vasprun_%s.xml\n" % dftgga
        cmd += "rm -f CHG POSCAR KPOINTS CONTCAR DOSCAR EIGENVAL INCAR REPORT vaspout.h5 XDATCAR\n"
        write_runsh(wdir+"/run.sh", cmd)
    
    # =================================================================================
    # opt
    # =================================================================================
    def dft_opt(incar=incar):
        incar.update(incar_opt)
        update_incar(incar, params, incar_tag)
        incar = ch_functional(incar)

        # enable WAVECAR
        incar["LWAVE"] = ".TRUE."
        incar["LCHARG"] = ".TRUE."

        # default ISIF = 3
        isif = {"ISIF": "3"}
        if "ISIF" in params:
            isif["ISIF"] = params["ISIF"]

        vasp_kpoints_gen(poscar.return_dict(), 
                         kspacing=float(params["kmesh"]) if "kmesh" in params else 0.3, 
                         kmesh=params["oddeven"] if "oddeven" in params else "odd", 
                         fpath=wdir, 
                         fname="KPOINTS_opt")

        # for multi-isif: need inputs from kwargs: ["params2", "params3"]
        if "mulisif" in params:
            mulstep = len(params["mulisif"]) + 1
            step = 1
            
            # write first inputs
            update_incar(incar, isif, ["ISIF"])
            write_incar(incar, 
                        fpath=wdir, 
                        fname="INCAR_%s_%s_isif%s_%d" \
                               %(dft, gga, incar["ISIF"], step))
            poscar.write_struct(fpath=wdir, 
                                fname="POSCAR_opt_init")
            dftgga = "%s_%s_isif%s" % (dft, gga, incar["ISIF"])
            cmd = "# opt isif=%s calculation\n" % incar["ISIF"]
            cmd += "cp POSCAR_opt_init POSCAR\n"
            cmd += "cp KPOINTS_opt KPOINTS\n"
            cmd += "for step in %d; do\n" % step
            cmd += "cp INCAR_%s_$step INCAR\n" % dftgga
            cmd += execode + "\n"
            cmd += "mv OUTCAR OUTCAR_%s_$step\n" % dftgga
            cmd += "mv vasprun.xml vasprun_%s_$step.xml\n" % dftgga
            cmd += "mv XDATCAR XDATCAR_%s_$step\n" % dftgga
            cmd += "mv OSZICAR OSZICAR_%s_$step\n" % dftgga
            cmd += "cp CONTCAR CONTCAR_%s_$step\n" % dftgga
            cmd += "mv CONTCAR POSCAR\n" 
            cmd += "done\n"
            write_runsh(wdir+"/run.sh", cmd)

            # write folowing steps
            for step in range(2, mulstep):
                params_tmp = parser_inputpara(kwargs["params%d" % step])
                update_incar(incar, params_tmp, incar_tag)

                write_incar(incar, 
                        fpath=wdir, 
                        fname="INCAR_%s_%s_isif%s_%d" \
                               %(dft, gga, incar["ISIF"], step))
            
                dftgga = "%s_%s_isif%s" % (dft, gga, incar["ISIF"])
                cmd = "# opt isif=%s calculation\n" % incar["ISIF"]
                cmd += "cp CONTCAR_%s_%d POSCAR\n" % (dftgga, step-1)
                cmd += "for step in %d; do\n" % step
                cmd += "cp INCAR_%s_$step INCAR\n" % dftgga
                cmd += execode + "\n"
                cmd += "mv OUTCAR OUTCAR_%s_$step\n" % dftgga
                cmd += "mv vasprun.xml vasprun_%s_$step.xml\n" % dftgga
                cmd += "mv XDATCAR XDATCAR_%s_$step\n" % dftgga
                cmd += "mv OSZICAR OSZICAR_%s_$step\n" % dftgga
                cmd += "cp CONTCAR CONTCAR_%s_$step\n" % dftgga
                cmd += "mv CONTCAR POSCAR\n" 
                cmd += "done\n"
                write_runsh(wdir+"/run.sh", cmd)







        else:
            update_incar(incar, isif, ["ISIF"])
            write_incar(incar, 
                        fpath=wdir, 
                        fname="INCAR_%s_%s_isif%s" %(dft, gga, isif["ISIF"]))
            poscar.write_struct(fpath=wdir, fname="POSCAR_opt_init")
            dftgga = dft+"_"+gga+"_isif"+isif["ISIF"]
            cmd = "# opt isif=%s calculation\n"
            cmd += "cp INCAR_%s INCAR\n" % dftgga
            cmd += "cp POSCAR_opt_init POSCAR\n"
            cmd += "cp KPOINTS_opt KPOINTS\n"
            cmd += "for step in 1; do\n"
            cmd += execode + "\n"
            cmd += "mv OUTCAR OUTCAR_%s_$step\n" % dftgga
            cmd += "mv vasprun.xml vasprun_%s_$step.xml\n" % dftgga
            cmd += "mv XDATCAR XDATCAR_%s_$step\n" % dftgga
            cmd += "mv OSZICAR OSZICAR_%s_$step\n" % dftgga
            cmd += "cp CONTCAR CONTCAR_%s_$step\n" % dftgga
            cmd += "mv CONTCAR POSCAR\n" 
            cmd += "done\n"
            cmd += "rm -f CHG KPOINTS DOSCAR EIGENVAL INCAR REPORT vaspout.h5\n"
            write_runsh(wdir+"/run.sh", cmd)
    
    # =================================================================================
    # convergence test: conv_encut conv_kmesh
    # =================================================================================
    def dft_conv_encut(incar=incar):
        """ generate input files for encut convergence test """
        incar.update(incar_scf)
        update_incar(incar, params, incar_tag)
        incar = ch_functional(incar)
        if "encut" in params:
            encut = params["encut"].split("-")
        else:
            encut = ["300", "350", "400", "450", "500", "550", "600", "650", "700", "750", "800"]
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
        cmd += "mv OUTCAR OUTCAR_%s_encut$encut\n" % (dftgga)
        cmd += "mv OSZICAR OSZICAR_%s_encut$encut\n" % (dftgga)
        cmd += "mv vasprun.xml vasprun_%s_encut$encut.xml\n" % (dftgga)
        cmd += "rm -f CHGCAR WAVECAR\n"
        cmd += "done\n"
        cmd += "rm -f CHG POSCAR KPOINTS CONTCAR DOSCAR EIGENVAL INCAR REPORT vaspout.h5 XDATCAR\n"
        write_runsh(wdir+"/run.sh", cmd)


    def dft_conv_kmesh(incar=incar):
        """ generate input files for kmesh convergence test """
        incar.update(incar_scf)
        update_incar(incar, params, incar_tag)
        incar = ch_functional(incar)

        if "encut" in params:
            kspacing = params["kspacing"].split("-")
        else:
            kspacing = [0.08, 0.1, 0.12, 0.14, 0.15, 0.2, 0.25, 0.3, 0.35]
        
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
        cmd += "mv vasprun.xml vasprun_%s_kmesh$kspacing.xml\n" % (dftgga)
        cmd += "done\n"
        write_runsh(wdir+"/run.sh", cmd)
    
    def dft_conv_vacuum(incar=incar): 
        """ generate input files for vacuum convergence test (2d matertials)"""
    
    def dft_conv_layer(incar=incar): 
        """ generate input files for atomic layers convergence test (2d matertials)"""

    # =================================================================================
    # dos, band, effective mass
    # =================================================================================
    
    def dft_band(incar=incar):
        """ """
        print("Make sure WAVECAR_scf and CHGCAR_scf are in the same directory, otherwise, the calculation will be perfomed with dense k-mesh and comsume more time.")
        print("Gnerate the k-path by yourself and named as KPOINTS_band")
        incar.update(incar_band)
        update_incar(incar, params, ["ENCUT", "PREC", "NBAND"])
        incar = ch_functional(incar)

        dftgga = dft+"_"+gga
        write_incar(incar, fpath=wdir, fname="INCAR_%s" % dftgga)
        cmd = "# band calculation\n"
        cmd += "cp INCAR_%s INCAR\n" % dftgga
        cmd += "cp KPOINTS_band KPOINTS\n"
        cmd += "cp WAVECAR_%s WAVECAR\n" % ("scf_"+gga)
        cmd += "cp CHGCAR_%s CHGCAR\n" % ("scf_"+gga)
        cmd += execode + "\n"
        cmd += "mv vasprun.xml vasprun_%s.xml\n" % dftgga
        write_runsh(wdir+"/run_band.sh", cmd)
    
    def dft_dos(incar=incar):
        """ """
        print("Make sure WAVECAR_scf and CHGCAR_scf are in the same directory, otherwise, the calculation will be perfomed with dense k-mesh and comsume more time.")
        print("Gnerate the k-path by yourself and named as KPOINTS_band")
        incar.update(incar_band)
        update_incar(incar, params, ["ENCUT", "PREC", "NBAND"])
        incar = ch_functional(incar)

        vasp_kpoints_gen(poscar.return_dict(), kspacing=float(params["kmesh"]) if "kmesh" in params else 0.1, kmesh=params["oddeven"] if "oddeven" in params else "odd", fpath=wdir, fname="KPOINTS_dos")

        dftgga = dft+"_"+gga
        write_incar(incar, fpath=wdir, fname="INCAR_%s" % dftgga)
        cmd = "# scf calculation\n"
        cmd += "cp INCAR_%s INCAR\n" % dftgga
        cmd += "cp KPOINTS_dos KPOINTS\n"
        cmd += "cp WAVECAR_%s WAVECAR\n" % ("scf_"+gga)
        cmd += "cp CHGCAR_%s CHGCAR\n" % ("scf_"+gga)
        cmd += execode + "\n"
        cmd += "cp vasprun.xml vasprun_%s.xml\n" % dftgga
        write_runsh(wdir+"/run_dos.sh", cmd)

    def dft_dp(incar=incar):
        incar.update(incar_scf)
        update_incar(incar, params, ["ENCUT", "PREC"])
        incar = ch_functional(incar)
        vasp_kpoints_gen(poscar.return_dict(), kspacing=float(params["kmesh"]) if "kmesh" in params else 0.15, kmesh=params["oddeven"] if "oddeven" in params else "odd", fpath=wdir, fname="KPOINTS_scf")

        # generate structure
        struct_dict_new = poscar.return_dict()
        struct_dict_lattice = poscar.return_dict()["lattice"]
        dp_step = np.linspace(1-float(params["range"]), 1+float(params["range"]), int(((1+float(params["range"]))-(1-float(params["range"])))/(float(params["step"]))+1))
        
        for step in dp_step:
            if "a" in params["direct"]:
                struct_dict_new["lattice"] = struct_dict_lattice*np.array([[step],[1],[1]])
                poscar.update_struct_dict(struct_dict_new)
                poscar.write_struct(fpath=wdir, fname="POSCAR_%s_%5.3f" % ("a", step))
            if "b" in params["direct"]:
                struct_dict_new["lattice"] = struct_dict_lattice*np.array([[1],[step],[1]])
                poscar.update_struct_dict(struct_dict_new)
                poscar.write_struct(fpath=wdir, fname="POSCAR_%s_%5.3f" % ("b", step))
            if "c" in params["direct"]:
                struct_dict_new["lattice"] = struct_dict_lattice*np.array([[1],[1],[step]])
                poscar.update_struct_dict(struct_dict_new)
                poscar.write_struct(fpath=wdir, fname="POSCAR_%s_%5.3f" % ("c", step))

        dftgga = dft+"_"+gga
        write_incar(incar, fpath=wdir, fname="INCAR_%s" % dftgga)
        cmd = "# deformation potential calculation\n"
        cmd += "for pos in POSCAR_*; do \n"
        cmd += "cp INCAR_%s INCAR\n" % dftgga
        cmd += "cp $pos POSCAR\n"
        cmd += "cp KPOINTS_scf KPOINTS\n"
        cmd += execode + "\n"
        cmd += "cp OUTCAR OUTCAR_$pos\n" 
        cmd += "cp vasprun.xml vasprun_$pos.xml\n" 
        cmd += "done\n"
        write_runsh(wdir+"/run.sh", cmd)
    
    # =================================================================================
    # machine learning part
    # =================================================================================
    def ml_heat(incar=incar):
        """ """

    # =================================================================================
    # AIMD NVT-NPT-NVE
    # =================================================================================
    def dft_nvt(incar=incar):
        incar.update(incar_nvt)
        update_incar(incar, params, incar_tag)
        incar = ch_functional(incar)

        vasp_kpoints_gen(poscar.return_dict(), kspacing=float(params["kmesh"]) if "kmesh" in params else 0.3, kmesh=params["oddeven"] if "oddeven" in params else "odd", fpath=wdir, fname="KPOINTS_md")

        dftgga = dft+"_"+gga
        write_incar(incar, fpath=wdir, fname="INCAR_%s" % dftgga)
        poscar.write_struct(fpath=wdir, fname="POSCAR_md_init")
        
        cmd = "# md nvt calculation\n"
        cmd += "cp INCAR_%s INCAR\n" % dftgga
        cmd += "cp POSCAR_md_init POSCAR\n"
        cmd += "cp KPOINTS_md KPOINTS\n"
        cmd += "for step in 1; do\n"
        cmd += execode + "\n"
        cmd += "mv OUTCAR OUTCAR_%s_$step\n" % dftgga
        cmd += "mv vasprun.xml vasprun_%s_$step.xml\n" % dftgga
        cmd += "mv XDATCAR XDATCAR_%s_$step\n" % dftgga
        cmd += "mv OSZICAR OSZICAR_%s_$step\n" % dftgga
        cmd += "mv REPORT REPORT_%s_$step\n" % dftgga
        cmd += "mv PCDAT PCDAT_%s_$step\n" % dftgga
        cmd += "cp CONTCAR CONTCAR_%s_$step\n" % dftgga
        cmd += "mv CONTCAR POSCAR\n" 
        cmd += "done\n"
        cmd += "rm -f CHG KPOINTS DOSCAR EIGENVAL INCAR vaspout.h5\n"
        write_runsh(wdir+"/run.sh", cmd)



    if dft == "scf":
        dft_scf()  
    elif dft == "opt":
        dft_opt()
    elif dft == "conv_encut":
        dft_conv_encut()
    elif dft == "conv_kmesh":
        dft_conv_kmesh()
    elif dft == "band":
        dft_band()
    elif dft == "dos":
        dft_dos()
    elif dft == "dp":
        dft_dp()
    elif dft == "md-nvt":
        dft_nvt()
    elif dft == "elecall":
        """ opt + scf + band + dos """
        vasp_gen_input(dft="opt", potpath=potpath, poscar=poscar, dryrun=False, wpath="./elecall", execode=execode, params=params)
        vasp_gen_input(dft="scf", potpath=potpath, poscar=poscar, dryrun=False, wpath="./elecall", execode=execode, params=params)
        vasp_gen_input(dft="band", potpath=potpath, poscar=poscar, dryrun=False, wpath="./elecall", execode=execode, params=params)
        vasp_gen_input(dft="dos", potpath=potpath, poscar=poscar, dryrun=False, wpath="./elecall", execode=execode, params=params)

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
    
    # =================================================================================
    # BORN effective charge
    # =================================================================================
    if dft == "born":
        """calculate BORN effective charge with DFPT method"""
        incar.update(incar_glob)
        incar.update(incar_prop["born"])
        update_incar(incar, params, incar_tag)
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
        cmd += "mv OUTCAR OUTCAR_%s\n" % dftgga
        cmd += "mv vasprun.xml vasprun_%s.xml\n" % dftgga
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
    :param struct1, struct2: string, the name of structures; or dictionary 
    """
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
        if type == "add": 
            vacuum = float(vacuum)
        elif type == "fit": 
            vacuum = float(vacuum) - lattice_para[2,2]
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

    vbm_cbm_kpoints_distance = np.dot(klist[cbm_kindex[0], :3] - klist[vbm_kindex[0], :3], klist[cbm_kindex[0], :3] - klist[vbm_kindex[0], :3])
    if vbm_cbm_kpoints_distance < 0.001:
        gaptype = "Direct"
    else:
        gaptype = "Indirect"

    return {
        "vbm_kindex": vbm_kindex,
        "vbm_kpoint": klist[vbm_kindex, :3],
        "vbm_ene": eigen1[vbm_kindex, vbm_index-1, 0],
        "cbm_kindex": cbm_kindex,
        "cbm_kpoint": klist[cbm_kindex, :3],
        "cbm_ene": eigen1[cbm_kindex, cbm_index-1, 0],
        "ene_gap": eigen1[cbm_kindex[0], cbm_index-1, 0] - eigen1[vbm_kindex[0], vbm_index-1, 0],
        "gap_type": gaptype
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
    

def gen_gnu_bandcar(xml:str, bandcar:str="BANDCAR", kpoints:str="KPOINTS_band"):
    """"""
    nelect = vaspxml_parser("parameters", xml)["NELECT"]
    efermi = vaspxml_parser("parameters", xml)["efermi"]
    incar = vaspxml_parser("incar", xml)
    if "LSORBIT" in incar:
        valence_index = int(nelect)
        conduction_index = valence_index + 1
    else:
        valence_index = int(nelect/2)
        conduction_index = valence_index + 1
    print(valence_index)
    print(conduction_index)

    with open(bandcar, "r") as f:
        bandcar_lines = f.readlines()
    with open(kpoints, "r") as f:
        kpoints_lines = f.readlines() 
    
    highpoints_value = [float(i) for i in bandcar_lines[1][1:].split()]
    print(highpoints_value)
    highpoints_label = []
    for i in range(4, len(kpoints_lines)):
        if len(kpoints_lines[i].split()) != 0:
            highpoints_label.append(kpoints_lines[i].split()[-1])
    highpoints_label = del_list_dupli_neighbor(highpoints_label)
    print(highpoints_label)

    xtics = ""
    for i in range(len(highpoints_value)):
        xtics += '"%s" %s, ' % (highpoints_label[i], str(highpoints_value[i]))
    
    midhigh = " ".join([str(i) for i in highpoints_value[1:-1]])

    # get band gap
    vbm_index = 0
    cbm_index = 0
    for i in range(len(bandcar_lines)):
        if "# "+str(valence_index) in bandcar_lines[i]:
            vbm_index = i + 1
        elif "# "+str(conduction_index) in bandcar_lines[i]:
            cbm_index = i + 1
    vbm_band = np.loadtxt(bandcar, skiprows=vbm_index, max_rows=cbm_index-vbm_index-3, usecols=[1])
    cbm_band = np.loadtxt(bandcar, skiprows=cbm_index, max_rows=cbm_index-vbm_index-3, usecols=[1])
    band_gap = np.min(cbm_band) - np.max(vbm_band)
    # print(band_gap)


    gnuscripts = gnubandcar % (bandcar, highpoints_value[-1]+0.0001, xtics, midhigh, bandcar, str(efermi))
    
    with open("plot.gnu", "w") as f:
        f.write(gnuscripts)
    
        