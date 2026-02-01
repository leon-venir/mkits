import numpy as np
import matplotlib.pyplot as plt
import spglib as spg
import os
import sys
import logging
import xml.etree.ElementTree as ET
from mkits import structure
from mkits import functions
from mkits import database
import pandas as pd



def vasp_gap_scissor(
        inp = "vasprun.xml",
        gap = 1.0, # in eV
        eperband = 2
):
    """
    DESCRIPTION:
    ------------
    Expand band gap to specific value, even for metal
    Use it carefully
    
    PARAMETERS:
    -----------
    inp:
        input
    gap:
    eperband:

    RETURN:
    -------
    """

    if "xml" in inp:
        # number of electrons
        nvalence = float(next(
            functions.find_lines_with_keyword(
            file_path=inp,
            keyword="NELECT",
            returnall=False
            )
        )[1][17:-5].replace(" ", "")) / eperband
        # get the number of band
        # if half occupied, treat it as a valence band
        nvalence = int(np.ceil(nvalence))

        # number of energy states
        #nkpoints = next(functions.find_lines_with_keyword(inp, "weights", False))[0] - \
        #           next(functions.find_lines_with_keyword(inp, "kpointlist", False))[0] - 2
        # get eigenvalue lineindex
        eigenindex = [int(i) for i, j in functions.find_lines_with_keyword(inp, '<set comment="kpoint')]
        
        # eigenstates
        # only support kpoints >= 2
        nstates = eigenindex[1] - eigenindex[0] - 2

        # scissors
        lines = open(inp, "r").readlines()
        for i in eigenindex:
            for j in range(i+nvalence, i+nstates):
                print(float(lines[j][10:20])+gap, lines[j][20:])
        #        lines[j] = "       <r>%10.4f%s" % (float(lines[j][10:20])+gap, lines[20:])

        # write 
        #open("%s_gap1" % inp, "w").writelines(lines)

    elif "EIGENVAL" in inp:
        pass
    elif "OUTCAR" in inp:
        pass
    else:
        print("Only xml, OUTCAR, EIGENVAL are supported.")


def contcar_statistics(
    wkdir="./",
    contcar="CONTCAR_opt"
):
    """
    
    """
    lattice6 = np.array([0, 0, 0, 0, 0, 0])
    func = os.listdir(wkdir)
    for _ in func:
        cell = structure.struct(wkdir+"/"+_+"/"+contcar)

        lattice6 = np.vstack((lattice6, cell.lattice6))
    lattice6 = pd.DataFrame(
        lattice6,
        columns=[r"$a$", r"$b$", r"$c$", r"$\alpha$", r"$\beta$", r"\gamma$"]
    )
    return lattice6



def extract_conv_test(
    wdir="./",
    plot=False
):
    """
    Description:
    ------------
    Extract data from convergence test calculation, 
    [input files are generated from func vasp_gen_input, ie encut, kmesh

    Parameters:
    -----------
    :param wdir: working directory
    """
    conv_name_list = []
    for _ in os.listdir(wdir):
        if "vasprun.xml_conv_" in _:
            conv_name_list.append(_)

    # check the type of convergence: vasprun_conv_[kmesh, encut]_[gga]
    if "vasprun.xml_conv_encut_" in conv_name_list[0]:
        
        ene = [
            vaspxml_parser(
                wdir+"/"+_, 
                "final_ene"
            )
            for _ in conv_name_list
        ]
            
        encut = [
            float((_.replace(".xml", "")).split("_encut_")[-1]) for _ in conv_name_list
            ]
        ene = np.array(ene)[np.argsort(encut)]
        encut = np.array(encut)[np.argsort(encut)]
        
        if plot:
            plt.plot(encut, ene, marker="x", lw=1.5)
            plt.xlabel("Kinetic Energy Cutoff (eV)")
            plt.ylabel("Energy (eV)")
        else:
            np.savetxt(
                wdir+"/encut-ene.dat", 
                (np.vstack((np.array(encut), np.array(ene))).T)[np.argsort(encut)[::-1], :], 
                header="encut(eV) ene(eV)"
            )
            with open(wdir+"/encut-ene.gnu", "w", newline="\n") as f:
                f.write(gnu2dline % ("encut-ene.png", 
                "Convergence test of kinetic energy cutoff", 
                "Kinetic energy cutoff(eV)", 
                "Total energy(eV)", 
                "encut-ene.dat"))
    
    if "vasprun.xml_conv_kmesh_" in conv_name_list[0]:
        ene = [
            vaspxml_parser(
                wdir+"/"+_,
                "final_ene"
            ) 
            for _ in conv_name_list
        ]
        kmesh = [
            float((_.replace(".xml", "")).split("_kmesh_")[-1]) for _ in conv_name_list
        ]
        ene = np.array(ene)[np.argsort(kmesh)]
        kmesh = np.array(kmesh)[np.argsort(kmesh)]

        if plot:
            plt.plot(kmesh, ene, marker="x", lw=1.5)
            plt.xlabel(r"K-spacing ($\AA^{-1}$)")
            plt.ylabel("Energy (eV)")
        else:
            np.savetxt(wdir+"/kmesh-ene.dat", (np.vstack((np.array(kmesh), np.array(ene))).T)[np.argsort(kmesh)[::-1], :], header="k-mesh ene(eV)")
            with open(wdir+"/kmesh-ene.gnu", "w", newline="\n") as f:
                f.write(gnu2dline % ("kmesh-ene.png", 
                                     "Convergence test of k-mesh", 
                                     "K-Spacing", 
                                     "Total energy(eV)", 
                                     "kmesh-ene.dat"))


def vasp_kpoints_gen(
        lattice9,
        kspacing=0.3, 
        kmesh="none", 
        fname="KPOINTS",
        kfix="-1 -1 -1"
):
    """
    Generate kpoints

    Parameters
    ----------
    poscar_dict: 
        dictionary, structure dict
    kspacing: 
        float
    kmesh: 
        string: odd, even, none(default); generate even or odd kmesh
    fix:
        Fix k-mesh in specific direction to a given number
        -1 -1 -1 means no fix
        -1 -1  1 means fix the k-mesh in z-axis direction to 1
    """
    
    a1 = np.sqrt(
        np.dot(lattice9[0,:], lattice9[0,:])
    )
    a2 = np.sqrt(
        np.dot(lattice9[1,:], lattice9[1,:])
    )
    a3 = np.sqrt(
        np.dot(lattice9[2,:], lattice9[2,:])
    )
    b1 = 2*np.pi/a1
    b2 = 2*np.pi/a2
    b3 = 2*np.pi/a3
    if kmesh=="none":
        oddeven = -1
    elif kmesh=="odd" or kmesh=="ODD":
        oddeven = 1
    elif kmesh=="even" or kmesh=="EVEN":
        oddeven = 0
    oddeven = -1

    n1 = int(
        np.max(
            np.array([1, functions.round_even_odd(
                b1/(kspacing), 
                oddeven
            )])
        )
    )
    n2 = int(
        np.max(
            np.array([1, functions.round_even_odd(
                b2/(kspacing), 
                oddeven
            )])
        )
    )
    n3 = int(
        np.max(
            np.array([1, functions.round_even_odd(
                b3/(kspacing), 
                oddeven
            )])
        )
    )

    # check fix
    kfix = [int(i) for i in kfix.split()]
    if kfix[0] != -1:
        n1 = kfix[0]
    if kfix[1] != -1:
        n2 = kfix[1]
    if kfix[2] != -1:
        n3 = kfix[2]
    
    with open(fname, "w", newline="\n") as f:
        f.write("mesh auto\n0\nG\n%s %s %s\n0 0 0" % (str(n1), str(n2), str(n3)))


def vasp_extract_band(xmlfile,
                      out="./band.dat",
                      thresholdfactor=5):
    """
    """
    klist = vaspxml_parser(xmlfile, "klist")
    eigen1, eigen2 = vaspxml_parser(xmlfile, "eigenvalues")
    recp_lattice = vaspxml_parser(xmlfile, "reciprocal_parameter")

    lines = functions.extractband(
        recp_lattice,
        klist,
        eigen1,
        thresholdfactor
    )

    if out == "none":
        return lines
    else:
        with open(out, "w") as f:
            f.writelines(lines)


def vasp_dos_extractor(
        xmlfile="vasprun.xml",
        shift_fermi=True,
        datatype="tot",
        fpath="./",
        **kwargs):
    """
    datatype:
        Options: [tot, ldos, pdos]
    """

    dos_tot, dos_partial, ion_num, orbit = vaspxml_parser(
        xmlfile=xmlfile,
        select_attrib="dos"
    )

    energy = dos_tot[:, 0].reshape(-1,1)
    if shift_fermi:
        energy -= vaspxml_parser(
            xmlfile=xmlfile,
            select_attrib="parameters"
        )["efermi"]
    
    atomic_types = vaspxml_parser(
            xmlfile=xmlfile,
            select_attrib="atoms"
        )

    try:
        ispin = vaspxml_parser(
            xmlfile=xmlfile,
            select_attrib="incar"
            )["ISPIN"]
    except:
        ispin = "1"
    
    if datatype == "tot":
        with open(
            fpath+"/dos_tot.dat",
            "w",
            newline="\n"
        ) as f:
            f.write("# tot DOS \n")
            f.write("# energy up up-integ (dn dn-integ if any) \n")
            np.savetxt(
                fname=f,
                X=np.hstack((
                    energy,
                    dos_tot[:, 1:],
                )),
                fmt="%15.6f"
            )
    elif datatype == "pdos":
        shape = dos_partial.shape
        # ene s p d f px py pz
        atomic_type = list(set(atomic_types))
        atomic_types = np.array(atomic_types)
        
        orbit = np.array([_[0] for _ in orbit])
        pdos_s = dos_partial[:, :, :, np.where(orbit=="s")[0]]
        pdos_p = dos_partial[:, :, :, np.where(orbit=="p")[0]]
        pdos_d = dos_partial[:, :, :, np.where(orbit=="d")[0]]

        pdos_p_sum = np.sum(pdos_p, axis=3).reshape(2, shape[1], shape[2], 1)
        pdos_d_sum = np.sum(pdos_d, axis=3).reshape(2, shape[1], shape[2], 1)

        if True: # "2" in ispin:
            with open(fpath+"/pdos_up.dat", "w", newline="\n") as f:
                head = []
                for i in ["s", "p", "d"]:
                    for j in atomic_type:
                        head.append(i+"-"+j)
                f.write("# pdos-up %s\n"  % " ".join(head))
                pdos_up = energy
                for _ in atomic_type:
                    pdos_up = np.hstack((
                        pdos_up, 
                        np.sum(pdos_s[0, np.where(atomic_types==_)[0], :], 
                               axis=0)
                    ))
                    pdos_up = np.hstack((
                        pdos_up, 
                        np.sum(pdos_p_sum[0, np.where(atomic_types==_)[0], :], 
                               axis=0)
                    ))
                    pdos_up = np.hstack((
                        pdos_up, 
                        np.sum(pdos_d_sum[0, np.where(atomic_types==_)[0], :], 
                               axis=0)
                    ))
                np.savetxt(f, pdos_up, "%15.6f")
            with open(fpath+"/pdos_dn.dat", "w", newline="\n") as f:
                head = []
                for i in ["s", "p", "d"]:
                    for j in atomic_type:
                        head.append(i+"-"+j)
                f.write("# pdos-dn %s\n"  % " ".join(head))
                pdos_up = energy
                for _ in atomic_type:
                    pdos_up = np.hstack((
                        pdos_up, 
                        np.sum(pdos_s[1, np.where(atomic_types==_)[0], :], 
                               axis=0)
                    ))
                    pdos_up = np.hstack((
                        pdos_up, 
                        np.sum(pdos_p_sum[1, np.where(atomic_types==_)[0], :], 
                               axis=0)
                    ))
                    pdos_up = np.hstack((
                        pdos_up, 
                        np.sum(pdos_d_sum[1, np.where(atomic_types==_)[0], :], 
                               axis=0)
                    ))
                np.savetxt(f, pdos_up, "%15.6f")
        
    
    elif datatype == "ldos":
        shape = dos_partial.shape
        # ene s p d f px py pz
        atomic_type = list(set(atomic_types))
        atomic_types = np.array(atomic_types)
        
        orbit = np.array([_[0] for _ in orbit])
        pdos_s = dos_partial[:, :, :, np.where(orbit=="s")[0]]
        pdos_p = dos_partial[:, :, :, np.where(orbit=="p")[0]]
        pdos_d = dos_partial[:, :, :, np.where(orbit=="d")[0]]
        pdos_p_sum = np.sum(pdos_p, axis=3).reshape(2, shape[1], shape[2], 1)
        pdos_d_sum = np.sum(pdos_d, axis=3).reshape(2, shape[1], shape[2], 1)

        select_idx = [int(_) for _ in kwargs["index"].split(",")]

        if "2" in ispin:
            with open(fpath+"/ldos_up.dat", "w", newline="\n") as f:
                f.write("# ldos-up s, p, d %s\n"  % " ".join(kwargs["index"].split(",")))
                pdos_up = energy
                pdos_up = np.hstack((
                    pdos_up,
                    np.sum(pdos_s[0, select_idx, :],
                        axis=0)
                ))
                pdos_up = np.hstack((
                    pdos_up,
                    np.sum(pdos_p_sum[0, select_idx, :],
                        axis=0)
                ))
                pdos_up = np.hstack((
                    pdos_up,
                    np.sum(pdos_d_sum[0, select_idx, :],
                        axis=0)
                ))             
                np.savetxt(f, pdos_up, "%15.6f")
            with open(fpath+"/ldos_dn.dat", "w", newline="\n") as f:
                f.write("# ldos-up s, p, d %s\n"  % " ".join(kwargs["index"].split(",")))
                pdos_up = energy
                pdos_up = np.hstack((
                    pdos_up,
                    np.sum(pdos_s[1, select_idx, :],
                        axis=0)
                ))
                pdos_up = np.hstack((
                    pdos_up,
                    np.sum(pdos_p_sum[1, select_idx, :],
                        axis=0)
                ))
                pdos_up = np.hstack((
                    pdos_up,
                    np.sum(pdos_d_sum[1, select_idx, :],
                        axis=0)
                ))             
                np.savetxt(f, pdos_up, "%15.6f")


def band_gap_overlapping(
        xml_file
):
    """
    
    """
    incar = read_vaspxml(
        xml_file,
        "incar"
    )
    eperband = 2
    if "LSORBIT" in incar.keys():
        eperband = 1
    electron = read_vaspxml(
        xml_file,
        "nelectron"
    )


def vasp_get_vbm_cbm(
        xml_file
):
    """
    
    """
    incar = read_vaspxml(
        xml_file,
        "incar"
    )
    eperband = 2
    if "LSORBIT" in incar.keys():
        eperband = 1
    electron = read_vaspxml(
        xml_file,
        "nelectron"
    )
    eigenvalues = read_vaspxml(
        xml_file,
        "eigenvalues"
    )
    if electron % 2 == 0:
        vbm_idx = int(electron / eperband)
        vbm = np.max(eigenvalues[:, vbm_idx-1, 0])
        cbm = np.min(eigenvalues[:, vbm_idx, 0])
        return vbm, cbm
    else:
        vbm_idx = int(electron / eperband)
        vbm = np.max(eigenvalues[:, vbm_idx, 0])
        cbm = np.min(eigenvalues[:, vbm_idx, 0])
        return vbm, cbm
    

def vasp_get_vbm_near(
        xml_file,
        idx = 1
):
    """
    
    """
    incar = read_vaspxml(
        xml_file,
        "incar"
    )
    eperband = 2
    if "LSORBIT" in incar.keys():
        eperband = 1
    electron = read_vaspxml(
        xml_file,
        "nelectron"
    )
    eigenvalues = read_vaspxml(
        xml_file,
        "eigenvalues"
    )
    if electron % 2 == 0:
        vbm_idx = int(electron / eperband)
        vbm = np.sort(eigenvalues[:, vbm_idx-1, 0])[-idx-1]
        return vbm
    else:
        vbm_idx = int(electron / eperband)
        vbm = np.sort(eigenvalues[:, vbm_idx, 0])[-idx-1]
        return vbm


def vasp_get_band_gap_overlapping(
        xml_file,
        overlapping=False
):
    """
    
    """
    
    vbm, cbm = vasp_get_vbm_cbm(
        xml_file
    )
    gap_overlapping = cbm - vbm
    # treat overlapping as 0
    if overlapping and gap_overlapping < 0:
        return gap_overlapping
    elif gap_overlapping < 0:
        return 0
    elif gap_overlapping >= 0:
        return gap_overlapping


def read_vaspxml(
    xml_file,
    select_entry
):
    """
    
    select_entry:
        lattice_basis -> 3x3 numpy array
        atom_position -> nx3 numpy array
        rec_basis
        nelectron
        efermi
        eigenvalues   -> numpy array (2, NKPOINTS, NBANDS, 2)
                         [Spin_Index, Kpoint_Index, Band_Index, (Energy, Occ)]
                         If ISPIN=1, Spin 2 is a copy of Spin 1.
        totaldos      -> numpy array [Energy, DOS_Spin1, DOS_Spin2]
                         (If ISPIN=1, DOS_Spin2 is a copy of DOS_Spin1)
        partialdos    -> numpy array [Energy, Ion1_Spin1, Ion1_Spin2, Ion2_Spin1...]
                         (If ISPIN=1, Spin2 columns are copies of Spin1)
        atominfo      -> list of strings ['Fe', 'Fe', 'O', ...] (Length = N_ions)

    """
    tree = ET.parse(xml_file)
    root = tree.getroot()

    if select_entry == "rec_basis":
        rec_basis = []
        for varray in root.findall(
            ".//structure[@name='primitive_cell']//crystal//varray[@name='rec_basis']"
        ):
            for v in varray.findall('v'):
                rec_basis.append(
                    [float(x) for x in v.text.split()]
                )
        return np.array(rec_basis)
    
    elif select_entry == "incar":
        incar = {}
        for key in root.findall(
            ".//incar"
        ):
            for i in key.findall('i'):
                incar[i.attrib.get("name")] = i.text
        return incar
    
    # ==========================================
    # Added: symbol (Extract atom types list)
    # ==========================================
    elif select_entry == "symbol":
        symbols = []
        # Find the 'atoms' array in 'atominfo'
        # Structure: <atominfo><array name="atoms"><set><rc><c>Element</c><c>Type</c>...
        atoms_ref = root.find(".//atominfo/array[@name='atoms']/set")
        
        if atoms_ref is not None:
            for rc in atoms_ref.findall("rc"):
                # Each row (rc) has columns (c). The first column is the element symbol.
                # Remove whitespace just in case (e.g., "Fe " -> "Fe")
                c_tags = rc.findall("c")
                if c_tags:
                    element_str = c_tags[0].text.strip()
                    symbols.append(element_str)
        return symbols

    elif select_entry == "nelectron":
        return float(
            root.findall(
                ".//parameters//separator[@name='electronic']//i[@name='NELECT']"
            )[0].text.split()[0]
        )

    elif select_entry == "efermi":
        # Usually found in <dos><i name="efermi">
        efermi_tag = root.find(".//dos/i[@name='efermi']")
        if efermi_tag is not None:
            return float(efermi_tag.text.strip())
        else:
            # Fallback
            tag = root.find(".//calculation/dos/i[@name='efermi']")
            return float(tag.text.strip()) if tag is not None else 0.0

    elif select_entry == "lattice_basis":
        lattice_basis = []
        for varray in root.findall(
            ".//structure[@name='primitive_cell']//crystal//varray[@name='basis']"
        ):
            for v in varray.findall('v'):
                lattice_basis.append([float(x) for x in v.text.split()])
        return np.array(lattice_basis)

    elif select_entry =="kpointlist":
        kpointlist = []
        for varray in root.findall(
            ".//kpoints//varray[@name='kpointlist']"
        ):
            for v in varray.findall('v'):
                kpointlist.append([float(x) for x in v.text.split()])
        return np.array(kpointlist)

    elif select_entry =="atom_position":
        atom_position = []
        for varray in root.findall(
            ".//structure[@name='primitive_cell']//varray[@name='positions']"
        ):
            for v in varray.findall('v'):
                atom_position.append([float(x) for x in v.text.split()])
        return np.array(atom_position)
    
    # ==========================================
    # Modified: eigenvalues (Merged & Auto-Copy)
    # ==========================================
    elif select_entry == "eigenvalues":
        # Get number of kpoints to reshape correctly
        kpointlist = read_vaspxml(xml_file, "kpointlist")
        kpoint_num = len(kpointlist) if kpointlist is not None else 0
        
        # 1. Parse Spin 1
        eigenvalues_1 = []
        spin1_node = root.find(".//calculation/eigenvalues/array/set/set[@comment='spin 1']")
        
        if spin1_node is not None:
            for setkpoint in spin1_node.findall(".//set"):
                for r in setkpoint.findall('r'):
                    # Each r.text is like "energy occupation"
                    eigenvalues_1.append([float(x) for x in r.text.split()])
            # Reshape: (NKPOINTS, NBANDS, 2)
            data_spin1 = np.array(eigenvalues_1).reshape(kpoint_num, -1, 2)
        else:
            return np.array([]) # Return empty if no eigenvalues found

        # 2. Parse Spin 2
        eigenvalues_2 = []
        spin2_node = root.find(".//calculation/eigenvalues/array/set/set[@comment='spin 2']")
        
        if spin2_node is not None:
            for setkpoint in spin2_node.findall(".//set"):
                for r in setkpoint.findall('r'):
                    eigenvalues_2.append([float(x) for x in r.text.split()])
            # Reshape: (NKPOINTS, NBANDS, 2)
            data_spin2 = np.array(eigenvalues_2).reshape(kpoint_num, -1, 2)
        else:
            # COPY Spin 1 if Spin 2 is missing
            data_spin2 = data_spin1.copy()

        # 3. Stack into a single array
        # Shape: (2, NKPOINTS, NBANDS, 2)
        # Index 0 -> Spin 1, Index 1 -> Spin 2
        return np.stack((data_spin1, data_spin2), axis=0)

    # ==========================================
    # Modified: totaldos (Force 3 columns)
    # ==========================================
    elif select_entry == "totaldos":
        # Find the total DOS block
        dos_total = root.find(".//calculation/dos/total/array/set")
        if dos_total is None: return np.array([])

        # Spin 1
        spin1_node = dos_total.find("./set[@comment='spin 1']")
        data_spin1 = []
        if spin1_node is not None:
            for r in spin1_node.findall("r"):
                data_spin1.append([float(x) for x in r.text.split()])
        data_spin1 = np.array(data_spin1) # [Energy, DOS, IntDOS]
        
        # Spin 2
        spin2_node = dos_total.find("./set[@comment='spin 2']")
        if spin2_node is not None:
            data_spin2 = []
            for r in spin2_node.findall("r"):
                data_spin2.append([float(x) for x in r.text.split()])
            data_spin2 = np.array(data_spin2)
            dos_down = data_spin2[:, 1]
        else:
            # COPY Spin 1 as Spin 2 if not enabled
            dos_down = data_spin1[:, 1]
            
        # Return: [Energy, Up, Down]
        return np.column_stack((data_spin1[:, 0], data_spin1[:, 1], dos_down))

    # ==========================================
    # Modified: partialdos (Force shape consistency)
    # ==========================================
    elif select_entry == "partialdos":
        # Find the partial DOS block
        dos_partial = root.find(".//calculation/dos/partial/array/set")
        if dos_partial is None: return np.array([])
        
        all_data = []
        energy_col = None
        
        # Iterate over ions
        for ion_set in dos_partial.findall("set"):
            # ion_comment = ion_set.get("comment") # e.g. "ion 1"
            
            # Spin 1
            spin1_node = ion_set.find("./set[@comment='spin 1']")
            s1_data = []
            if spin1_node is not None:
                for r in spin1_node.findall("r"):
                    s1_data.append([float(x) for x in r.text.split()])
            s1_data = np.array(s1_data)
            
            # Capture energy from the first ion
            if energy_col is None:
                energy_col = s1_data[:, 0]
            
            # Append Spin 1 orbitals (s, p, d...)
            all_data.append(s1_data[:, 1:])
            
            # Spin 2
            spin2_node = ion_set.find("./set[@comment='spin 2']")
            if spin2_node is not None:
                s2_data = []
                for r in spin2_node.findall("r"):
                    s2_data.append([float(x) for x in r.text.split()])
                s2_data = np.array(s2_data)
                # Append Spin 2 orbitals
                all_data.append(s2_data[:, 1:])
            else:
                # COPY Spin 1 orbitals as Spin 2
                all_data.append(s1_data[:, 1:])

        # Flatten: [Energy, Ion1_S1, Ion1_S2, Ion2_S1, Ion2_S2 ...]
        if all_data:
            flat_orbitals = np.hstack(all_data)
            return np.column_stack((energy_col, flat_orbitals))
        else:
            return np.array([])
    
    # ==========================================
    # Added: orbital (Get labels for partialdos columns)
    # ==========================================
    elif select_entry == "orbital":
        # 1. Get Atom Symbols
        symbols = read_vaspxml(xml_file, "symbol")
        if not symbols: return []

        # 2. Get Orbital Names from XML definition
        # Location: calculation -> dos -> partial -> array -> field
        partial_array = root.find(".//calculation/dos/partial/array")
        if partial_array is None: return []

        raw_orbitals = []
        for field in partial_array.findall("field"):
            name = field.text.strip()
            if "energy" not in name.lower(): # Exclude energy column
                raw_orbitals.append(name)
        
        # 3. Construct the list matching partialdos data structure
        # Structure: For each Ion -> Spin1(Orbs) -> Spin2(Orbs)
        labels = []
        for i, sym in enumerate(symbols):
            # Spin 1
            for orb in raw_orbitals:
                labels.append(f"{sym}{i+1}-{orb}-spin1")
            
            # Spin 2 (Always present in partialdos output due to copy logic)
            for orb in raw_orbitals:
                labels.append(f"{sym}{i+1}-{orb}-spin2")
        
        return labels
    
    else:
        print("Available Entry: lattice_basis, atom_position, rec_basis, nelectron, efermi, eigenvalues, eigenvalues2, totaldos, partialdos")
        return None


def band_extractor(
        xmlfile="vasprun.xml", 
        shift_fermi=True, 
        write_data=True
):
    """
    Extract band structure data (K-path and Eigenvalues).

    Returns:
    --------
    np.array
        Columns: [Kpath, Spin1_Band1, Spin1_Band2..., Spin2_Band1, Spin2_Band2...]
        Units: 1/Angstrom for Kpath, eV for Energy.
    """
    
    # 1. Read Data
    lattice_basis = read_vaspxml(xmlfile, "lattice_basis") # Real space basis (Angstrom)
    kpointlist = read_vaspxml(xmlfile, "kpointlist")       # Fractional coordinates
    eigenvalues = read_vaspxml(xmlfile, "eigenvalues")     # Shape: (2, NK, NB, 2)
    efermi = read_vaspxml(xmlfile, "efermi")
    
    if kpointlist is None or len(kpointlist) == 0:
        print("Error: No kpoints found.")
        return None

    # 2. Calculate Reciprocal Lattice (Physics convention with 2*pi)
    # B = 2*pi * (A^-1).T
    # This ensures units are 1/Angstrom
    rec_metric = 2 * np.pi * np.linalg.inv(lattice_basis).T
    
    # 3. Calculate K-path (Cumulative Distance)
    kpath = [0.0]
    high_sym_points = [0.0] # Always include start point
    
    # Convert all kpoints to Cartesian coordinates first
    # k_cart = k_frac @ rec_metric
    k_cart = kpointlist @ rec_metric
    
    for i in range(1, len(k_cart)):
        # Distance between current and previous k-point
        dist = np.linalg.norm(k_cart[i] - k_cart[i-1])
        kpath.append(kpath[-1] + dist)
        
        # Detect High Symmetry Points (Kinks)
        # Check if the direction of the path changes
        if i < len(k_cart) - 1:
            vec_prev = k_cart[i] - k_cart[i-1]
            vec_next = k_cart[i+1] - k_cart[i]
            
            # Normalize vectors to check angle
            norm_prev = np.linalg.norm(vec_prev)
            norm_next = np.linalg.norm(vec_next)
            
            if norm_prev > 1e-5 and norm_next > 1e-5:
                # Cosine similarity
                cos_angle = np.dot(vec_prev, vec_next) / (norm_prev * norm_next)
                # If angle deviates significantly from 0 degrees (cos < 0.999), it's a kink
                if cos_angle < 0.999:
                    high_sym_points.append(kpath[-1])
    
    # Add the final point as a high symmetry point
    high_sym_points.append(kpath[-1])
    
    kpath = np.array(kpath)

    # 4. Process Eigenvalues
    # eigenvalues shape: (Spin, Kpoints, Bands, 2) -> Last dim is [Energy, Occ]
    # We only need Energy (index 0)
    
    # Extract Spin 1 Energies
    # Shape: (NK, NB)
    bands_spin1 = eigenvalues[0, :, :, 0]
    
    # Extract Spin 2 Energies
    # read_vaspxml guarantees shape consistency (duplicates spin1 if ispin=1)
    bands_spin2 = eigenvalues[1, :, :, 0]
    
    # Shift Fermi Level
    if shift_fermi:
        bands_spin1 -= efermi
        bands_spin2 -= efermi

    # 5. Construct Final Array
    # Columns: Kpath, S1_B1, S1_B2..., S2_B1, S2_B2...
    final_data = np.column_stack((kpath, bands_spin1, bands_spin2))
    
    # 6. Write Data
    if write_data:
        nbands = bands_spin1.shape[1]
        
        # Format Comment Lines
        # Line 1: High Symmetry Points locations
        high_sym_str = " ".join([f"{x:.5f}" for x in high_sym_points])
        comment_line1 = f"# High Symmetry Points (Kpath): {high_sym_str}"
        
        # Line 2: Column Headers
        headers = ["Kpath"] + \
                  [f"Spin1_Band{i+1}" for i in range(nbands)] + \
                  [f"Spin2_Band{i+1}" for i in range(nbands)]
        comment_line2 = "# " + " ".join([f"{h:>12s}" for h in headers])
        
        filename = "band.dat"
        with open(filename, "w") as f:
            f.write(comment_line1 + "\n")
            f.write(comment_line2 + "\n")
            np.savetxt(f, final_data, fmt="%15.8f")
            
        print(f"Band structure saved to {filename}")
        print(f"High symmetry points found at: {high_sym_str}")

    return final_data


def dos_extractor(
        xmlfile="vasprun.xml", 
        shift_fermi=True, 
        write_data=True,
        write_to="dos.dat",
        ldos=True
):
    """
    Extract DOS data from vasprun.xml.
    
    Parameters:
    -----------
    ldos : bool
        If True: Output DOS for every single ion (e.g., Fe1_s, Fe2_s...).
        If False: Sum DOS by element type (e.g., Fe_s, O_s...).
    
    Output columns structure (if ldos=False):
    Energy, 
    Spin1_Total, Spin1_Elem1_Total, Spin1_Elem1_s..., Spin1_Elem2_Total...
    """
    
    # 1. Read Data
    efermi = read_vaspxml(xmlfile, "efermi")
    total_dos = read_vaspxml(xmlfile, "totaldos")     
    partial_dos = read_vaspxml(xmlfile, "partialdos") 
    symbols = read_vaspxml(xmlfile, "symbol")         
    raw_labels = read_vaspxml(xmlfile, "orbital")     
    
    if total_dos is None or len(total_dos) == 0:
        print("Error: No DOS data found.")
        return None

    # 2. Parse Orbital Structure
    orb_groups = {
        's': ['s'],
        'p': ['p', 'px', 'py', 'pz'],
        'd': ['d', 'dxy', 'dyz', 'dz2', 'dxz', 'x2-y2'],
        'f': ['f', 'f-3', 'f-2', 'f-1', 'f0', 'f1', 'f2', 'f3']
    }
    
    # Identify which orbitals are present in the raw data
    n_ions = len(symbols)
    if n_ions > 0:
        n_cols_per_spin_ion = len(raw_labels) // (2 * n_ions)
        sample_labels = raw_labels[:n_cols_per_spin_ion]
        raw_orb_names = [l.split('-')[1] for l in sample_labels]
    else:
        raw_orb_names = []
        n_cols_per_spin_ion = 0

    # 3. Define Output Columns Order (e.g., s, p, px, py...)
    output_orb_order = []
    
    if 's' in raw_orb_names: output_orb_order.append('s')
        
    p_sub = [o for o in raw_orb_names if o in orb_groups['p']]
    if p_sub:
        if 'p' not in p_sub and len(p_sub) > 1: output_orb_order.append('p') 
        p_sorted = sorted(p_sub, key=lambda x: ['p', 'px', 'py', 'pz'].index(x) if x in ['p', 'px', 'py', 'pz'] else 99)
        output_orb_order.extend(p_sorted)

    d_sub = [o for o in raw_orb_names if o in orb_groups['d']]
    if d_sub:
        if 'd' not in d_sub and len(d_sub) > 1: output_orb_order.append('d')
        output_orb_order.extend([o for o in raw_orb_names if o in d_sub])

    f_sub = [o for o in raw_orb_names if o in orb_groups['f']]
    if f_sub:
        if 'f' not in f_sub and len(f_sub) > 1: output_orb_order.append('f')
        output_orb_order.extend([o for o in raw_orb_names if o in f_sub])

    # 4. Construct Data Array
    energy = total_dos[:, 0]
    if shift_fermi:
        energy -= efermi
        
    final_data_columns = [energy]
    final_headers = ["Energy"]
    
    # --- Inner Helper Functions ---
    def get_raw_col(ion_idx, spin_idx, orb_name):
        block_len = n_cols_per_spin_ion
        start_idx = 1 + ion_idx * (2 * block_len) + spin_idx * block_len
        if orb_name in raw_orb_names:
            rel_idx = raw_orb_names.index(orb_name)
            return partial_dos[:, start_idx + rel_idx]
        else:
            return np.zeros_like(energy)

    def get_derived_col(ion_idx, spin_idx, out_orb_name):
        if out_orb_name in raw_orb_names:
            return get_raw_col(ion_idx, spin_idx, out_orb_name)
        
        target_group = None
        for g_name, g_members in orb_groups.items():
            if out_orb_name == g_name:
                target_group = g_members
                break
        
        if target_group:
            sum_data = np.zeros_like(energy)
            found_any = False
            for member in target_group:
                if member in raw_orb_names and member != out_orb_name:
                    sum_data += get_raw_col(ion_idx, spin_idx, member)
                    found_any = True
            return sum_data if found_any else np.zeros_like(energy)
        return np.zeros_like(energy)

    # --- Build Columns ---
    spin_names = ["spin1", "spin2"]
    
    for s_idx, s_name in enumerate(spin_names):
        # 1. System Total DOS
        tdos_col = total_dos[:, s_idx + 1]
        final_data_columns.append(tdos_col)
        final_headers.append(f"{s_name}_Total")
        
        # ---------------------------------------------------------
        # Branch 1: LDOS = True (Per Ion)
        # ---------------------------------------------------------
        if ldos:
            for i_idx, sym in enumerate(symbols):
                prefix = f"{s_name}_{sym}{i_idx+1}"
                
                # Ion Total
                ion_total = np.zeros_like(energy)
                for raw_o in raw_orb_names:
                    ion_total += get_raw_col(i_idx, s_idx, raw_o)
                
                final_data_columns.append(ion_total)
                final_headers.append(f"{prefix}_Total")
                
                # Ion Orbitals
                for orb in output_orb_order:
                    data_col = get_derived_col(i_idx, s_idx, orb)
                    final_data_columns.append(data_col)
                    final_headers.append(f"{prefix}_{orb}")

        # ---------------------------------------------------------
        # Branch 2: LDOS = False (Per Element) [NEW]
        # ---------------------------------------------------------
        else:
            # Get unique elements preserving order (e.g. Fe, O)
            unique_elements = sorted(list(set(symbols)), key=symbols.index)
            
            for elem in unique_elements:
                prefix = f"{s_name}_{elem}"
                
                # Identify indices of this element
                # indices = [0, 1] for Fe, [2, 3, 4] for O
                ion_indices = [i for i, s in enumerate(symbols) if s == elem]
                
                # Element Total (Sum of all orbitals for all ions of this element)
                elem_total = np.zeros_like(energy)
                for i_idx in ion_indices:
                    for raw_o in raw_orb_names:
                        elem_total += get_raw_col(i_idx, s_idx, raw_o)
                
                final_data_columns.append(elem_total)
                final_headers.append(f"{prefix}_Total")
                
                # Element Orbitals (Sum specific orbital across all ions of this element)
                for orb in output_orb_order:
                    elem_orb_sum = np.zeros_like(energy)
                    for i_idx in ion_indices:
                        elem_orb_sum += get_derived_col(i_idx, s_idx, orb)
                    
                    final_data_columns.append(elem_orb_sum)
                    final_headers.append(f"{prefix}_{orb}")

    # Combine
    final_array = np.column_stack(final_data_columns)
    
    # 5. Write Data
    if write_data:
        indices_line = "# " + " ".join([f"{i:>15d}" for i in range(len(final_headers))])
        headers_line = "# " + " ".join([f"{h:>15s}" for h in final_headers])
        
        with open(write_to, "w") as f:
            f.write(indices_line + "\n")
            f.write(headers_line + "\n")
            np.savetxt(f, final_array, fmt="%16.8e")
            
    return final_array


def vaspxml_parser(xmlfile,
                   select_attrib):
    """
    XML parser

    PARAMETERS:
    -----------
    xmlfile:
        The input xml file.
    select_attrib: str
        Options: [
            incar, cell_parameter, reciprocal_parameter, parameters,
            final_ene, eigenvalues, klist, 
        ]
    """

    try:
        xml = ET.parse(xmlfile)
        xmlroot = xml.getroot()
        incar = xmlroot.find("incar")
        atominfo = xmlroot.find("atominfo")
        structure = xmlroot.find("structure")
        
        # [修改] 使用 findall 找到所有 calculation，并取最后一个 [-1]
        # 注意：如果文件不完整（正在计算中），这也可能报错，但这通常是我们想要的行为（读取最新一步）
        calculations = xmlroot.findall("calculation")
        if calculations:
            calculation = calculations[-1]
        else:
            print("Error: No calculation blocks found.")
            return None
            
        parameters = xmlroot.find("parameters")      
    except Exception as e:
        print("Error, cannot parse the xmlfile:", e)
        return None
    
    # incar
    incartag = {}
    for tag in incar:
        incartag[tag.attrib["name"]] = tag.text
    
    # cell parameters
    basis = []
    rec_basis = []
    for _ in range(len(structure[0][0])):
        basis.append(
            [float(i) for i in structure[0][0][_].text.split()]
        )
    for _ in range(len(structure[0][2])):
        rec_basis.append(
            [float(i) * 2 * np.pi 
             for i in structure[0][2][_].text.split()]
        )
    
    # klist
    kpoints = xmlroot.find("kpoints")
    klist = []
    kpoints = kpoints.findall("varray")
    attributes = [_.attrib["name"] for _ in kpoints]
    kpoints = kpoints[attributes.index("kpointlist")]
    for _ in range(len(kpoints)):
        klist.append([float(i) for i in kpoints[_].text.split()])
    klist = np.array(klist)
    
    if select_attrib == "incar":    
        return incartag
    
    elif select_attrib == "cell_parameter":
        return basis
    
    elif select_attrib == "reciprocal_parameter":
        return np.array(rec_basis)
    
    elif select_attrib == "parameters":
        parameters = parameters.findall("separator")
        attributes = [_.attrib["name"] for _ in parameters]
        parameters_electronic = parameters[
            attributes.index("electronic")
        ]
        attributes_electronic = [
            _.attrib["name"] for _ in parameters_electronic
        ]

        dos = calculation.find("dos")
        efermi = float(dos[0].text)

        parameters_dict = {}
        parameters_dict["NELECT"] = float(
            parameters_electronic[
                attributes_electronic.index("NELECT")
            ].text
        )
        parameters_dict["NBANDS"] = float(
            parameters_electronic[
                attributes_electronic.index("NBANDS")
            ].text
        )
        parameters_dict["efermi"] = efermi
        return parameters_dict
  
    elif select_attrib == "final_ene":
        # 获取 calculation 块中直接子节点的 energy (即该离子步的收敛能量)
        energy_block = calculation.find("energy")
        
        # 在 energy 块中查找名为 "e_0_energy" 的项 (Sigma -> 0)
        # 如果找不到，再尝试找 e_fr_energy (Free energy) 或回退到 [-1]
        e_0 = energy_block.find("./i[@name='e_0_energy']")
        if e_0 is not None:
            return float(e_0.text)
        else:
            return float(energy_block[-1].text)
    
    elif select_attrib == "eigenvalues":
        eigenvalues = calculation.find("eigenvalues")
        eigen_spin1 = []
        eigen_spin2 = []
        
        # 
        if "ISPIN" in incartag and int(incartag["ISPIN"]) == 2:
            pass
        else:
            for k in range(len(eigenvalues[0][5][0])):
                for e in range(len(eigenvalues[0][5][0][k])):
                    eigen_spin1.append(
                        [float(i) 
                         for i in eigenvalues[0][5][0][k][e].text.split()]
                    )
        eigen_spin1 = np.array(eigen_spin1).reshape(
            (len(klist), -1)
        )
        eigen_spin2 = np.array(eigen_spin1).reshape(
            (len(klist), -1)
        )

        return eigen_spin1, eigen_spin1

    elif select_attrib == "klist":   
        return np.array(klist)
    
    elif select_attrib == "atoms":
        ion_num = int(atominfo[0].text)
        atoms = []
        for _ in range(ion_num):
            atoms.append(atominfo[2][3][:][_][0].text)
        return atoms
    
    elif select_attrib == "dos":
        dos = calculation.find("dos")
        # fermi : dos[0]
        # total dos: dos[1] -> data: dos[1][0][5][0]
        # partial dos: dos[2] -> data: dos[2][0][5][0]
        dos_lines = "".join(dos[1][0][5][0].itertext())
        dos_lines = np.array([_.split() 
                                  for _ in dos_lines.split("\n")[1:-1]])
        dos_lines_tot = dos_lines.astype(np.float64)

        # ISPIN =2
        if "ISPIN" in incartag.keys():
            if "2" in incartag["ISPIN"]:
                dos_lines = "".join(dos[1][0][5][1].itertext())
                dos_lines = np.array([_.split()
                                      for _ in dos_lines.split("\n")[1:-1]])
                dos_lines_tot = np.hstack((
                    dos_lines_tot,
                    dos_lines.astype(np.float64)[:, 1:]
                ))

        # partial
        dos_lines_partial = np.array([])
        col_field = ''
        
        if "LORBIT" in incartag.keys():
            if "11" in incartag["LORBIT"]:
                # ion num
                ion_num = int(atominfo[0].text)
                # column field: energy s py pz ...
                for _ in dos[2][0][3:]: 
                    col_field += "%s " % _.text
                col_field = col_field.split()
                
                dos_lines = [
                    "".join(dos[2][0][-1][_][0].itertext()) 
                    for _ in range(ion_num)
                ]

                for i in dos_lines:
                    dos_lines_partial = np.append(
                        dos_lines_partial, 
                        np.array([_.split() 
                                  for _ in i.split("\n")[1:-1]])
                    )

                dos_lines_partial = dos_lines_partial.astype(np.float64)
                dos_lines_partial = dos_lines_partial.reshape((
                    ion_num, -1, len(col_field)
                ))

                # ispin = 2
                if "ISPIN" in incartag.keys():
                    if "2" in incartag["ISPIN"]:
                        dos_lines = [
                            "".join(dos[2][0][-1][_][1].itertext()) 
                            for _ in range(ion_num)
                        ]
                        for i in dos_lines:
                            dos_lines_partial = np.append(
                                dos_lines_partial, 
                                np.array([_.split() 
                                        for _ in i.split("\n")[1:-1]])
                            )
                        dos_lines_partial = dos_lines_partial.astype(np.float64)
                        dos_lines_partial = dos_lines_partial.reshape((
                            2, ion_num, -1, len(col_field)
                ))
                # ispin = 1 return up = dn
                else:
                    dos_lines_partial = np.array([dos_lines_partial, dos_lines_partial])

        return dos_lines_tot, dos_lines_partial, ion_num, col_field

    elif select_attrib == "none":
        kvector = np.array([0])


def initvaspsh(dftlabel, execode, step=0, iswave=False):
    """
    Write a bash scripts to prepare the vasp files.
    Parameters
    ----------
    dftlabel:
        The label of POSCAR

    Return
    ------
    String 
    """
    initstring = ""
    
    if dftlabel == "opt":
        if step == 0:
            initstring = "  cp POSCAR_opt_init POSCAR\n"
            initstring += "  cp INCAR_%s INCAR\n" % dftlabel
            initstring += "  cp KPOINTS_%s KPOINTS\n" % dftlabel
        elif step == 1:
            initstring = "  cp POSCAR_opt_init POSCAR\n"
            initstring += "  cp INCAR_%s_%d INCAR\n" % (dftlabel, step)
            initstring += "  cp KPOINTS_%s_%d KPOINTS\n" % (dftlabel, step)
        else:
            initstring = "  cp CONTCAR_%s_%d POSCAR\n" % (dftlabel, step-1)
            initstring += "  cp INCAR_%s_%d INCAR\n" % (dftlabel, step)
            initstring += "  cp KPOINTS_%s_%d KPOINTS\n" % (dftlabel, step)

    elif dftlabel in ["scf", "freq", "nvt"]:
        if step == 0:
            initstring = "  cp POSCAR_init POSCAR\n"
            initstring += "  cp INCAR_%s INCAR\n" % dftlabel
            initstring += "  cp KPOINTS_%s KPOINTS\n" % dftlabel
        elif step == 1:
            initstring = "  cp POSCAR_init POSCAR\n"
            initstring += "  cp INCAR_%s_%d INCAR\n" % (dftlabel, step)
            initstring += "  cp KPOINTS_%s_%d KPOINTS\n" % (dftlabel, step)
        else:
            initstring = "  cp CONTCAR_%s_%d POSCAR\n" % (dftlabel, step-1)
            initstring += "  cp INCAR_%s_%d INCAR\n" % (dftlabel, step)
            initstring += "  cp KPOINTS_%s_%d KPOINTS\n" % (dftlabel, step)
        if iswave:
            initstring += "  cp WAVECAR_%s_%d WAVECAR\n" % (dftlabel, step-1)
            initstring += "  cp CHGCAR_%s_%d CHGCAR\n" % (dftlabel, step-1)
    
    elif dftlabel == "dos" or dftlabel == "band":
        initstring = "  cp POSCAR_init POSCAR\n"
        initstring += "  cp INCAR_%s INCAR\n" % dftlabel
        initstring += "  cp KPOINTS_%s KPOINTS\n" % dftlabel
        initstring += "  cp WAVECAR_%s WAVECAR\n" % (dftlabel)
        initstring += "  cp CHGCAR_%s CHGCAR\n" % (dftlabel)
    
    initstring += "\n"
    initstring += "  " + execode + " >& vasp.out 2>&1\n\n"
    
    return initstring


def savevaspsh(savelist, savelabel):
    """
    Write a bash scripts to save the vasp files with provided labels.
    Parameters:
    -----------

    Return:
    -------
    String 
    """
    
    savelist2string = ""
    breakline = 0
    for i in savelist:
        if breakline > 3:
            savelist2string += " \\\n          "
            breakline = 0
        savelist2string += "'" + i + "' "
        breakline += 1
        
    
    savestring = '''  ffout=( %s )
  for ifout in ${ffout[@]}; do
    mv $ifout "$ifout"_%s
  done
    ''' % (savelist2string, savelabel)
    
    return savestring


def vasp_single_atom_ground(
        atom = "C",
        kinetic = 500,
):
    """
    Description:
    ------------
    

    Parameters:
    -----------
    
    """
    pass


def vasp_gen_input(
        dft="scf", 
        potpath="./", 
        poscar="POSCAR", 
        dryrun=False, 
        metal=True, 
        mag=False, 
        soc=False,
        wpath="./", 
        wname:str="none", 
        execode="srun --mpi=pmi2 vasp_std", 
        params="gga=pbelda", 
        **kwargs
):
    """
    
    """

    func_help = """

"""
    if dryrun:
        print(func_help)
        functions.lexit("Show vasp_gen help.")

    #==========================================================================
    # global setting
    # =========================================================================
    val_electron = 0
    incar = {}
    params = functions.parser_inputpara(params)

    gga = "pbesol"
    if "gga" in params.keys():
        gga = params["gga"]

    # get structure info
    poscar = structure.struct(poscar)
    atomic_indexs = poscar.position[1:, 0].tolist()
    atomic_index = np.array(
        list(
            int(_) for _ in set(atomic_indexs)
        )
    )
    atomic_index = atomic_index[np.argsort(atomic_index)].tolist()
    atomic_num = [
        atomic_indexs.count(_) for _ in atomic_index
    ]
    atomic_type = [
        database.atom_data[_][1] for _ in atomic_index
    ]

    # gen directory
    if wpath == "./": 
        wpath = os.path.abspath("./")
    wkdir = wpath+"/"+dft+"/"
    if wname != "none":
        wkdir = wpath+"/"+wname+"/"
    if not os.path.exists(wpath):
        os.mkdir(wpath)
    if not os.path.exists(wkdir):
        os.mkdir(wkdir)

    # gen POTCAR and get electron of valence
    if os.path.exists(wkdir+"/POTCAR"):
        os.remove(wkdir+"/POTCAR")
    potcar_lines = []
    if potpath == "recommend":
        potpath_recommend = ""
    else:
        for _ in range(len(atomic_type)):
            with open("%s/POTCAR_%s" %(potpath, atomic_type[_])) as f:
                lines = f.readlines()
            val_electron += float(lines[1][:-1])*atomic_num[_]
            potcar_lines += lines
        with open(wkdir+"/POTCAR", "w", newline="\n") as f:
            f.writelines(potcar_lines)
    # charged system with charged=-1 in parameters
    if "charged" in params:
        val_electron += - int(params["charged"])
        incar["NELECT"] = str(val_electron)

    # global incar
    incar.update(database.incar_glob)
    # update setting: "d" for d-elements, "f" for f-elements; set the LMAXMIX=2/4/6
    if 71 >= max(atomic_index) >= 57 or \
       103 >= max(atomic_index) >= 89 :
        incar["LMAXMIX"] = "6"
    elif 30 >= max(atomic_index) >= 21 or \
         48 >= max(atomic_index) >= 39 or \
         80 >= max(atomic_index) >= 72 or \
         112 >= max(atomic_index) >= 104 :
        incar["LMAXMIX"] = "4"
    else:
        incar["LMAXMIX"] = "2"
    
    # chage functionals
    incar.update(database.incar_functionals[gga])
    
    # set GGA+U, need inputs from kwargs "ggau" and "ggaj"
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
            val_u = functions.parser_inputpara(kwargs["ggau"])
        except:
            pass
        try:
            val_j = functions.parser_inputpara(kwargs["ggaj"])
        except:
            pass
        for atom in atomic_type:
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
                
    def add_U(incar, kwargs, atomic_type):
        val_u = {}
        val_j = {}
        incar["# gga+u setting "] = " "
        incar["LDAU"] = ".TRUE."
        incar["LDAUTYPE"] = "2"
        incar["LDAUL"] = ""
        incar["LDAUU"] = ""
        incar["LDAUJ"] = ""
        try:
            val_u = functions.parser_inputpara(kwargs["ggau"])
        except:
            pass
        try:
            val_j = functions.parser_inputpara(kwargs["ggaj"])
        except:
            pass
        for atom in atomic_type:
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
        return incar
    
    # copy vdw kernel to working folder
    if params["gga"] in ["rev-vdW-DF2", "optB88", "optPBE"]:
        try:
            with open(kwargs["vdw_kernel_path"], "rb") as f:
                vdwlines = f.readlines()
            with open(wkdir + "/vdw_kernel.bindat", "wb") as f:
                f.writelines(vdwlines)
        except:
            print("Cannot find the vdw kernel, remember to copy it to \n" +
                  "the working folder. Or you can enable the path to the\n" +
                  "file with vdw_kernel_path=/path/to/vdw_kernel\n")

    # nbands, specify the NBANDS based on the valence electrons
    try:
        nbands = float(params["nbands"])
        if nbands > 4 or nbands <=0:
            print("""
The number of bands typically out of reasonable range.
The reasonable range is 0 ~ 4.
""")
        incar["NBANDS"] = str(int(val_electron * nbands *1.5))
        params["NBANDS"] = incar["NBANDS"]
        #print(val_electron)
    except:
        pass

    # dynrange or not
    dyn = False
    if "dynrange" in kwargs.keys():
        dyn = True
        selec_dyn = kwargs["dynrange"]
        selec_dyn = functions.parser_inputpara(selec_dyn)
        dynrange = {"xmin": -1e8, "ymin": -1e8, "zmin": -1e8,
                    "xmax": 1e8,  "ymax": 1e8,  "zmax": 1e8}
        if "fix" in selec_dyn.keys():
            #print("fix")
            _fix_atom = selec_dyn["fix"]
            del(selec_dyn["fix"])
        else:
            _fix_atom = "none"
        if "move" in selec_dyn.keys():
            #print("move")
            _move_atom = selec_dyn["move"]
            del(selec_dyn["move"])
        else:
            _move_atom = "none"
        for key in selec_dyn.keys():
            dynrange[key] = float(selec_dyn[key])
        #print(dynrange)
        poscar.add_dyn(
            xmin=dynrange["xmin"],
            xmax=dynrange["xmax"],
            ymin=dynrange["ymin"],
            ymax=dynrange["ymax"],
            zmin=dynrange["zmin"],
            zmax=dynrange["zmax"],
            fix=_fix_atom,
            move=_move_atom
        )
        
    
    # WARNING of very long lattice vectors (>50 A) 
    # decrease AMIN to a smaller value (e.g. 0.01)
    if any ([poscar.lattice6[0] > 50,     # lattice a
             poscar.lattice6[1] > 50,     # lattice b
             poscar.lattice6[2] > 50,]):  # lattice c
        incar["AMIN"] = "0.01"
    
    # save file list
    save_output_list = [
        "OUTCAR", "OSZICAR", "REPORT", "IBZKPT", 
        "CONTCAR", "EIGENVAL", "vasprun.xml", 
        "CHGCAR", "WAVECAR", "DOSCAR", "vasp.out"
    ]
    def update_save_list(save_output_list, dft, incar):
        if gga == "opt":
            save_output_list.remove("WAVECAR")
            save_output_list.remove("CHGCAR")
        if "LVHAR" in incar:
            save_output_list.append("LOCPOT")
            save_output_list.append("POT")
        if "LORBIT" in incar:
            save_output_list.append("PROCAR")
        return save_output_list

    # kpoints parameters
    def update_kpoints_para(params):
        kspacing = float(
            params["kspacing"]
        ) if "kspacing" in params else 0.25
        kmesh=params[
            "oddeven"
        ] if "oddeven" in params else "odd"
        kfix=params[
            "kfix"
        ] if "kfix" in params else "-1 -1 -1"
        return kspacing, kmesh, kfix
    
    # update incar
    def update_incar(incar, params):
        """ add params to incar, if the key is duplicated then updates it """
        for key in params.keys():
            if key in database.incar_tag:
                incar.update(
                    {
                        key: params[key]
                    }
                )
        return incar
    
    # write incar
    def write_incar(incar, fname="INCAR"):
        """ write INCAR """
        with open(fname, "w", newline="\n") as f:
            for item in incar.keys():
                f.write(item+"="+incar[item]+"\n")
    
    # change functional from parameters
    incar = update_incar(
        incar,
        database.incar_functionals[gga]
    )

    # metal or not
    if metal:
        incar.update(database.incar_metal)
    else:
        incar.update(database.incar_semi)
    
    # set mag， soc need inputs from kwargs "magmom"
    if mag and not soc:
        incar["ISPIN"] = "2"
        magtag = ""
        if "MAGMOM" in params:
            magtag = params["MAGMOM"]
        else:
            for i in range(len(atomic_index)):
                magtag += "%d*%.2f " % (int(atomic_num[i]), 
                                        database.atom_data[int(atomic_index[i])][4])
        incar["MAGMOM"] = magtag
    elif not mag and soc:
        incar["LSORBIT"] = ".TRUE."
    elif mag and soc:
        pass
    else:
        pass

    # set SOC, need inputs from kwargs "ggau" and "ggaj"
    """
    metal=True, 
        mag=False, 
        soc=False,
    """
    
    # =========================================================================
    # scf
    # =========================================================================
    def dft_scf(incar=incar):
        """ write """
        poscar.write_struct(
            fpath=wkdir, 
            fname="POSCAR_init",
            calculator="poscar",
            dyn=dyn
        )
        incar.update(database.incar_scf)

        # for multi-scf: need inputs from kwargs: ["params2", "params3"]
        if "mulscf" in params:
            mulstep = int(params["mulscf"]) + 1
            step = 1

            # incar
            incar1 = update_incar(incar, params)
            if "ggau" in params:
                incar1 = add_U(incar1, kwargs, atomic_type)
            write_incar(
                incar=incar1,
                fname="%s/INCAR_%s_%d" % (wkdir, dft, step)
            )

            # kpoints
            kspacing, kmesh, kfix = update_kpoints_para(params)
            vasp_kpoints_gen(
                poscar.lattice9,
                kspacing,
                kmesh,
                wkdir + "KPOINTS_scf_%d" % step, 
                kfix
            )

            # run
            cmd = "# scf step %d\n" % step
            cmd += initvaspsh(dft, execode, 1, False)
            cmd += savevaspsh(
                update_save_list(
                    save_output_list,
                    dft,
                    incar
                ),
                dft+"_%d" % step
            )
            functions.write_runsh(
                wkdir+"/run_scf_%d.sh" % step,
                cmd
            )

            # write following steps
            for step in range(2, mulstep):
                params2 = functions.parser_inputpara(
                    kwargs["params%d" % step]
                )
                incar2 = update_incar(incar, params2)
                write_incar(
                    incar=incar2,
                    fname="%s/INCAR_%s_%d" % (wkdir, dft, step)
                )

                kspacing, kmesh, kfix = update_kpoints_para(params2)
                vasp_kpoints_gen(
                    poscar.lattice9,
                    kspacing,
                    kmesh,
                    wkdir + "KPOINTS_scf_%d" % step, 
                    kfix
                )
                cmd = "# scf step %d\n" % step
                cmd += initvaspsh(dft, execode, step, True)
                cmd += savevaspsh(
                    update_save_list(
                        save_output_list,
                        dft,
                        incar
                    ),
                    dft+"_%d" % step
                )
                functions.write_runsh(
                    wkdir+"/run_scf_%d.sh" % step,
                    cmd
                )
        else:
            incar = update_incar(incar, params)
            
            write_incar(
                incar=incar,
                fname="%s/INCAR_%s" % (wkdir, dft)
            )
            kspacing, kmesh, kfix = update_kpoints_para(params)
            vasp_kpoints_gen(
                poscar.lattice9,
                kspacing,
                kmesh,
                wkdir+"/KPOINTS_scf", 
                kfix
            )
            cmd = "# scf\n"
            cmd += initvaspsh(dft, execode, 0, False)
            cmd += savevaspsh(
                update_save_list(
                    save_output_list,
                    dft,
                    incar
                ),
                dft
            )
            functions.write_runsh(
                wkdir+"/run_scf.sh",
                cmd
            )
    
    # =========================================================================
    # opt
    # =========================================================================
    def dft_opt(incar=incar):
        # enable WAVECAR by default
        incar["LWAVE"] = ".TRUE."
        incar["LCHARG"] = ".TRUE."

        poscar.write_struct(
            fpath=wkdir, 
            fname="POSCAR_opt_init",
            calculator="poscar",
            dyn=dyn
        )
        incar.update(database.incar_opt)
        incar = update_incar(incar, params)


        # for multiisif=243
        # need inputs from kwargs: ["params2", "params3"]
        if "mulisif" in params:
            mulstep = len(params["mulisif"]) + 1
            step = 1

            incar1 = update_incar(incar, params)
            write_incar(
                incar=incar1,
                fname="%s/INCAR_%s_%d" % (wkdir, dft, step)
            )

            kspacing, kmesh, kfix = update_kpoints_para(params)
            vasp_kpoints_gen(
                poscar.lattice9,
                kspacing,
                kmesh,
                wkdir+"/KPOINTS_opt_%d" % step, 
                kfix
            )
            cmd = "# opt step %d\n" % step
            cmd += initvaspsh(dft, execode, 1, True)
            cmd += savevaspsh(
                update_save_list(
                    save_output_list,
                    dft,
                    incar
                ),
                dft+"_%d" % step
            )
            functions.write_runsh(
                wkdir+"/run_opt_%d.sh" % step,
                cmd
            )

            # write following steps
            for step in range(2, mulstep):
                params2 = functions.parser_inputpara(
                    kwargs["params%d" % step]
                )

                incar2 = update_incar(incar, params2)
                write_incar(
                    incar=incar2,
                    fname="%s/INCAR_%s_%d" % (wkdir, dft, step)
                )

                kspacing, kmesh, kfix = update_kpoints_para(params2)
                vasp_kpoints_gen(
                    poscar.lattice9,
                    kspacing,
                    kmesh,
                    wkdir + "KPOINTS_opt_%d" % step, 
                    kfix
                )
                cmd = "# opt step %d\n" % step
                cmd += initvaspsh(dft, execode, 1, True)
                cmd += savevaspsh(
                    update_save_list(
                        save_output_list,
                        dft,
                        incar
                    ),
                    dft+"_%d" % step
                )
                functions.write_runsh(
                    wkdir+"/run_opt_%d.sh" % step,
                    cmd
                )
        else:
            incar = update_incar(incar, params)
            write_incar(
                incar=incar,
                fname="%s/INCAR_%s" % (wkdir, dft)
            )
            kspacing, kmesh, kfix = update_kpoints_para(params)
            vasp_kpoints_gen(
                poscar.lattice9,
                kspacing,
                kmesh,
                wkdir+"/KPOINTS_opt", 
                kfix
            )
            cmd = "# opt\n"
            cmd += initvaspsh(dft, execode, 0, True)
            cmd += savevaspsh(
                update_save_list(
                    save_output_list,
                    dft,
                    incar
                ),
                dft
            )
            functions.write_runsh(
                wkdir+"/run_opt.sh",
                cmd
            )
    
    # =========================================================================
    # convergence test: conv_encut conv_kmesh
    # =========================================================================
    def dft_conv_encut(incar=incar):
        """generate input files for encut convergence test"""
        poscar.write_struct(
            fpath=wkdir, 
            fname="POSCAR_init",
            calculator="poscar",
            dyn=dyn
        )

        incar.update(database.incar_scf)
        incar = update_incar(incar, params)

        if "encut" in params:
            encut = params["encut"].split("-")
        else:
            encut = ["300", "350", "400", "450", "500", "550", "600", "650", 
                     "700", "750", "800"]
        incar["ENCUT"] = "xxx"

        write_incar(
            incar=incar,
            fname="%s/INCAR_%s" % (wkdir, dft)
        )
        kspacing, kmesh, kfix = update_kpoints_para(params)
        vasp_kpoints_gen(
            poscar.lattice9,
            kspacing,
            kmesh,
            wkdir+"/KPOINTS_scf", 
            kfix
        )

        cmd = "# encut convergence test\n"
        cmd += "  cp POSCAR_init POSCAR\n"
        cmd += "  cp KPOINTS_scf KPOINTS\n"
        cmd += "  for encut in %s; do\n" % " ".join(encut)
        cmd += "  sed -e s/xxx/$encut/g INCAR_%s > INCAR\n" % dft
        cmd += "  " + execode + " >& vasp.out 2>&1\n\n"
        save_output_list.remove("WAVECAR")
        save_output_list.remove("CHGCAR")
        cmd += savevaspsh(
            save_output_list,
            dft+"_${encut}"
        )
        cmd += "done\n"
        functions.write_runsh(
            wkdir+"/run_conv_encut.sh",
            cmd
        )
    def dft_conv_kmesh(incar=incar):
        " generate input files for kmesh convergence test "
        poscar.write_struct(
            fpath=wkdir, 
            fname="POSCAR_init",
            calculator="poscar",
            dyn=dyn
        )

        incar.update(database.incar_scf)
        incar = update_incar(incar, params)
        write_incar(
            incar=incar,
            fname="%s/INCAR_%s" % (wkdir, dft)
        )

        kspacing = []
        if "encut" in params:
            kspacing = params["kspacing"].split("-")
        else:
            kspacing = [0.08, 0.1, 0.12, 0.14, 0.15, 0.2, 0.25, 0.3, 0.35]
        for _ in kspacing:
            vasp_kpoints_gen(
                poscar.lattice9, 
                _, 
                "none", 
                wkdir+"/KPOINTS_kconv_k%.2f" % _,
                "-1 -1 -1"
            )
        
        cmd = "# encut convergence test\n"
        cmd += "  cp POSCAR_init POSCAR\n"
        cmd += "  cp INCAR_%s INCAR\n" %dft
        cmd += "  for kspacing in %s; do\n" % " ".join([str(_) for _ in kspacing])
        cmd += "  cp KPOINTS_kconv_k$kspacing KPOINTS\n"
        cmd += "  " + execode + " >& vasp.out 2>&1\n\n"
        save_output_list.remove("WAVECAR")
        save_output_list.remove("CHGCAR")
        cmd += savevaspsh(
            save_output_list,
            dft+"_${kspacing}"
        )
        cmd += "done\n"
        functions.write_runsh(
            wkdir+"/run_conv_kmesh.sh",
            cmd
        )
    
    # =========================================================================
    # dos, band, effective mass
    # =========================================================================
    def dft_band(incar=incar):
        """ """
        print("Make sure WAVECAR_scf and CHGCAR_scf are in the same directory, otherwise, the calculation will be perfomed with dense k-mesh and comsume more time.")
        print("Gnerate the k-path by yourself and named as KPOINTS_band")
        
        incar.update(database.incar_band)
        incar = update_incar(incar, params)
        write_incar(
            incar=incar,
            fname="%s/INCAR_%s" % (wkdir, dft)
        )

        cmd = "# band calculation\n"
        cmd += "  cp INCAR_%s INCAR\n" % dft
        cmd += "  cp KPOINTS_band KPOINTS\n"
        cmd += "  cp WAVECAR_scf WAVECAR\n"
        cmd += "  cp CHGCAR_scf CHGCAR\n"
        cmd += "  " + execode + " >& vasp.out 2>&1\n\n"
        save_output_list.remove("WAVECAR")
        save_output_list.remove("CHGCAR")
        cmd += savevaspsh(
            save_output_list,
            dft
        )
        functions.write_runsh(
            wkdir+"/run_band.sh",
            cmd
        )
    
    def dft_dos(incar=incar):
        """ """
        incar.update(database.incar_dos)
        incar = update_incar(incar, params)
        
        kspacing, kmesh, kfix = update_kpoints_para(params)
        vasp_kpoints_gen(
            poscar.lattice9,
            kspacing,
            kmesh,
            wkdir + "KPOINTS_dos", 
            kfix
        )
        write_incar(
            incar=incar,
            fname="%s/INCAR_%s" % (wkdir, dft)
        )

        cmd = "# dos calculation\n"
        cmd += "  cp INCAR_%s INCAR\n" % dft
        cmd += "  cp KPOINTS_dos KPOINTS\n"
        cmd += "  cp WAVECAR_scf WAVECAR\n"
        cmd += "  cp CHGCAR_scf CHGCAR\n"
        cmd += "  " + execode + " >& vasp.out 2>&1\n\n"
        save_output_list.remove("WAVECAR")
        save_output_list.remove("CHGCAR")
        cmd += savevaspsh(
            save_output_list,
            dft
        )
        functions.write_runsh(
            wkdir+"/run_dos.sh",
            cmd
        )
    
    # =========================================================================
    # frequency: zpe and entropy
    # =========================================================================
    def dft_freq(incar=incar):
        # enable WAVECAR by default

        poscar.write_struct(
            fpath=wkdir, 
            fname="POSCAR_init",
            calculator="poscar",
            dyn=dyn
        )

        incar.update(database.incar_freq)
        incar = update_incar(incar, params)
        incar.pop("NCORE")
        write_incar(
            incar=incar,
            fname="%s/INCAR_%s" % (wkdir, dft)
        )

        kspacing, kmesh, kfix = update_kpoints_para(params)
        vasp_kpoints_gen(
            poscar.lattice9,
            kspacing,
            kmesh,
            wkdir+"/KPOINTS_freq", 
            kfix
        )
        cmd = "# freq \n"
        cmd += initvaspsh(dft, execode, 0, False)
        cmd += savevaspsh(
            update_save_list(
                save_output_list,
                dft,
                incar
            ),
            dft
        )
        functions.write_runsh(
            wkdir+"/run_freq.sh",
            cmd
        )


    # =========================================================================
    # aimd: NVT-NPT-NVE
    # =========================================================================
    def dft_aimd(incar=incar, aimd="nvt"):
        """ write """
        poscar.write_struct(
            fpath=wkdir, 
            fname="POSCAR_init",
            calculator="poscar",
            dyn=dyn
        )
        if aimd == "nvt":
            incar.update(database.incar_nvt)
        elif aimd == "nve":
            incar.update(database.incar_nve)

        # for multi-scf: need inputs from kwargs: ["params2", "params3"]
        if "mulmd" in params:
            mulstep = int(params["mulmd"]) + 1
            step = 1

            incar1 = update_incar(incar, params)
            write_incar(
                incar=incar1,
                fname="%s/INCAR_%s_%d" % (wkdir, dft, step)
            )
            kspacing, kmesh, kfix = update_kpoints_para(params)
            vasp_kpoints_gen(
                poscar.lattice9,
                kspacing,
                kmesh,
                wkdir + "KPOINTS_%s_%d" % (dft, step), 
                kfix
            )
            cmd = "# md step %d\n" % step
            cmd += initvaspsh(dft, execode, 1, False)
            cmd += savevaspsh(
                update_save_list(
                    save_output_list,
                    dft,
                    incar
                ),
                dft+"_%d" % step
            )
            functions.write_runsh(
                wkdir+"/run_%s_%d.sh" % (dft, step),
                cmd
            )

            # write following steps
            for step in range(2, mulstep):
                params2 = functions.parser_inputpara(
                    kwargs["params%d" % step]
                )
                incar2 = update_incar(incar, params2)
                write_incar(
                    incar=incar2,
                    fname="%s/INCAR_%s_%d" % (wkdir, dft, step)
                )

                kspacing, kmesh, kfix = update_kpoints_para(params2)
                vasp_kpoints_gen(
                    poscar.lattice9,
                    kspacing,
                    kmesh,
                    wkdir + "KPOINTS_md_%d" % step, 
                    kfix
                )
                cmd = "# md step %d\n" % step
                cmd += initvaspsh(dft, execode, step, True)
                cmd += savevaspsh(
                    update_save_list(
                        save_output_list,
                        dft,
                        incar
                    ),
                    dft+"_%d" % step
                )
                functions.write_runsh(
                    wkdir+"/run_%s_%d.sh" % (dft, step),
                    cmd
                )
        else:
            incar = update_incar(incar, params)
            
            write_incar(
                incar=incar,
                fname="%s/INCAR_%s" % (wkdir, dft)
            )
            kspacing, kmesh, kfix = update_kpoints_para(params)
            vasp_kpoints_gen(
                poscar.lattice9,
                kspacing,
                kmesh,
                wkdir+"/KPOINTS_md", 
                kfix
            )
            cmd = "# md\n"
            cmd += initvaspsh(dft, execode, 0, False)
            cmd += savevaspsh(
                update_save_list(
                    save_output_list,
                    dft,
                    incar
                ),
                dft
            )
            functions.write_runsh(
                wkdir+"/run_%s.sh" % dft,
                cmd
            )


    #print(incar)
    if dft == "scf":
        dft_scf()
    elif dft == "opt":
        dft_opt()
    elif dft == "band":
        dft_band()
    elif dft == "dos":
        dft_dos()
    elif dft == "conv_encut":
        dft_conv_encut()
    elif dft == "conv_kmesh":
        dft_conv_kmesh()
    elif dft == "freq":
        dft_freq()
    elif dft == "nvt":
        dft_aimd(aimd="nvt")
        

