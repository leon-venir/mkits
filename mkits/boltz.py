import numpy as np
from mkits.sysfc import *
from mkits.globle import *


"""
:func split_boltz_dope: 
"""


def split_boltz_dope(dopefile, para_dict={"key1":"attrib1"}):
    """
    :para dopefile: input file from BoltzTraP2
    :para para_dict: a dictionary 
                     for case.dope.trace and case.dope.condtens
                     {
                         "nmin": "17", # n-type minimum
                         "nmax": "22", # n-type maximum
                         "nnum": "51", # number of n-type doping level
                         "pmin": "17", # p-type minimum
                         "pmax": "22", # p-type maximum
                         "pnum": "51", # number of p-type doping level
                    }

    """
    
    # check inputs
    if ".trace" in dopefile:
        condtens = False
    elif ".condtens" in dopefile:
        condtens = True
    else:
        lexit("Make sure that the input file contains dope.trace or dope.condtens extension.")
    try:
        data = np.loadtxt(dopefile, skiprows=1)
    except: lexit("Cannot find "+dopefile)

    # data processing
    # add power factor columns
    if condtens:
        # Sigma                       Seebeck/tau                   
        # xx xy xz yx yy yz zx zy zz  xx xy xz yx yy yz zx zy zz
        # 3  4  5  6  7  8  9  10 11  12 13 14 15 16 17 18 18 20
        pf = data[:,3:12]**2*data[:,12:21]
        data = np.hstack((data, pf))
        unitheads = "# doping[cm-3]  T[K]  N[e/uc]  sigma/tau0[1/(ohm*m*s)] xx xy xz yx yy yz zx zy zz  S[V/K] xx xy xz yx yy yz zx zy zz  kappae/tau0[W/(m*K*s)]  xx xy xz yx yy yz zx zy zz  powerfactor/tau0 xx xy xz yx yy yz zx zy zz\n" 
    else:
        # Sigma  Seebeck/tau                   
        # 4      5
        pf = data[:,4]**2*data[:,5]
        data = np.hstack((data, pf.reshape((-1,1))))
        unitheads = "# doping[cm-3]  T[K]  N[e/uc]  DOS(ef)[1/(Ha*uc)]  S[V/K]  sigma/tau0[1/(ohm*m*s)]  RH[m**3/C]  kappae/tau0[W/(m*K*s)]   cv[J/(mol*K)]   chi[m**3/mol]   powerfactor/tau0\n"

    # get temperature number
    temperature_num = 1
    for i in data[1:,1]:
        if i != data[0,1]: 
            temperature_num += 1
        else:
            break
    # reshape 
    data = np.reshape(data, (-1,temperature_num, len(data[1,:])))
    # get temperature and doping
    temperature = data[0,:,1]

    if "dope" in dopefile:
        try:
            doping = np.hstack((np.linspace(float(para_dict["nmax"]),float(para_dict["nmin"]),int(para_dict["nnum"])), np.linspace(float(para_dict["pmin"]),float(para_dict["pmax"]),int(para_dict["pnum"]))))
            npsplit = int(para_dict["nnum"])
        except:
            lexit("Error in additional parameters.")
        for i in range(len(temperature)):
            data[:,i,0] = doping
        boltzname = "dope"
    else:
        # npsplit = int(len(data)/2)
        npsplit = 0
        boltzname = "inte"

    if condtens:
        # write p type
        with open("ptype_%s.condtens" % boltzname, "w") as f:
            f.write(unitheads)
            for i in range(len(temperature)):
                f.write("# "+str(temperature[i])+" K\n")
                np.savetxt(f, data[npsplit:, i, :], fmt="%10.7f %10.2f %12.5f %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e")
                f.write("\n\n")
        # write n type
        with open("ntype_%s.condtens" % boltzname, "w") as f:
            f.write(unitheads)
            for i in range(len(temperature)):
                f.write("# "+str(temperature[i])+" K\n")
                np.savetxt(f, data[:npsplit, i, :], fmt="%10.7f %10.2f %12.5f %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e")
                f.write("\n\n")
    else:
        # write p type
        with open("ptype_%s.trace" % boltzname, "w") as f:
            f.write(unitheads)
            for i in range(len(temperature)):
                f.write("# "+str(temperature[i])+" K\n")
                np.savetxt(f, data[npsplit:, i, :], fmt="%10.7f %10.2f %12.5f %12.3f %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e")
                f.write("\n\n")
        # write n type
        with open("ntype_%s.trace" % boltzname, "w") as f:
            f.write(unitheads)
            for i in range(len(temperature)):
                f.write("# "+str(temperature[i])+" K\n")
                np.savetxt(f, data[:npsplit, i, :], fmt="%10.7f %10.2f %12.5f %12.3f %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e")
                f.write("\n\n") 
        # fmt="%10.3f %10.2f %12.5f %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e"

    return 0


def shift_eigin_wien(energyso, shift, val_ind, direction="cb"):
    # shift = 0.2 eV
    # direction = cb for conduction bands and vb for valence bands
    try:
        with open(energyso) as f: 
            ene_lines = f.readlines()
    except:
        lexit("Cannot find the energyso file.")
    
    with open("new.energyso", "w") as f:
        for line in ene_lines:
            if line[:8] == "        " and int(line[:12]) > val_ind:
                new_states = float(line[12:])+shift*0.0734986176
                f.write(line[:12]+"{:>19.15f}".format(new_states)+"     \n")
            else:
                f.write(line)
    return 0

