# Strcuture builder 

import numpy as np
from mkits.globle import *


class tiling(object):
    """

    """
    def __init__(self) -> None:
        pass
    



def carbontube(m:int=5, n:int=5, l:int=1, out:str="out", fmt:str="xyz", **kwargs):
    """
    A scripts to build carbon nano tube

    Properties -----
    fmt: str vasp xyz

    Output --------- 

    """

    # bond length
    cc_bond = 1.42165
    thick = cc_bond * np.sqrt(3)

    # parameters
    radium = 0
    vacuum = 15
    coor_tran = 0

    if fmt == "xyz":
        pass
    elif fmt == "vasp":
        print(kwargs)
        try: 
            vacuum  = float(kwargs["vacuum"])
        except:
            print("No vacuum specified, use 15 angstrom. If you want other value, specify with vaccum=15")
            exit()

    if m == n:
        """ """
        alpha = (np.pi/m)
        gamma = np.arctan((np.sin(alpha/2))/(1/2+np.cos(alpha/2)))
        beta = alpha-2*gamma
        radium = cc_bond / 2 / np.sin(gamma)
        coor_tran = vacuum/2 + radium

        
        cood_alpha = 0
        cood_xy = np.array([[coor_tran,coor_tran+radium,thick/2]])
        for i in range(m):
            # 1st point
            cood_alpha += 2*gamma
            x = np.sin(cood_alpha) * radium
            y = np.cos(cood_alpha) * radium
            cood_xy = np.vstack((cood_xy, np.array([[x,y,thick/2]])))
            # 2nd point
            cood_alpha += beta
            x = np.sin(cood_alpha) * radium
            y = np.cos(cood_alpha) * radium
            cood_xy = np.vstack((cood_xy, np.array([[x,y,0]])))
            # 3rd point
            cood_alpha += 2*gamma
            x = np.sin(cood_alpha) * radium
            y = np.cos(cood_alpha) * radium
            cood_xy = np.vstack((cood_xy, np.array([[x,y,0]])))
            # 4th point
            cood_alpha += beta
            x = np.sin(cood_alpha) * radium
            y = np.cos(cood_alpha) * radium
            cood_xy = np.vstack((cood_xy, np.array([[x,y,thick/2]])))


    # write to file
    with open("./tube_m%dn%dl%d.%s" % (m, n, l, fmt), "w") as f:
        if fmt == "vasp":
            lattice = np.array([[vacuum + 2*radium + 0.0001, 0, 0],
                                [0, vacuum + 2*radium + 0.0001, 0],
                                [0, 0, thick]])
            f.write("tube_m%dn%dl%d\n" % (m, n, l))
            f.write("1.000\n")
            np.savetxt(f, lattice, fmt="%20.15f%20.15f%20.15f")
            f.write("C\n")
            f.write("%d\n" % len(cood_xy[1:,:]))
            f.write("Cartesian\n")
            np.savetxt(f, cood_xy[1:,:] + np.array([coor_tran, coor_tran, 0]), fmt="%20.15f")
        elif fmt == "xyz":
            f.write("%d\n" % len(cood_xy[1:]))
            f.write("tube_m%dn%dl%d\n" % (m, n, l))
            for i in range(1, len(cood_xy)):
                f.write("%s%20.15f%20.15f%20.15f\n" % ("C", cood_xy[i][0], cood_xy[i][1], cood_xy[i][2]))
        else:
            pass