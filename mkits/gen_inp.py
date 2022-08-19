import sys
import numpy as np

# python3 gen_inp.py ../../CONTCAR 10
poscar_file = sys.argv[1]
lattice = np.loadtxt(poscar_file, max_rows=3, skiprows=2)
coef = np.loadtxt(poscar_file, max_rows=1, skiprows=1)
a, c = lattice[0,0]*coef, lattice[-1,-1]*coef
lattice_param = np.array([a, a, c, 90, 90, 120])

with open(poscar_file, "r") as f:
  poscar_lines = f.readlines()

elements_map = {
  "Ti" : 22,
  "S"  : 16 
}

atoms_type = [elements_map[i] for i in poscar_lines[5].split()]
atoms_numb = [int(i) for i in poscar_lines[6].split()]
atoms_list = []
for i in range(len(atoms_numb)):
  atoms_list += [atoms_type[i]]*atoms_numb[i]

frac_pos = np.loadtxt(poscar_file, skiprows=8, max_rows=sum(atoms_numb))
frac_pos = np.hstack((np.array(atoms_list).reshape((-1,1)), frac_pos))

cal_index = sys.argv[2]
inp_file = "tis2_inp_%s.txt" % cal_index

inp_temp_1 = """! Fdmnes indata file
! Calculation for the ti K-edge in tis2
! Finite difference method calculation with convolution

 Filout
./out/tis2_%s
!   /home/mcsete/work/wma/fdmnes/parallel_fdmnes/

 Range                              ! Energy range of calculation (eV)
  -8. 0.5  10. 1.  18. 2. 60.        ! first energy, step, intermediary energy, step ..., last energy
 !-8. 0.2  13. 0.5 18. 1. 50. 2 120. ! first energy, step, intermediary energy, step ..., last energy    

 Radius                             ! Radius of the cluster where final state calculation is performed
   5.0                              ! For a good calculation, this radius must be increased up to 6 or 7 Angstroems

 Crystal                            ! Periodic material description (unit cell)
""" % cal_index

inp_temp_2 = """! get the signal coming from each atom
! Allsite

 SCF

 SCFexc

 Screening
    .6

! Green
 
 Quadrupole


! Convolution keyword : broadening with a width increasing versus energy as an arctangent

 Convolution 

 End
"""

with open(inp_file, "w") as f:
  f.write(inp_temp_1)
  f.write("%15.10f%15.10f%15.10f%6.1f%6.1f%6.1f\n" % (a, a, c, 90, 90, 120))
  np.savetxt(f, frac_pos, fmt="%3d%15.10f%15.10f%15.10f")
  f.write(inp_temp_2)
