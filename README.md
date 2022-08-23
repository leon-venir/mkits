## welcome to my DFT kits: mkits

mkits is a python written tool containing many helpful initial- or post-processing commands for some popular first-principles calculation codes.

# 1. installation

+ install from pip
  ```
  pip install mkits
  ```
  

+ build and install from source
  
  ```
  sudo python3 setup.py bdist_wheel
  pip3 install dist/thebuild.whl
  ```


 # 2. available functionals
 
 ```
 mkits: my DFT helper
 |
 |----vasp_init 
 |        |--------param            : additional parameters
 |        |--------kgen             : generate odd KPOINTS
 |        |--------vasp_gen_input   : generate input files for VASP
 |        |--------dryrun           : generate inputs for VASP, show details
 |        |--------dft              : choose the functional
 |        |--------potpath          : path to directory containing POTCAR
 |        |--------poscar           : specify the structure file
 |        |--------prec             : specify the calculation precision
 |        |--------wpath            : specify the working path
 |        |--------execode          : specify the mpirun code
 |        |--------gen2d            : generate 2dimensional structures
 |
 |----vasp_post
 |        |--------param            : 
 |        |--------fname            :
 |        |--------extract_band     :
 |        |--------extract_dos      :
 |        |--------structdiff       :
 |        |--------extract_conv_test:
 |        |--------wpath            :
 |        |--------extract_xdatcar  :
 |        |--------mse_xdatcar      :
 |
 |----wien_init
 |
 |
 |----wien_post
 |
 |
 |----structgen
 |
 |----boltz2
 |
 |----critic2
 |
 |----fdmnes
 ```


# 3. release history

+ version 0.2 (August 23rd 2022)

  1. add layered structures generator: structgen

+ version 0.1 (August 19th 2022)
  
  1. add WIEN2k initial interface: wien_init;
  2. add XDATCAR extrctor to vasp_post.

+ version 0.0 (April 1st 2022)
  
  1. vasp_init: add VASP initial interface
  2. VASP post-processing interface: vasp_post
