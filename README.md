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
 |        |----❀---kgen             : generate odd KPOINTS
 |        |----❀---vasp_gen_input   : generate input files for VASP
 |        |--------dryrun           : generate inputs for VASP, show details
 |        |--------dft              : choose the functional
 |        |--------potpath          : path to directory containing POTCAR
 |        |--------poscar           : specify the structure file
 |        |--------prec             : specify the calculation precision
 |        |--------wpath            : specify the working path
 |        |--------execode          : specify the mpirun code
 |        |----❀---gen2d            : generate 2dimensional structures
 |
 |----vasp_post
 |        |--------param            : additional parameters
 |        |--------fname            : input file name
 |        |----❀---extract_band     : extract eigen value from vasprun.xml
 |        |----❀---extract_dos      : extract total and partial DOS from vasprun.xml
 |        |----❀---structdiff       : check the difference between 2 structures in angstrom
 |        |----❀---extract_conv_test: extract data from convergence test (vasp_init)
 |        |--------wpath            : specify the working path
 |        |----❀---extract_xdatcar  : extract the i-th structure from XDATCAR
 |        |----❀---mse_xdatcar      : calculate the mean squared error
 |
 |----wien_init
 |        |--------fpath            : specify the working directory
 |        |--------fname            : input file name
 |        |--------struct           : input structure name
 |        |--------params           : additional parameters
 |        |----❀---gen_born_struct  : generate input files for BORN effective charge
 |
 |----wien_post
 |
 |----structgen
 |        |--------node_num         : maximum nodes number of the layer
 |        |--------node_type        : trigonal, layered structure
 |        |--------block_num        : give the block number
 |        |--------fix_block        : fix block
 |        |--------block_type       : specify the block
 |        |--------random_kation    : give the kation
 |        |--------random_anion     : give the anion
 |        |--------test_strgen      : do some test
 |
 |----boltz2
 |
 |----critic2
 |        |----❀---extract_gvh      : extract the data from
 |        |--------inp              : input files
 |        |--------out              : specify the name of the output file
 |
 |----fdmnes
 |        |----❀---gen_input        : generate the input file for FDMNES
 |        |--------struct_file      : path to the structure file
 |        |--------fpath            : path to save output files
 |        |--------findex           : specify the index
 |        |--------params           : additional parameters
 ```


# 3. release history

+ version 0.5 (September 14th 2022)
  1. fix the error when write wien2k structure

+ version 0.3 (August 27th 2022)

  1. fix some errors

+ version 0.2 (August 23rd 2022)

  1. add layered structures generator: structgen
  2. fix some errors

+ version 0.1 (August 19th 2022)
  
  1. add WIEN2k initial interface: wien_init;
  2. add XDATCAR extrctor to vasp_post.

+ version 0.0 (April 1st 2022)
  
  1. vasp_init: add VASP initial interface
  2. VASP post-processing interface: vasp_post


# 4. Introduction

## QE


## VASP

### Generate input files

1. convergence test
2. scf
3. band & dos