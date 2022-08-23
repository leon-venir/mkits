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
 mkits
 |
 |----vasp_init 
 |      |
 |
 |----vasp_post
 |
 |
 |----wien_init
 |
 |
 |----fdmnes
 |
 ```

# 3. release history

+ version 0.1 (August 19th 2022)
  
  1. add WIEN2k initial interface: wien_init;
  2. add XDATCAR extrctor to vasp_post.

+ version 0.0 (April 1st 2022)
  
  1. vasp_init: add VASP initial interface
  2. VASP post-processing interface: vasp_post
