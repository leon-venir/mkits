# Welcome to my DFT kits: mkits

mkits is a python written tool containing many helpful initial- or post-processing commands for some popular first-principles calculation codes.


# INSTALLATION

``` pip install mkits ```

# VASP

## External 
### vasp_gen_input

This is a function to generate VASP input files. 


``` --potpath ```

``` --params ```

1. ```charged``` 
   
   Set a charged system. charged=-1 means one additional electron in the system.

2. ```mulisif```  

   if mulisif specified, another 

   need inputs from kwargs: ["params2", "params3"]

3. ```ggau```
   
4. ```INCAR-tag```