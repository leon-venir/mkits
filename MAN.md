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

1. ```mulmd```
   
   Type: int
   
   If mulmd specified, multistep of molecular dynamic calculation will enabled.
   
   Need additional inputs from kwargs: ["params1", "params2", ...], and the variable of params will act as shared parameters.
   

2. ```charged``` 
   
   Set a charged system. charged=-1 means one additional electron in the system.

3. ```mulisif```  

   if mulisif specified, another 

   need inputs from kwargs: ["params2", "params3"]

4. ```ggau```
   
   
5. ```INCAR-tag```

```kwargs```

1. ```vdw_kernel_path```
   The path to the file of Van de Waals kernel table, which is copied to the calculation folder if enable the vdw corrections. 

2. ```nbands``` 
   Specify the number of the bands, which is identical to NBANDS in --params. But the input value provided here is the multiple of the valence electrons. The default variable is 0.5 for non-spin orbit coupling calculation. 

3. ``` ```