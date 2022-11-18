from mkits.globle import *
from mkits.qe import *
import numpy as np

inp = "/Users/weiliangma/Works/Personal/tb/order/o081122_Cu_Cl_youan/cal_vasp/cu_bulk/CONTCAR_opt_rev-vdW-DF2_isif3_4"
upf_path = "/Users/weiliangma/Scripts/dft_code/qe/upf/SSSP_1.1.2_PBE_efficiency/"

#qe_geninput(calculation="scf", wpath="/Users/weiliangma/Desktop/", struct_inp=inp, upf_path=upf_path)

#print(qe_kgen(struct_inp=inp, kspacing=0.15, oddeven="even"))


#print(qe_output_extractor("/Users/weiliangma/Works/Personal/tb/order/o081122_Cu_Cl_youan/cal_qe/conv_test/conv_ecut/ecutwfc_50.out"))

#qe_conv_extractor("/Users/weiliangma/Works/Personal/tb/order/o081122_Cu_Cl_youan/cal_qe/conv_test/conv_ecut/")


#

#qe_geninput(calculation="conv_ecutwfc", wpath="/Users/weiliangma/Works/Personal/tb/order/o071122_CuO_absorption/conv_test", struct_inp="/Users/weiliangma/Works/Personal/tb/order/o071122_CuO_absorption/struct/CONTCAR_bulk_opted.vasp", upf_path=upf_path, params_file="/Users/weiliangma/Works/Personal/tb/order/o071122_CuO_absorption/calparams")


qe_conv_extractor(wkdir="/Users/weiliangma/Works/Personal/tb/order/o071122_CuO_absorption/conv_test/conv_ecutwfc/")