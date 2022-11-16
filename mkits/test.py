from mkits.globle import *
from mkits.qe import *

inp = "/Users/weiliangma/Works/Personal/tb/order/o081122_Cu_Cl_youan/cal_vasp/cu_bulk/CONTCAR_opt_rev-vdW-DF2_isif3_4"

a = qe_get_struct(inp)

print(a)
