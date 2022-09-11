import numpy as np


"""
:constants
"""
cons_hbar = 1.0545718e-34  # m^2kg/s
cons_kb = 1.38064852e-23  # m^2kg/(Ks^2)
cons_pi = 3.14159265359
cons_emass = 9.10938356e-31  # electron mass kg
cons_echarge = 1.60217662e-19  # electron charge



"""
:units convertor
"""
uc_ev2j = 1.602176634e-19  # electonic volt to Joule
uc_ry2j = 2.179874099E-18  # Ry
uc_gpa2nm2 = 1e9  # GPa to N/m^2
uc_bohr2ang = 0.529177 # bohr to angstrom
uc_ang2m = 1e-10


"""
:database 
"""

# :database vasp: VASP input tempelates ====================== #
incar_glob = {
    "PREC": "Normal",
    "ENCUT": "450",
    "EDIFF": "1e-6",
    "NCORE": "4",
    "NELM": "150"
}
incar_opt = {
    "IBRION": "2",
    "ISIF": "3",
    "NSW": "500",
    "ISMEAR": "0",
    "SIGMA": "0.1",
    "POTIM": "0.5"
}
incar_scf = {
    "IBRION": "-1",
    "NSW": "0",
    "ISMEAR": "0",
    "SIGMA" : "0.05"
}
incar_dos = {
    "ICHARG": "11", 
    "LORBIT": "11", 
    "IBRION": "-1", 
    "NSW": "0", 
    "ISMEAR": "-5", 
    "SIGMA": "0.1", 
    "NEDOS": "2500", 
    "NBANDS": "36", 
    "LMAXMIX": "4"
}
incar_md = {}


incar_prop = {
    "xanes": {
        "#XANES": "",
        "ALGO": "FAST",
        "ISMEAR": "0",
        "SIGMA": "0.05",
        "ICORELEVEL": "2",
        "CLNT": "xxx",
        "CLN": "1",
        "CLL": "0",
        "CLZ": "1.0",
        "CH_LSPEC": ".TRUE.",
        "CH_SIGMA": "0.3",
        "CH_NEDOS": "1000",
        "LREAL": "Auto"
    },
    "dfpt": {
        "#DFPT": "input for DFPT",
        "PREC": "Accurate",
        "ENCUT": "450",
        "EDIFF": "1e-8",
        "IBRION": "-1",
        "NSW": "0",
        "NELM": "200",
        "ISMEAR": "0",
        "SIGMA": "0.01",
        "IALGO": "38",
        "LREAL": ".FALSE." 
    },
    "born": {
        "#BORN": "effective charge with DFPT",
        "PREC": "Accurate",
        "ISMEAR": "0",
        "SIGMA": "0.01",
        "EDIFF": "1.E-8",
        "IBRION": "8",
        "LEPSILON": ".TRUE."
    }
}


incar_functionals = {
    "pbelda": {
        "#GGA_or_LDA": "make sure use proper POTCARs"
    },
    "pbesol": {
        "#pbesol": "INCAR for PBEsol, make sure use PBE POTCAR",
        "GGA": "PS"
    },
    "rev-vdW-DF2": {
        "#vdw-rev-vdW-DF2": "INCAR for rev-vdW-DF2, make sure use PBE POTCAR",
        "GGA": "MK",
        "LUSE_VDW": ".TRUE.",
        "AGGAC": "0.0",
        "LASPH": ".TRUE.",
        "PARAM1": "0.1234",
        "PARAM2": "0.711357",
        "ZAB_VDW": "-1.8867"
    },
    "optB88": {
        "#vdw-optB88": "INCAR for optB88, make sure use PBE POTCAR",
        "GGA": "BO",
        "PARAM1": "0.1833333333",
        "PARAM2": "0.22",
        "LUSE_VDW": ".TRUE.",
        "AGGAC": "0.0",
        "LASPH": ".TRUE."
    },
    "optPBE": {
        "#vdw-optPBE": "INCAR for optPBE, make sure use PBE POTCAR",
        "GGA": "OR",
        "LUSE_VDW": ".TRUE.",
        "AGGAC": "0.0",
        "LASPH": ".TRUE."
    },
    "hsesol": {
        "#hsesol": "INCAR for HSEsol, make sure use PBE POTCAR",
        "LHFCALC": ".TRUE.",
        "GGA": "PS",
        "HFSCREEN": "0.2",
        "ALGO": "Damped",
        "AEXX": "0.25",
        "AGGAX": "0.75",
        "AGGAC": "1.0",
        "ALDAC": "1.0"
    },
    "hse06" : {
        "#hse06": "INCAR for HSE06, make sure use PBE POTCAR",
        "LHFCALC": ".TRUE.",
        "GGA": "PE",
        "HFSCREEN": "0.2",
        "ALGO": "Damped",
        "AEXX": "0.25",
        "AGGAX": "0.75",
        "AGGAC": "1.0",
        "ALDAC": "1.0"
    }
}

atom_data = [ 
    [  0, "N", "Name", "atomic mass", ""], # 0
    [  1, "H", "Hydrogen", 1.00794, "1s1"], # 1
    [  2, "He", "Helium", 4.002602, "1s2"], # 2
    [  3, "Li", "Lithium", 6.941], # 3
    [  4, "Be", "Beryllium", 9.012182], # 4
    [  5, "B", "Boron", 10.811], # 5
    [  6, "C", "Carbon", 12.0107], # 6
    [  7, "N", "Nitrogen", 14.0067], # 7
    [  8, "O", "Oxygen", 15.9994], # 8
    [  9, "F", "Fluorine", 18.9984032], # 9
    [ 10, "Ne", "Neon", 20.1797], # 10
    [ 11, "Na", "Sodium", 22.98976928], # 11
    [ 12, "Mg", "Magnesium", 24.3050], # 12
    [ 13, "Al", "Aluminium", 26.9815386], # 13
    [ 14, "Si", "Silicon", 28.0855], # 14
    [ 15, "P", "Phosphorus", 30.973762], # 15
    [ 16, "S", "Sulfur", 32.065], # 16
    [ 17, "Cl", "Chlorine", 35.453], # 17
    [ 18, "Ar", "Argon", 39.948], # 18
    [ 19, "K", "Potassium", 39.0983], # 19
    [ 20, "Ca", "Calcium", 40.078], # 20
    [ 21, "Sc", "Scandium", 44.955912], # 21
    [ 22, "Ti", "Titanium", 47.867, "1s2,2s2p6,3s2p6d2,4s2"], # 22
    [ 23, "V", "Vanadium", 50.9415], # 23
    [ 24, "Cr", "Chromium", 51.9961], # 24
    [ 25, "Mn", "Manganese", 54.938045], # 25
    [ 26, "Fe", "Iron", 55.845], # 26
    [ 27, "Co", "Cobalt", 58.933195], # 27
    [ 28, "Ni", "Nickel", 58.6934], # 28
    [ 29, "Cu", "Copper", 63.546], # 29
    [ 30, "Zn", "Zinc", 65.38], # 30
    [ 31, "Ga", "Gallium", 69.723], # 31
    [ 32, "Ge", "Germanium", 72.64], # 32
    [ 33, "As", "Arsenic", 74.92160], # 33
    [ 34, "Se", "Selenium", 78.96], # 34
    [ 35, "Br", "Bromine", 79.904], # 35
    [ 36, "Kr", "Krypton", 83.798], # 36
    [ 37, "Rb", "Rubidium", 85.4678], # 37
    [ 38, "Sr", "Strontium", 87.62], # 38
    [ 39, "Y", "Yttrium", 88.90585], # 39
    [ 40, "Zr", "Zirconium", 91.224], # 40
    [ 41, "Nb", "Niobium", 92.90638], # 41
    [ 42, "Mo", "Molybdenum", 95.96], # 42
    [ 43, "Tc", "Technetium", None], # 43
    [ 44, "Ru", "Ruthenium", 101.07], # 44
    [ 45, "Rh", "Rhodium", 102.90550], # 45
    [ 46, "Pd", "Palladium", 106.42], # 46
    [ 47, "Ag", "Silver", 107.8682], # 47
    [ 48, "Cd", "Cadmium", 112.411], # 48
    [ 49, "In", "Indium", 114.818], # 49
    [ 50, "Sn", "Tin", 118.710], # 50
    [ 51, "Sb", "Antimony", 121.760], # 51
    [ 52, "Te", "Tellurium", 127.60], # 52
    [ 53, "I", "Iodine", 126.90447], # 53
    [ 54, "Xe", "Xenon", 131.293], # 54
    [ 55, "Cs", "Caesium", 132.9054519], # 55
    [ 56, "Ba", "Barium", 137.327], # 56
    [ 57, "La", "Lanthanum", 138.90547], # 57
    [ 58, "Ce", "Cerium", 140.116], # 58
    [ 59, "Pr", "Praseodymium", 140.90765], # 59
    [ 60, "Nd", "Neodymium", 144.242], # 60
    [ 61, "Pm", "Promethium", None], # 61
    [ 62, "Sm", "Samarium", 150.36], # 62
    [ 63, "Eu", "Europium", 151.964], # 63
    [ 64, "Gd", "Gadolinium", 157.25], # 64
    [ 65, "Tb", "Terbium", 158.92535], # 65
    [ 66, "Dy", "Dysprosium", 162.500], # 66
    [ 67, "Ho", "Holmium", 164.93032], # 67
    [ 68, "Er", "Erbium", 167.259], # 68
    [ 69, "Tm", "Thulium", 168.93421], # 69
    [ 70, "Yb", "Ytterbium", 173.054], # 70
    [ 71, "Lu", "Lutetium", 174.9668], # 71
    [ 72, "Hf", "Hafnium", 178.49], # 72
    [ 73, "Ta", "Tantalum", 180.94788], # 73
    [ 74, "W", "Tungsten", 183.84], # 74
    [ 75, "Re", "Rhenium", 186.207], # 75
    [ 76, "Os", "Osmium", 190.23], # 76
    [ 77, "Ir", "Iridium", 192.217], # 77
    [ 78, "Pt", "Platinum", 195.084], # 78
    [ 79, "Au", "Gold", 196.966569], # 79
    [ 80, "Hg", "Mercury", 200.59], # 80
    [ 81, "Tl", "Thallium", 204.3833], # 81
    [ 82, "Pb", "Lead", 207.2], # 82
    [ 83, "Bi", "Bismuth", 208.98040], # 83
    [ 84, "Po", "Polonium", None], # 84
    [ 85, "At", "Astatine", None], # 85
    [ 86, "Rn", "Radon", None], # 86
    [ 87, "Fr", "Francium", None], # 87
    [ 88, "Ra", "Radium", None], # 88
    [ 89, "Ac", "Actinium", None], # 89
    [ 90, "Th", "Thorium", 232.03806], # 90
    [ 91, "Pa", "Protactinium", 231.03588], # 91
    [ 92, "U", "Uranium", 238.02891], # 92
    [ 93, "Np", "Neptunium", None], # 93
    [ 94, "Pu", "Plutonium", None], # 94
    [ 95, "Am", "Americium", None], # 95
    [ 96, "Cm", "Curium", None], # 96
    [ 97, "Bk", "Berkelium", None], # 97
    [ 98, "Cf", "Californium", None], # 98
    [ 99, "Es", "Einsteinium", None], # 99
    [100, "Fm", "Fermium", None], # 100
    [101, "Md", "Mendelevium", None], # 101
    [102, "No", "Nobelium", None], # 102
    [103, "Lr", "Lawrencium", None], # 103
    [104, "Rf", "Rutherfordium", None], # 104
    [105, "Db", "Dubnium", None], # 105
    [106, "Sg", "Seaborgium", None], # 106
    [107, "Bh", "Bohrium", None], # 107
    [108, "Hs", "Hassium", None], # 108
    [109, "Mt", "Meitnerium", None], # 109
    [110, "Ds", "Darmstadtium", None], # 110
    [111, "Rg", "Roentgenium", None], # 111
    [112, "Cn", "Copernicium", None], # 112
    [113, "Uut", "Ununtrium", None], # 113
    [114, "Uuq", "Ununquadium", None], # 114
    [115, "Uup", "Ununpentium", None], # 115
    [116, "Uuh", "Ununhexium", None], # 116
    [117, "Uus", "Ununseptium", None], # 117
    [118, "Uuo", "Ununoctium", None], # 118
    ]

symbol_map = {
    "H":1,
    "He":2,
    "Li":3,
    "Be":4,
    "B":5,
    "C":6,
    "N":7,
    "O":8,
    "F":9,
    "Ne":10,
    "Na":11,
    "Mg":12,
    "Al":13,
    "Si":14,
    "P":15,
    "S":16,
    "Cl":17,
    "Ar":18,
    "K":19,
    "Ca":20,
    "Sc":21,
    "Ti":22,
    "V":23,
    "Cr":24,
    "Mn":25,
    "Fe":26,
    "Co":27,
    "Ni":28,
    "Cu":29,
    "Zn":30,
    "Ga":31,
    "Ge":32,
    "As":33,
    "Se":34,
    "Br":35,
    "Kr":36,
    "Rb":37,
    "Sr":38,
    "Y":39,
    "Zr":40,
    "Nb":41,
    "Mo":42,
    "Tc":43,
    "Ru":44,
    "Rh":45,
    "Pd":46,
    "Ag":47,
    "Cd":48,
    "In":49,
    "Sn":50,
    "Sb":51,
    "Te":52,
    "I":53,
    "Xe":54,
    "Cs":55,
    "Ba":56,
    "La":57,
    "Ce":58,
    "Pr":59,
    "Nd":60,
    "Pm":61,
    "Sm":62,
    "Eu":63,
    "Gd":64,
    "Tb":65,
    "Dy":66,
    "Ho":67,
    "Er":68,
    "Tm":69,
    "Yb":70,
    "Lu":71,
    "Hf":72,
    "Ta":73,
    "W":74,
    "Re":75,
    "Os":76,
    "Ir":77,
    "Pt":78,
    "Au":79,
    "Hg":80,
    "Tl":81,
    "Pb":82,
    "Bi":83,
    "Po":84,
    "At":85,
    "Rn":86,
    "Fr":87,
    "Ra":88,
    "Ac":89,
    "Th":90,
    "Pa":91,
    "U":92,
    "Np":93,
    "Pu":94,
    "Am":95,
    "Cm":96,
    "Bk":97,
    "Cf":98,
    "Es":99,
    "Fm":100,
    "Md":101,
    "No":102,
    "Lr":103,
    "Rf":104,
    "Db":105,
    "Sg":106,
    "Bh":107,
    "Hs":108,
    "Mt":109,
    "Ds":110,
    "Rg":111,
    "Cn":112,
    "Uut":113,
    "Uuq":114,
    "Uup":115,
    "Uuh":116,
    "Uus":117,
    "Uuo":118,
    }

""" :database encut: recommended kinetic cut off of each element. """
kinetic_cutoff = [
    [  0, "vasp",   "qe", "", ""], # 0
    [ 14,    450], # 14 Si
    [ 22,    450], # 22 Ti
    [ 23,    500], # 
    [ 32,    450], # Ge
    [ 51,    500], # Sb
    [ 52,    500]  # Te
]

""" :database smearing: recommended smearing, ie sigma in VASP, """
smearing = [
    [  0, "vasp",   "qe", ""], # 0
    [ 14,   0.05,      0], # Si
    [ 22,   0.05], # Ti
    [ 23,   0.05], # 
    [ 32,   0.05], # Ge
    [ 51,   0.05], # Sb
    [ 52,   0.05]  # Te
]

""" :database VASP recommended PAW """
vasp_paw = [
    [  0, "PBE", "val_num",  "LDA"],
    [ 14, "Si", 6, "Si", 6],
    [ 32, "Ge_h", 14, "Ge_h", 14],
    [ 51, "Sb", 5, "Sb", 5],
    [ 52, "Te", 6, "Te", 6]
]


""" :database kpath_sg: recommended kpath for each space group. """
kpath_sg = {
    164: np.array([[0,0,0], [1,2,3]]),
    164+500: ["Gamma", ""]
}


""" :script """
""" :script gnu2dlines: """
gnu2dline = """set terminal pngcairo enhanced font "Times,24" lw 5 size 1000, 900
set output '%s'
unset key
set title '%s'
set xlabel "%s" font "Times,28"
set ylabel "%s" font "Times,28"
set grid
set colorsequence default
set style line 1 dt 1 lt 1 lw 0.8 pt 2 ps 3
plot "%s" with linespoints ls 1 notitle
"""

""" :script extract_vasp_xanes_outcar: """
extract_xanes_vasp = """#!/bin/bash

parallel=-1
normal=-1
all=-1
tauc=-1
trace=-1
while [[ $# -gt 0 ]]
do
   key="$1"
   case $key in
      -parallel) parallel=0
      ;;
      -normal) normal=0
      ;;
      -trace) trace=0
      ;;
      -tauc) tauc=0
      ;;
   esac
   shift
done


cat > helpscript.perl  <<EOF
#!/bin/perl

use strict;
use warnings;
my \$mode=shift;

while(<>)
{
   chomp;
   if(\$_ =~ /frequency dependent IMAGINARY DIELECTRIC FUNCTION/)
   {
      \$_=<>;
      \$_=<>;
      while (<>)
      {
         my \$sum=0;
         if (\$_ !~ /[0-9]/) {last;}
         chomp;
         \$_=~s/^/ /;
         my @help=split(/[\t,\s]+/);
         if (\$help[2]=~/NaN/||\$help[3]=~/NaN/||\$help[4]=~/NaN/) {next;}
         if (\$help[5]=~/NaN/||\$help[6]=~/NaN/||\$help[4]=~/NaN/) {next;}
         if (\$mode==0) {\$sum=\$help[2]+\$help[3]+\$help[4]+\$help[5]+\$help[6]+\$help[7];}
         if (\$mode==1) {\$sum=\$help[4];}
         if (\$mode==2) {\$sum=\$help[2]+\$help[3];}
         if (\$mode==3) {\$sum=\$help[2]+\$help[3]+\$help[4];}
         if (\$mode==4) {\$sum=(\$help[1]*\$help[1]*(\$help[2]+\$help[3]+\$help[4]+\$help[5]+\$help[6]+\$help[7]))**0.5;}
         if (\$mode==5) {\$sum=(\$help[1]*\$help[1]*(\$help[4]))**0.5;}
         if (\$mode==6) {\$sum=(\$help[1]*\$help[1]*(\$help[2]+\$help[3]))**0.5;}
         if (\$mode==7) {\$sum=(\$help[1]*\$help[1]*(\$help[2]+\$help[3]+\$help[4]))**0.5;}
         print \$help[1]," ",\$sum,"\n";
      }
   }
   last if eof;
}
EOF

if [[ $normal -eq 0 ]]; then
   if [[ $tauc -eq 0 ]]; then
      perl helpscript.perl 4 OUTCAR > xanes.dat
   else
      perl helpscript.perl 1 OUTCAR > xanes.dat
   fi
else
   if [[ $parallel -eq 0 ]]; then
      if [[ $tauc -eq 0 ]]; then
         perl helpscript.perl 5 OUTCAR > xanes.dat
      else
         perl helpscript.perl 2 OUTCAR > xanes.dat
      fi
   else
      if [[ $trace -eq 0 ]]; then
         if [[ $tauc -eq 0 ]]; then
            perl helpscript.perl 6 OUTCAR > xanes.dat
         else
            perl helpscript.perl 3 OUTCAR > xanes.dat
         fi
      else
         if [[ $tauc -eq 0 ]]; then
            perl helpscript.perl 7 OUTCAR > xanes.dat
         else
            perl helpscript.perl 0 OUTCAR > xanes.dat
         fi
      fi
   fi
fi
rm helpscript.perl
"""


""" :scipt monitor_cluster """
monitor_cluster = """import subprocess
import time
from datetime import datetime

"""
