import argparse
from mkits.fdmnes import fdmnes_gen_inp
from mkits.vasp import *
from mkits.globle import *
from mkits.sysfc import *
from mkits.critic2 import *
from mkits.wien import *
from mkits.structgen import *
from mkits.boltz import *

def vasp_init_parser(args):
    try:
        parameter = parser_inputpara(args.param)
    except:
        parameter = {}
    if args.gen2d:
        vasp_build_low_dimension(parameter[0], parameter[1], parameter[2], parameter[3])
    elif args.pot:
        poscar_dict = struct("POSCAR").return_dict()
        vasp_potcar_gen(poscar_dict, args.pot_path)
    elif args.kgen != 0:
        poscar_dict = struct("POSCAR").return_dict()
        vasp_kpoints_gen(poscar_dict, args.kgen, parameter["kmesh"])
    elif args.vasp_gen_input:
        vasp_gen_input(args.dft if args.dft else "scf", args.potpath if args.potpath else "./", args.poscar if args.poscar else "POSCAR", args.dryrun if args.dryrun else False, args.prec if args.prec else "Normal", args.wpath if args.wpath else "./", args.execode if args.execode else "mpirun -np $SLURM_NTASKS vasp_std", args.param if args.param else "gga=pbelda")
    elif args.gen_arb_klist:
        arbitrary_klist(kend=args.gen_arb_klist, kmesh=parameter["kmesh"] if "kmesh" in parameter else "7", fpath=parameter["fpath"] if "fpath" in parameter else "./", fname=parameter["fname"] if "fname" in parameter else "KPOINTS")

def vasp_post_parser(args):
    try:
        parameter = parser_inputpara(args.param)
    except:
        parameter = {}
    if args.extract_band:
        #vasp_band_extractor(eign=parameter["eignval"] if "eignval" in parameter else "EIGENVAL", poscar=parameter["poscar"] if "poscar" in parameter else "POSCAR", kpoint=parameter["kpoints"] if "kpoints" in parameter else "KPOINTS")
        vasp_band_extractor(fpath=parameter["fpath"] if "fpath" in parameter else "./", fname=parameter["fname"] if "fname" in parameter else "BANDCAR", xmlfile=args.extract_band)
    elif args.extract_dos:
        vasp_dos_extractor(args.extract_dos)
    elif args.structdiff:
        print(struct_diff(args.structdiff.split(",")[0],args.structdiff.split(",")[1]))
    elif args.extract_conv_test:
        extract_conv_test(args.wpath)
    elif args.extract_xdatcar:
        extract_xdatcar(xdatcar=args.extract_xdatcar, fpath=args.wpath if args.wpath else "./", idx=[int(_) for _ in parameter["list"].split("-")])
    elif args.mse_xdatcar:
        mse_xdatcar(xdatcar=args.mse_xdatcar, crys_struct=args.fname, fpath=args.wpath)
    elif args.getvbmcbm:
        [print(key, " : ", value) for key, value in getvbmcbm(args.getvbmcbm).items()]

def structgen_parser(args):
    if args.node_num and args.node_type and args.fix_block:
        gen_struct_fix_block(args.node_num, args.node_type, args.block_type)
    elif args.test_strgen:
        test_struct()

def boltz2_parser(args):
    if args.gnuplot and args.argpara:
        split_boltz_dope(args.inputdata, parser_inputpara(args.argpara))
    elif args.gnuplot:
        split_boltz_dope(args.inputdata)

def critic2_parser(args):
    if args.extract_gvh:
        critc2_extract_data(args.inp if args.inp else "more.cro", args.out if args.out else "more.crx")

def fdmnes_parser(args):
    if args.gen_input:
        fdmnes_gen_inp(args.struct_file, args.fpath, args.findex, args.params if args.params else "")

def wien_init_parser(args):
    try:
        parameters = parser_inputpara(args.params)
    except:
        pass
    if args.gen_born_struct:
        print("The struct file is: ", args.struct, "\nThe index of moved atoms are: ", [int(_) for _ in parameters["atomlist"].split("-")], "\nThe directions are: ", [*parameters["direction"]], "\nNew structure files are written into: ", args.fpath if args.fpath else "./", "\nThe displacement is: ", parameters["displacement"] if "displacement" in parameters else 0.01)
        gen_wien_born(structure=args.struct, atom_list_index=[int(_) for _ in parameters["atomlist"].split("-")], direction=[*parameters["direction"]], fpath=args.fpath if args.fpath else "./", displacement=float(parameters["displacement"]) if "displacement" in parameters else 0.01)

def wien_post_parser(args):
    pass

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="DFT helper"
    )
    parser.add_argument(
        "--version",
        action="version",
        version="0.5",
        help="print version information"
    )
    subparser = parser.add_subparsers(
        help="sub command:"
    )
    # ==========================================================================
    # subparser
    # ==========================================================================
    parser_vasp_init = subparser.add_parser("vasp_init", help="enable VASP initial interface")
    parser_vasp_post = subparser.add_parser("vasp_post", help="enable VASP post-processing interface")
    parser_wien_init = subparser.add_parser("wien_init", help="enable WIEN2k initial interface")
    parser_wien_post = subparser.add_parser("wien_post", help="enable WIEN2k post-processing interface")
    parser_structgen = subparser.add_parser("structgen", help="generate layer structure")
    parser_boltz2     = subparser.add_parser("boltz2", help="enable BoltzTraP2 post-processing iterface")
    parser_critic2    = subparser.add_parser("critic2", help="enable Critic2 post-processing iterface")
    parser_fdmnes     = subparser.add_parser("fdmnes", help="enable FDMNES init-/post-processing iterface")

    # option2function
    parser_vasp_init.set_defaults(func=vasp_init_parser)
    parser_vasp_post.set_defaults(func=vasp_post_parser)
    parser_wien_init.set_defaults(func=wien_init_parser)
    parser_wien_post.set_defaults(func=wien_post_parser)
    parser_structgen.set_defaults(func=structgen_parser)
    parser_boltz2.set_defaults(func=boltz2_parser)
    parser_critic2.set_defaults(func=critic2_parser)
    parser_fdmnes.set_defaults(func=fdmnes_parser)

    # ==========================================================================
    # vasp init argument
    # ==========================================================================
    parser_vasp_init.add_argument(
        "--param",
        action="store",
        type=str,
        help="The additional options, more details can be found in specific argument."
    )
    parser_vasp_init.add_argument(
        "--kgen",
        action="store",
        default=0,
        type=float,
        help="Generate odd KPOINTS following with kspacing, avail optional parameter: kmesh:[odd,even,none]"
    )
    parser_vasp_init.add_argument(
        "--vasp_gen_input",
        action="store_true",
        help="generate inputs for vasp, use --dryrun with --vasp_gen_input to show more details."
    )
    parser_vasp_init.add_argument(
        "--dryrun",
        action="store_true",
        help="generate inputs for vasp, show details."
    )
    parser_vasp_init.add_argument(
        "--dft",
        action="store",
        default="scf",
        type=str,
        help="Choose the calculation type: optional [opt, scf, scfband, conv_encut, conv_kmesh, conv_sigma]"
    )
    parser_vasp_init.add_argument(
        "--potpath",
        action="store",
        help="Path to directory containing POTCAR, rename all POTCARs to POTCAR_Te."
    )
    parser_vasp_init.add_argument(
        "--poscar",
        action="store",
        default="POSCAR",
        type=str,
        help="Specify the strucutre file: POSCAR"
    )
    parser_vasp_init.add_argument(
        "--prec",
        action="store",
        default="Normal",
        type=str,
        help="Choose the calculation precision: optional [Normal, Accurate, Low]"
    )
    parser_vasp_init.add_argument(
        "--wpath",
        action="store",
        default="./",
        type=str,
        help="Specify the working directory. default: ./"
    )
    parser_vasp_init.add_argument(
        "--execode",
        action="store",
        default="mpirun -np $SLURM_NTASKS vasp_std",
        type=str,
        help="Specify the execute code. default: mpirun -np $SLURM_NTASKS vasp_std"
    )
    parser_vasp_init.add_argument(
        "--gen2d",
        action="store_true",
        default=False,
        help="generate the 2-dimensional structure, specify parameters with --parameter: [x,y,z],[20],[fix,add] the direction,  vaccum in ang, fit or add vaccum"
    )
    parser_vasp_init.add_argument(
        "--pot",
        action="store_true",
        default=False,
        help="write POTCAR based on POSCAR, specify path to POTCAR database with --potpath"
    )
    parser_vasp_init.add_argument(
        "--gen_arb_klist",
        action="store",
        type=str,
        help="Generate arbitrary klist for effective mass calculation. --gen_arb_klist [0,0,0,0.03,0.03,0.0 for cuboid center_x_y_z and length of side, for line start and end coordinates] --param kmesh=7-7-7,fpath=./,fname=KPOINTS."
    )
    



    # ==========================================================================
    # vasp post argument
    # ==========================================================================
    parser_vasp_post.add_argument(
        "--param",
        action="store",
        type=str,
        help="The additional options, more details can be found in specific argument."
    )
    parser_vasp_post.add_argument(
        "--fname",
        action="store",
        type=str,
        help="Input file."
    )
    parser_vasp_post.add_argument(
        "--extract_band",
        action="store",
        help="Extract eigen value from vasprun, and export to gnuplot format data. eg --extract_band vasprun_band.xml or specify the input files with --param fpath=./,fname=BANDCAR"
    )
    parser_vasp_post.add_argument(
        "--extract_dos",
        action="store",
        help="Extract total and partial DOS from vasprun.xml, and export to plotable data."
    )
    parser_vasp_post.add_argument(
        "--structdiff",
        action="store",
        help="Check the difference between 2 structures after opt: structname1,structname2"
    )
    parser_vasp_post.add_argument(
        "--extract_conv_test",
        action="store_true",
        help="Extract data from vasp_init --vasp_gen_input --dft [conv_encut, conv_kmesh, ...], and generate gnuplot files"
    )
    parser_vasp_post.add_argument(
        "--wpath",
        action="store",
        default="./",
        type=str,
        help="Specify the working directory. default: ./"
    )
    parser_vasp_post.add_argument(
        "--extract_xdatcar",
        action="store",
        type=str,
        help="Extract the i-th strucutre from XDATCAR, need to specify the working directory with --wpath, and the index list of target structure. eg: --wpath ./ --param list=500-1000"
    )
    parser_vasp_post.add_argument(
        "--mse_xdatcar",
        action="store",
        type=str,
        help="Calculate the mean squared error for MD simulation. need to specify the original crytal with --fname. eg: --mse_xdatcar XDATCAR_2000fs --fname POSCAR_init --wpath ./"
    )
    parser_vasp_post.add_argument(
        "--getvbmcbm",
        action="store",
        type=str,
        help="Extract VBM and CBM from vasprun.xml file. eg --getvbmcbm vasprun_dos_pbesol.xml"
    )

    # ==========================================================================
    # vasp init argument
    # ==========================================================================
    parser_wien_init.add_argument(
        "--fpath",
        action="store",
        type=str,
        help="Specify the path to read or to write files."
    )
    parser_wien_init.add_argument(
        "--fname",
        action="store",
        type=str,
        help="Specify the file name."
    )
    parser_wien_init.add_argument(
        "--struct",
        action="store",
        type=str,
        help="Specify the input structure file."
    )
    parser_wien_init.add_argument(
        "--params",
        action="store",
        type=str,
        help="The additional options, more details can be found in specific argument."
    )
    parser_wien_init.add_argument(
        "--gen_born_struct",
        action="store_true",
        help="generate the structure file for BORN effective charge calculation; use with --struct and --params; parameters need to be specified are\n1. direction=[xyz],\n2. atomlist=[1-4-12],\n the optional parameter is displacement=0.01\n example: --params direction=xz,atomlist=1-3-5,displacement=0.01."
    )


    # ==========================================================================
    # structure generation argument
    # ==========================================================================
    parser_structgen.add_argument(
        "--node_num",
        action="store",
        type=str,
        help="The maximum nodes number of the layer structure: 2,3,4,5,6,7"
    )

    parser_structgen.add_argument(
        "--node_type",
        action="store",
        choices=["trigonal"],
        help="trigonal: layered structure\n"
    )

    parser_structgen.add_argument(
        "--block_num",
        type=str,
        action="store",
        help="give the block number seperated with comma without space: 2,3,5"
    )

    parser_structgen.add_argument(
        "--fix_block",
        action="store_true",
        help="Use the specific block, should give the block with —block_type"
    )

    parser_structgen.add_argument(
        "--block_type",
        type=str,
        action="store",
        help="Specify the block seperated with comma without space: Bi-Bi,Te-Bi-Te,Te-Bi-Te-Bi-Te"
    )

    parser_structgen.add_argument(
        "--random_kation",
        type=str,
        action="store",
        help="Give the kation seperated with comma: Pb,Bi"
    )

    parser_structgen.add_argument(
        "--random_anion",
        type=str,
        action="store",
        help="Give the kation seperated with comma: Te,Se"
    )

    parser_structgen.add_argument(
        "--test_strgen",
        action="store_true"
    )

    # ==========================================================================
    # BoltzTraP2 post-processing argument
    # ==========================================================================
    parser_boltz2.add_argument(
        "--gnuplot",
        action="store_true",
        help='Export data with gnuplot format. Need following arguments: --inputdata, --argpara. The additional parameters of doped file are give as following: "nmin": "17", # n-type minimum; "nmax": "22", # n-type maximum; "nnum": "51", # number of n-type doping level; "pmin": "17", # p-type minimum; "pmax": "22", # p-type maximum; "pnum": "51", # number of p-type doping level. And The additional parameters for integrated file are: '
    )

    parser_boltz2.add_argument(
        "--inputdata",
        action="store",
        help="input file: trace or condtens for --gnuplot."
    )

    parser_boltz2.add_argument(
        "--argpara",
        action="store",
        help='Additional parameters.'
    )


    # ==========================================================================
    # critic2 post-processing argument
    # ==========================================================================
    parser_critic2.add_argument(
        "--extract_gvh",
        action="store_true",
        help="Extract the data from [auto, qtree, pointprop gtf_kir/vtf_kir/htf_kir, cpreport verylong]."
    )
    parser_critic2.add_argument(
        "--inp",
        type=str,
        action="store",
        help="Input files"
    )
    parser_critic2.add_argument(
        "--out",
        type=str,
        action="store",
        help="The name of output files."
    )


    # ==========================================================================
    # FDMNES init/post-processing argument
    # ==========================================================================
    parser_fdmnes.add_argument(
        "--gen_input",
        action="store_true",
        help="Generate the input file for FDMNES, additional parameters: --struct_file [path/to/structure], --fpath [path/to/store/inputfiles], --findex [index], --params [other setting parameters]"
    )
    parser_fdmnes.add_argument(
        "--struct_file",
        type=str,
        action="store",
        help="Path to the structure file."
    )
    parser_fdmnes.add_argument(
        "--fpath",
        type=str,
        action="store",
        help="Path to save some files."
    )
    parser_fdmnes.add_argument(
        "--findex",
        type=int,
        action="store",
        help="Index of fdmnes input files: [compounds]_inp_[index].txt"
    )
    parser_fdmnes.add_argument(
        "--params",
        type=str,
        action="store",
        help="Additional setting parameters."
    )


    # ==========================================================================
    # option2function
    # ==========================================================================

    args = parser.parse_args()
    args.func(args)


def mkits_main():
    """Entry point for the mkits script."""
    parse_arguments()
