import sys

sys.path.insert(0, '../../')
from fastsemsim.SemSim import *
import constants
import argparse
from evaluation.homogeneity import main

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--pf', dest='pf', default=10)
    parser.add_argument('--base_folder', dest='base_folder', default=os.path.join(constants.OUTPUT_GLOBAL_DIR, "oob"))
    parser.add_argument('--sim_method', dest='sim_method', default="Resnik")
    parser.add_argument('--file_format', dest='file_format', default="emp_diff_modules_{}_{}_passed_oob.tsv")
    parser.add_argument('--cutoffs', dest='cutoffs',
                        default="1.0,2.0,3.0,4.0")
    parser.add_argument('--recalc_module_report', dest='recalc_module_report',
                        default="true")

    args = parser.parse_args()
    base_folder = args.base_folder
    sim_method = args.sim_method
    file_format = args.file_format
    pf = int(args.pf)
    cutoffs = np.array(args.cutoffs.split(','), dtype=float)
    recalc_module_report = args.recalc_module_report == "true"

    datasets = ["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"]
    algos = ["DOMINO2", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    prefix = "GE"
    main(prefix, base_folder, sim_method, file_format, pf, datasets, algos, cutoffs, recalc_module_report)

    datasets = ["brca", "crh", "scz", "tri", "t2d", "cad", "bmd", "hgt", "amd", "af"]
    algos = ["DOMINO2", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    prefix = "PASCAL_SUM"
    main(prefix, base_folder, sim_method, file_format, pf, datasets, algos, cutoffs, recalc_module_report)


