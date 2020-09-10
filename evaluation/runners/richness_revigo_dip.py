import sys
sys.path.insert(0, '../../')

import os
import constants
import argparse
from evaluation.richness_revigo import compute_redundancy_for_solutions

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--pf', dest='pf', default=10)
    parser.add_argument('--base_folder', dest='base_folder', default=os.path.join(constants.OUTPUT_GLOBAL_DIR,"oob"))
    parser.add_argument('--file_format', dest='file_format', default="emp_diff_modules_{}_{}_passed_oob.tsv")
    parser.add_argument('--sim_method', dest='sim_method', default="Resnik")
    args = parser.parse_args()

    file_format = args.file_format
    base_folder= args.base_folder
    sim_method = args.sim_method
    cutoffs = [1.0, 2.0, 3.0, 4.0]
    pf=int(args.pf)

    datasets=["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"]
    algos=["DOMINO2", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    prefix="GE"
    compute_redundancy_for_solutions(prefix=prefix, datasets=datasets, algos=algos, pf=pf, base_folder=base_folder, file_format=file_format, sim_method=sim_method, cutoffs=cutoffs)

    datasets=["brca", "crh", "scz", "tri", "t2d", "cad", "bmd", "hgt", "amd", "af"]
    algos=["DOMINO2", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    prefix="PASCAL_SUM"
    compute_redundancy_for_solutions(prefix=prefix, datasets=datasets, algos=algos, pf=pf, base_folder=base_folder, file_format=file_format, sim_method=sim_method, cutoffs=cutoffs)
