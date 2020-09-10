import sys

sys.path.insert(0, '../../')
from evaluation.robustness_f1 import main as main_f1
from evaluation.robustness_aupr import main as main_aupr
from evaluation.robustness_empty_solutions import main as main_empty
import constants
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--n_start', help="number of iterations (total n permutation is pf*(n_end-n_start))",
                        dest='n_start', default=0)
    parser.add_argument('--n_end', help="number of iterations (total n permutation is pf*(n_end-n_start))",
                        dest='n_end', default=100)
    parser.add_argument('--ss_ratios', help="ss_ratios", dest='ss_ratios', default="0.4,0.3,0.2,0.1")
    parser.add_argument('--pf', help="parallelization_factor", dest='pf', default=10)
    parser.add_argument('--hg_th', help="hg_th", dest='hg_th', default=0.05)
    parser.add_argument('--base_folder', help="base_folder", dest='base_folder', default=constants.OUTPUT_GLOBAL_DIR)

    args = parser.parse_args()

    n_start = int(args.n_start)
    n_end = int(args.n_end)
    ss_ratios = [float(a) for a in args.ss_ratios.split(",")]
    hg_th = float(args.hg_th)
    base_folder = args.base_folder
    parallelization_factor = int(args.pf)

    datasets = ["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"]
    algos = ["DOMINO3", "netbox3"]
    prefix = "GE"
    main_f1(prefix, datasets, algos, parallelization_factor, n_start, n_end, ss_ratios, hg_th, base_folder)
    main_aupr(prefix,datasets,algos,n_start,n_end,ss_ratios)
#     main_empty(prefix, datasets, algos, parallelization_factor, n_start, n_end, ss_ratios)
