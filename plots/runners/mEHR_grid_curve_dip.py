import matplotlib

matplotlib.use('Agg')

import sys
import os
sys.path.insert(0, '../../')

import argparse
import seaborn as sns
import constants

sns.set(color_codes=True)
from plots.mEHR_grid_curve import plot_modules_ehr_summary

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets',
                        default="Breast_Cancer.G50,Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50")
    parser.add_argument('--prefix', dest='prefix', default="PASCAL_SUM")
    parser.add_argument('--base_folder_format', dest='base_folder_format',
                        default=os.path.join(constants.OUTPUT_GLOBAL_DIR, "oob"))
    parser.add_argument('--terms_file_name_format', dest='terms_file_name_format',
                        default="emp_diff_modules_{}_{}_passed_oob.tsv")
    parser.add_argument('--algos', dest='algos',
                        default="jactivemodules_greedy,jactivemodules_sa,netbox,bionet,dcem")

    args = parser.parse_args()

    datasets = args.datasets.split(",")
    algos = args.algos.split(",")
    prefix = args.prefix
    base_folder_format = args.base_folder_format
    terms_file_name_format = args.terms_file_name_format

    datasets = ["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"]
    algos = ["DOMINO2", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    prefix = "GE"
    plot_modules_ehr_summary(prefix, datasets, algos)

    datasets = ["brca", "crh", "scz", "tri", "t2d", "cad", "bmd", "hgt", "amd", "af"]
    algos = ["DOMINO2", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    prefix = "PASCAL_SUM"
    plot_modules_ehr_summary(prefix, datasets, algos)

