import sys

sys.path.insert(0, '../')

import numpy as np
import pandas as pd
import os

import matplotlib

matplotlib.use("Agg")

from rpy2.robjects import pandas2ri

pandas2ri.activate()

import constants

import argparse

from utils.go_similarity import calc_similarity_matrix

def main(datasets, algos, prefix, pf=10):
    if not os.path.exists(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "ds_2_alg_scores")):
        try:
            os.makedirs(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "ds_2_alg_scores"))
        except Exception, e:
            print "error while creating ds_2_alg_scores folder: {}".format(e)

    for cur_ds in datasets:
        cur_full_ds="{}_{}".format(prefix, cur_ds)
        print "cur ds: {}".format(cur_full_ds)
        constants.update_dirs(DATASET_NAME_u=cur_full_ds)
        algos_signals = []

        for i_algo, cur_alg in enumerate(algos):

            print "current solution: {} {}".format(cur_ds,cur_alg)

            try:
                emp_results = pd.read_csv(
                    os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "MAX",
                                 "emp_diff_modules_{}_{}_passed_oob.tsv".format(cur_full_ds[len(prefix) + 1:], cur_alg)),
                    sep='\t', index_col=0)

            except:
                print "file {} wasn't found".format(
                    "emp_diff_modules_{}_{}_passed_oob.tsv".format(cur_full_ds[len(prefix) + 1:], cur_alg))
                continue

            emp_results = emp_results.sort_values(by='emp_rank')
            emp_results_fdr = emp_results.dropna().loc[emp_results.dropna()["passed_oob_permutation_test"].apply(
                lambda a: np.any(np.array(a[1:-1].split(", ")) == "True")).values, :]["GO name"]

            algos_signals.append(len(emp_results_fdr.index))
            all_go_terms = emp_results_fdr.index.values

            sim_method='Resnik'
            cache_file = os.path.join(constants.CACHE_GLOBAL_DIR,
                                      "similarity_cache_emp_{}_{}_{}.npy".format(cur_ds, cur_alg, sim_method))
            calc_similarity_matrix(all_go_terms, all_go_terms, pf, cache_file=cache_file, sim_method=sim_method)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="Breast_Cancer.G50,Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50")
    parser.add_argument('--prefix', dest='prefix', default="PASCAL_SUM")
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy,jactivemodules_sa,bionet,netbox")
    parser.add_argument('--pf', dest='pf', default=5)
    args = parser.parse_args()

    prefix = args.prefix
    datasets =  args.datasets.split(",")
    algos = args.algos.split(",")
    pf = int(args.pf)
    print "test"
    ds_summary = pd.DataFrame()
    main(datasets=datasets, algos=algos, prefix=prefix, pf=pf)
