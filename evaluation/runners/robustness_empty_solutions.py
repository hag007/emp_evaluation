import sys

sys.path.insert(0, '../')
import seaborn as sns

sns.set(color_codes=True)
import logging

sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from infra import *
import pandas as pd
import numpy as np
import argparse
from pandas.errors import EmptyDataError
from utils.daemon_multiprocessing import MyPool, func_star
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
from utils.add_GO_terms_metadata_agg import get_all_genes_for_term, vertices
import multiprocessing

def get_enriched_terms(algo, dataset, x, ss_ratio):
 
        df_go = pd.DataFrame(columns=['GO id', 'qval', 'pval'])
        df_go_pvals = pd.DataFrame()
        df_go_pvals.index.name = "GO id" 
    
        report_folder=os.path.join(constants.ROBUSTNESS_SOLUTIONS_DIR, "sol_{}_{}_robustness_{}_{}".format(algo, dataset, x, ss_ratio),"report")
        n_modules=len(pd.read_csv(os.path.join(report_folder, "modules_summary.tsv"), sep='\t').index)  
        go_results = [os.path.join(report_folder, cur_module) for cur_module in os.listdir(report_folder) if
                  "separated_modules" in cur_module and int(cur_module.split("_")[1]) < n_modules]
        for cur in go_results:
            try:
                df_go = pd.concat((df_go, pd.read_csv(cur, sep='\t')))
                df_go_pvals = pd.concat((df_go_pvals, pd.read_csv(cur, sep='\t').set_index("GO id")['pval']),
                                        axis=1)
            except (EmptyDataError, KeyError):
                print(EmptyDataError)  
                return False

        df_go_pvals[df_go_pvals.isna()] = 1

        df_go_pvals = df_go_pvals.min(axis=1)

        n_genes = [len(get_all_genes_for_term(vertices, cur_go_id, cur_go_id, cur_go_id == cur_go_id)) for i, cur_go_id
                   in
                   enumerate(df_go_pvals.index.values)]
        n_genes_series = pd.Series(n_genes, index=df_go_pvals.index)

        filtered_go_sets = n_genes_series.loc[
            np.logical_and.reduce([n_genes_series.values > 5, n_genes_series.values < 500])].index.values
        df_go_pvals = df_go_pvals[df_go_pvals.index.isin(filtered_go_sets)]

        print "total n_genes with pval:{}/{}".format(np.size(df_go_pvals), 7435)
        hg_pvals = np.append(df_go_pvals, np.ones(7435 - np.size(df_go_pvals)))
        print(hg_pvals)
        fdr_results = fdrcorrection0(hg_pvals, alpha=0.05, method='indep', is_sorted=False)[0]
        if np.sum(fdr_results) > 0:
            return True
        return False

def count_empties(prefix, dataset, cur, algo, ss_ratio=0.4, empties=None):
    print "starting iteration: {}, {}, {}, {}".format(prefix, dataset, ss_ratio, cur)
    try:
        empties.append(get_enriched_terms(algo, dataset, cur, ss_ratio))
    except Exception, e:
        print e
        empties.append(False)


def main(datasets, algos, prefix,parallelization_factor,n_start,n_end,ss_ratios):
    for ss_ratio in ss_ratios:

        df_empties = pd.DataFrame()
        for dataset in datasets:

            for algo in algos:
                empties = multiprocessing.Manager().list()
                p = MyPool(parallelization_factor)

                params = [[count_empties,
                           [prefix, dataset, x, algo, ss_ratio, empties]] for x
                          in np.arange(int(n_start), int(n_end))]
                p.map(func_star, params)
                p.close()

                df_empties.loc[algo, dataset] = sum(empties) / float(len(empties))

                print "empties for dataset {} and algo {}: {}".format(dataset, algo, df_empties.loc[algo, datasets])

        fname = os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation",
                             "robustness_{}_{}_{}_matrix_empty.tsv".format(prefix, n_end, ss_ratio))
        df_empties.to_csv(fname, sep='\t')
        print "save file to: {}".format(fname)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="TNFa_2")
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy,netbox,my_netbox_td")
    parser.add_argument('--network', dest='network', default="dip.sif")
    parser.add_argument('--pf', help="parallelization_factor", dest='pf', default=20)
    parser.add_argument('--n_start', help="number of iterations (total n permutation is pf*(n_end-n_start))",
                        dest='n_start', default=0)
    parser.add_argument('--n_end', help="number of iterations (total n permutation is pf*(n_end-n_start))",
                        dest='n_end', default=100)
    parser.add_argument('--ss_ratio', help="ss_ratio", dest='ss_ratio', default="0.4,0.3,0.2,0.1")
    parser.add_argument('--override_permutations', help="takes max or all samples", dest='override_permutations',
                        default="false")
    parser.add_argument('--base_folder', help="base_folder", dest='base_folder',
                        default=os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr/MAX"))

    args = parser.parse_args()

    datasets = args.datasets.split(",")
    algos = args.algos.split(",")
    prefix = args.prefix
    network_file_name = args.network
    parallelization_factor = int(args.pf)
    n_start = int(args.n_start)
    n_end = int(args.n_end)
    ss_ratios = [float(a) for a in args.ss_ratio.split(",")]

    main(datasets, algos, prefix,parallelization_factor,n_start,n_end,ss_ratios)

