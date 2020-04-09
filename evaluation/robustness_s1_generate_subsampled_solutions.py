import sys

sys.path.insert(0, '../')
import seaborn as sns

sns.set(color_codes=True)
import logging

sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from infra import *
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from utils.randomize_data import create_random_ds
from utils.randomize_data import permutation_output_exists
import argparse
from pandas.errors import EmptyDataError
from runners.datasets_multithread_runner import run_dataset
from utils.daemon_multiprocessing import MyPool, func_star
import shutil
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
from utils.add_GO_terms_metadata_agg import get_all_genes_for_term, vertices
import multiprocessing
import random
import codecs


def subsample_dataset(dataset_name, network_file_name, ss_ratio=0.4):
    if dataset_name.startswith("GE_"):
        data_file_path = os.path.join(constants.DATASETS_DIR, dataset_name, "data", "ge.tsv")
    else:
        data_file_path = os.path.join(constants.DATASETS_DIR, dataset_name, "data", "score.tsv")

    df = pd.read_csv(data_file_path, index_col=0, sep='\t')

    network_genes = set(
        np.unique(pd.read_csv(os.path.join(constants.NETWORKS_DIR, network_file_name), sep='\t').values))
    assay_genes = set(df.index.values)
    overlapped_genes = list(network_genes.intersection(assay_genes))
    # np.random.seed(int(random.random() * 1000))
    np.random.seed(int(codecs.encode(os.urandom(4), 'hex'), 16))
    genes_to_omit = np.random.choice(np.array(overlapped_genes), int(len(overlapped_genes) * ss_ratio), replace=False)

    df = df.loc[~df.index.isin(genes_to_omit), :]
    df.to_csv(os.path.join(constants.DATASETS_DIR, dataset_name, "data", "ge.tsv"), sep='\t')


def get_enriched_terms(algos, datasets, is_plot=False, empirical_th=None):
    for cur_algo in algos:
        algos_filter = cur_algo

        df_go = pd.DataFrame(columns=['GO id', 'qval', 'pval'])
        df_go_pvals = pd.DataFrame()
        df_go_pvals.index.name = "GO id"
        for cur_ds in datasets:
            go_results = [os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo, cur_module) for cur_algo in
                          os.listdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds))
                          if os.path.isdir(
                    os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo)) and cur_algo == algos_filter for
                          cur_module in os.listdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, cur_ds, cur_algo)) if
                          "separated_modules" in cur_module]

            for cur in go_results:
                try:
                    df_go = pd.concat((df_go, pd.read_csv(cur, sep='\t')))
                    df_go_pvals = pd.concat((df_go_pvals, pd.read_csv(cur, sep='\t').set_index("GO id")['pval']),
                                            axis=1)
                except (EmptyDataError, KeyError):
                    pass
        # print df_go_pvals
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
        fdr_results = fdrcorrection0(hg_pvals, alpha=0.05, method='indep', is_sorted=False)[0]
        if np.sum(fdr_results) > 0:
            fdr_th = np.max(hg_pvals[fdr_results])
        else:
            fdr_th = 0
        print "fdr_th :{}".format(fdr_th)
        HG_CUTOFF = -np.log10(fdr_th)
        print "HG cutoff: {}".format(HG_CUTOFF)

        df_go_pvals = df_go_pvals.loc[df_go_pvals.values <= fdr_th]

        pval = -np.log10(df_go["pval"].values)
        if np.size(pval) == 0:
            pval = np.array([0])

        return pval, df_go, df_go_pvals


def recovery_iteration(prefix, dataset, cur, algo, network_file_name="dip.sif",
                       base_folder='/home/hag007/Desktop/aggregate_report/oob', ss_ratio=0.4, precisions=None,
                       recalls=None, override_permutations=False):
    print "starting iteration: {}, {}, {}".format(prefix, dataset, cur)
    # recovered_dataset_name="{}_{}".format(prefix, dataset)
    recovered_dataset_name = "{}_{}_recovery_{}_{}_{}".format(prefix, dataset, algo, ss_ratio, cur)

    if not os.path.exists(os.path.join(constants.OUTPUT_GLOBAL_DIR, recovered_dataset_name)) or override_permutations:

       if os.path.exists(os.path.join(constants.DATASETS_DIR, recovered_dataset_name)):
           shutil.rmtree(os.path.join(constants.DATASETS_DIR, recovered_dataset_name))
       # if os.path.exists(os.path.join(constants.OUTPUT_GLOBAL_DIR, recovered_dataset_name)):
       #     shutil.rmtree(os.path.join(constants.OUTPUT_GLOBAL_DIR, recovered_dataset_name))

       if not os.path.exists(os.path.join(constants.DATASETS_DIR, recovered_dataset_name)):
            shutil.copytree(os.path.join(constants.DATASETS_DIR, prefix + "_" + dataset),
                            os.path.join(constants.DATASETS_DIR, recovered_dataset_name))
            shutil.rmtree(os.path.join(constants.DATASETS_DIR, recovered_dataset_name, "output"))
            os.makedirs(os.path.join(constants.DATASETS_DIR, recovered_dataset_name, "output"))
            shutil.rmtree(os.path.join(constants.DATASETS_DIR, recovered_dataset_name, "cache"))
            os.makedirs(os.path.join(constants.DATASETS_DIR, recovered_dataset_name, "cache"))

            subsample_dataset(recovered_dataset_name, network_file_name, ss_ratio=ss_ratio)

       permuted_network_file_name = network_file_name  # # _perm

       run_dataset(recovered_dataset_name, score_method=score_method,
                   algos=[algo], network_file_name=permuted_network_file_name)


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
    parser.add_argument('--ss_ratio', help="ss_ratio", dest='ss_ratio', default=0.4)
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
    ss_ratio = float(args.ss_ratio)
    override_permutations = args.override_permutations.lower() == "true"
    base_folder = args.base_folder
    df = pd.DataFrame()
    summary = []
    for dataset in datasets:

        score_method = constants.PREDEFINED_SCORE
        if prefix == "GE":
            score_method = constants.DEG_EDGER
            if dataset.startswith("IE"):
                score_method = constants.DEG_T

        df_all_terms = pd.DataFrame()
        cur_real_ds = "{}_{}".format(prefix, dataset)

        p_means = []
        p_stds = []
        r_means = []
        r_stds = []
        f1_means = []
        f1_stds = []
        for algo in algos:
            recalls = []
            recalls = multiprocessing.Manager().list()
            precisions = multiprocessing.Manager().list()
            prcs = []
            p = MyPool(parallelization_factor)

            # for cur in np.arange(n_start,n_end):
            #     recovered_dataset_name = "{}_{}_recovery_{}_{}".format(prefix, dataset, algo, cur)
            #     if os.path.exists(os.path.join(constants.DATASETS_DIR, recovered_dataset_name)):
            #         shutil.rmtree(os.path.join(constants.DATASETS_DIR, recovered_dataset_name))
            #     if os.path.exists(os.path.join(constants.OUTPUT_GLOBAL_DIR, recovered_dataset_name)):
            #         shutil.rmtree(os.path.join(constants.OUTPUT_GLOBAL_DIR, recovered_dataset_name))

            params = [[recovery_iteration,
                       [prefix, dataset, x, algo, network_file_name, base_folder, ss_ratio, precisions, recalls, override_permutations]] for x
                      in np.arange(int(n_start), int(n_end))]
            p.map(func_star, params)
            p.close()

            print "done recovery iteration for dataset {} and algo {}".format(
                dataset, algo)
