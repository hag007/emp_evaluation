import sys

sys.path.insert(0, '../')

from infra import *
import pandas as pd
import numpy as np
import argparse
from pandas.errors import EmptyDataError
from utils.daemon_multiprocessing import MyPool, func_star
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
from utils.add_GO_terms_metadata_agg import get_all_genes_for_term, vertices
import multiprocessing



def get_enriched_terms(algos, datasets):
    for cur_algo in algos:
        algos_filter = cur_algo

        df_go = pd.DataFrame(columns=['GO id', 'qval', 'pval'])
        df_go_pvals = pd.DataFrame()
        df_go_pvals.index.name = "GO id"
        for cur_ds in datasets:
            go_results = [os.path.join(constants.ROBUSTNESS_SOLUTIONS_DIR, cur_ds, "report", a) for a in os.listdir(os.path.join(constants.ROBUSTNESS_SOLUTIONS_DIR, cur_ds, "report"))
                          if "separated_modules" in a]

            for cur in go_results:
                try:
                    df_go = pd.concat((df_go, pd.read_csv(cur, sep='\t')))
                    df_go_pvals = pd.concat((df_go_pvals, pd.read_csv(cur, sep='\t').set_index("GO id")['pval']),
                                            axis=1)
                except (EmptyDataError, KeyError):
                    pass

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


def aggregate_solutions(dataset, cur, algo,
                        base_folder=None, ss_ratio=0.4, precisions=None,
                        recalls=None):
    print "starting iteration: {}, {}".format(dataset, cur)
    try:
        recovered_dataset_name = "sol_{}_{}_robustness_{}_{}".format(algo, dataset, cur, ss_ratio)

        cur_pval, df_terms, df_pval_terms = get_enriched_terms([algo], [recovered_dataset_name])

        df = pd.read_csv(os.path.join(base_folder, "oob", "emp_diff_modules_{}_{}_passed_oob.tsv".format(dataset, algo)),
                         sep='\t')

        full_sig_terms = df.loc[df["passed_oob_permutation_test"].dropna(axis=0).apply(
            lambda a: np.any(np.array(a[1:-1].split(", ")) == "True")).values, :].sort_values(by=["hg_pval_max"],
                                                                                              ascending=False)['GO id']
    except KeyError, e:
        print e
        full_sig_terms = np.array([])

    df_detailed_pr = df_pval_terms.to_frame().rename(columns= {0: 'pval'})
    df_detailed_pr.columns
    df_detailed_pr['is_significant'] = df_detailed_pr.index.isin(full_sig_terms)

    df_detailed_pr.to_csv(os.path.join(constants.ROBUSTNESS_SOLUTIONS_DIR, recovered_dataset_name, "df_detailed_pr.tsv"),
                         sep='\t')


    recall = len(set(full_sig_terms).intersection(df_pval_terms.index.values)) / float(max(1, len(full_sig_terms)))
    precision = len(set(full_sig_terms).intersection(df_pval_terms.index.values)) / float(
        max(1, df_pval_terms.shape[0]))
    print "precision: {}/{}={} ({}, {})".format(len(set(full_sig_terms).intersection(df_pval_terms.index.values)),
                                                float(df_pval_terms.shape[0]), precision, len(set(full_sig_terms)),
                                                len(set(df_pval_terms.index.values)))
    print "recall: {}/{}={} ({}, {})".format(len(set(full_sig_terms).intersection(df_pval_terms.index.values)),
                                             float(len(full_sig_terms)), recall, len(set(full_sig_terms)),
                                             len(set(df_pval_terms.index.values)))
    print "done iteration: {}, {}".format(dataset, cur)
    # if precisions is not None:
    precisions.append(precision)
    recalls.append(recall)
    return precision, recall


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--prefix', dest='prefix', default="GE") 
    parser.add_argument('--datasets', dest='datasets', default="brca")
    parser.add_argument('--algos', dest='algos', default="DOMINO")
    parser.add_argument('--network', dest='network', default="dip.sif")
    parser.add_argument('--pf', help="parallelization_factor", dest='pf', default=10)
    parser.add_argument('--n_start_r', help="number of iterations (total n permutation is pf*(n_end-n_start))",
                        dest='n_start_r', default=0)
    parser.add_argument('--n_end_r', help="number of iterations (total n permutation is pf*(n_end-n_start))",
                        dest='n_end_r', default=100)
    parser.add_argument('--ss_ratios', help="ss_ratios", dest='ss_ratios', default="0.4,0.3,0.2,0.1") # ,0.2,0.2,0.1
    parser.add_argument('--override_permutations', help="takes max or all samples", dest='override_permutations',
                        default="false")
    parser.add_argument('--base_folder', help="base_folder", dest='base_folder', default=constants.OUTPUT_GLOBAL_DIR)

    args = parser.parse_args()
    prefix=args.prefix
    datasets = args.datasets.split(",")
    algos = args.algos.split(",")
    network_file_name = args.network
    parallelization_factor = int(args.pf)
    n_start = int(args.n_start_r)
    n_end = int(args.n_end_r)
    ss_ratios = [float(a) for a in args.ss_ratios.split(",")]
    override_permutations = args.override_permutations.lower() == "true"
    base_folder = args.base_folder
    df = pd.DataFrame()
    summary = []

    for ss_ratio in ss_ratios:
        df = pd.DataFrame()
        summary = []
        for dataset in datasets:
            
            df_all_terms = pd.DataFrame()
            cur_real_ds = dataset

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

                params = [[aggregate_solutions,
                           [dataset, x, algo, base_folder, ss_ratio, precisions, recalls]] for x
                          in np.arange(int(n_start), int(n_end))]
                p.map(func_star, params)
                p.close()

                # [aggregate_solutions(dataset, x, algo, base_folder, ss_ratio, precisions, recalls, override_permutations) for x
                #           in np.arange(int(n_start), int(n_end))]

                f1s = [(2 * p * r) / float(max(p + r, 10e-5)) for p, r in zip(precisions, recalls)]
                r_means.append(np.mean(recalls))
                r_stds.append(np.std(recalls))
                p_means.append(np.mean(precisions))
                p_stds.append(np.std(precisions))
                f1_means.append(np.mean(f1s))
                f1_stds.append(np.std(f1s))
                print "recovery_fractions for dataset {} and algo {}:\nnp_mean:{}, p_std: {}\nr_means: {}, r_stds: {}".format(
                    dataset, algo, p_means[-1], p_stds[-1], r_means[-1], r_stds[-1])
                df = df.append({"p_mean": p_means[-1], "p_std": p_stds[-1],
                                "r_mean": r_means[-1], "r_std": r_stds[-1], "f1_mean": f1_means[-1], "f1_std": f1_stds[-1],
                                "algo": algo, "dataset": dataset, "ss_ratio": ss_ratio}, ignore_index=True)

                df.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation",
                                       "robustness_f1_{}_{}_{}.tsv".format(prefix, n_end, ss_ratio)), sep='\t')
                print "save file to: {}".format(os.path.join(constants.OUTPUT_GLOBAL_DIR, "robustness_results_{}_{}_{}.tsv".format(prefix, n_end, ss_ratio)))

