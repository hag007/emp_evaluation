import sys
sys.path.insert(0, '../')
import os
import constants

import pandas as pd
import numpy as np
import argparse

from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="brca")
    parser.add_argument('--algos', dest='algos', default="DOMINO")
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--n_start', help="number of iterations (total n permutation is pf*(n_end-n_start))",
                        dest='n_start', default=0)
    parser.add_argument('--n_end', help="number of iterations (total n permutation is pf*(n_end-n_start))",
                        dest='n_end', default=100)
    parser.add_argument('--ss_ratios', help="ss_ratios", dest='ss_ratios', default="0.4,0.3,0.2,0.1")

    args = parser.parse_args()

    datasets = args.datasets.split(",")
    algos = args.algos.split(",")
    n_start = int(args.n_start)
    n_end = int(args.n_end)
    prefix = args.prefix
    ss_ratios = [float(a) for a in args.ss_ratios.split(",")]

    for ss_ratio in ss_ratios:
        df_summary = pd.DataFrame()
        df_pr_auc=pd.DataFrame()
        summary=[]

        for dataset in datasets:

            p_means = []
            p_stds = []
            r_means = []
            r_stds = []
            f1_means = []
            f1_stds = []
            for algo in algos:
                df_detailed_pr_agg=pd.DataFrame()
                df_detailed_pr_is_sig_agg = pd.DataFrame()
                df_detailed_pr_pval_agg = pd.DataFrame()

                for cur in np.arange(int(n_start), int(n_end)):
                    recovered_dataset_name = "sol_{}_{}_robustness_{}_{}".format(algo, dataset, cur, ss_ratio)
                    try:
                        df_detailed_pr = pd.read_csv(
                            os.path.join(constants.ROBUSTNESS_SOLUTIONS_DIR, recovered_dataset_name, "df_detailed_pr.tsv"),
                            sep='\t', index_col=0)
                    except IOError,e:
                        print e
                        print "continue..."
                        continue

                    df_detailed_pr_pval_agg = pd.concat([df_detailed_pr_pval_agg, df_detailed_pr['pval']], axis=1)
                    df_detailed_pr_is_sig_agg = pd.concat([df_detailed_pr_is_sig_agg, df_detailed_pr['is_significant']], axis=1)

                try:
                    ehr_terms = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"oob", "emp_diff_modules_{}_{}_passed_oob.tsv".format(dataset, algo)), sep='\t')
                except Exception as e:
                    print("error: {}".format(e))
                    continue
                    
                try:
                    ehr_terms=ehr_terms.loc[ehr_terms["passed_oob_permutation_test"].dropna(axis=0).apply(
                lambda a: np.any(np.array(a[1:-1].split(", ")) == "True")).values, :].sort_values(by=["hg_pval_max"],  ascending=False)['GO id']
                except KeyError,e:
                    ehr_terms = pd.Series()

                print "intersected terms: {}/{}".format(len(set(ehr_terms).intersection(df_detailed_pr_pval_agg.index)), len(ehr_terms.index))
                df_detailed_pr_agg['pval_frequency']=df_detailed_pr_pval_agg.apply(lambda a: np.sum(~pd.isnull(a)), axis=1)
                df_detailed_pr_agg['is_sig']=df_detailed_pr_is_sig_agg.apply(lambda a: np.any(a), axis=1)

                missing_sig_terms=ehr_terms.loc[~ehr_terms.isin(df_detailed_pr_pval_agg.index)]
                print "n of missing terms: {}/{}".format(missing_sig_terms.shape[0], df_detailed_pr_pval_agg.shape[0])
                for cur_missing_term in missing_sig_terms.values:
                    df_detailed_pr_agg.loc[cur_missing_term, 'is_sig'] = False
                    df_detailed_pr_agg.loc[cur_missing_term, 'pval_frequency'] = 0
                df_detailed_pr_agg=df_detailed_pr_agg.sort_values(by=['pval_frequency'], ascending=False)


                y_test=df_detailed_pr_agg['is_sig'].values.astype(np.int)
                print "df shape: {}".format(df_detailed_pr_agg.shape)
                y_score=df_detailed_pr_agg['pval_frequency'].values/float(n_end-n_start)

                print "y_test (total={}):\n{}".format(np.sum(y_test), y_test)
                print "y_score:\n{}".format(y_score)

                df_pr_auc.loc[algo,dataset]=np.nan
                if len(y_score)!=0:
                    average_precision = average_precision_score(y_test,y_score)
                    if np.isnan(average_precision):
                        average_precision = 0
                    df_detailed_pr_agg.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "robustness_cache", "recovery_terms_frequency_{}_{}.tsv".format(n_end, ss_ratio)), sep='\t')

                    precision, recall, _ = precision_recall_curve(y_test, y_score)

                    print "average precision: {}".format(average_precision)
                    print "save curve in :{}".format(os.path.join(constants.OUTPUT_GLOBAL_DIR, "robustness_cache", "pr_curve_terms_recovery_{}_{}_{}_{}.png".format(dataset,algo, n_end, ss_ratio)))
                    df_pr_auc.loc[algo,dataset]=average_precision
        df_pr_auc.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation", "robustness_auc_{}_{}_{}.tsv".format(prefix,n_end,ss_ratio)),sep='\t')
