
import sys
sys.path.insert(0, '../')

import pandas as pd
import constants
import os

import argparse


def main(datasets, algos):

    df_summary=pd.DataFrame()
    for cur_ds in datasets:

        constants.update_dirs(DATASET_NAME_u=cur_ds)
        total_num_genes=[]
        avg_num_genes=[]
        std_num_genes=[]
        n_enriched_terms=[]
        n_modules=[]

        for i_algo, cur_algo in enumerate(algos):
            print "current aggregation: {}, {}".format(cur_ds,cur_algo)
            try:
                n_genes_per_modules=pd.read_csv(
                    os.path.join(constants.TRUE_SOLUTIONS_DIR, "{}_{}".format(cur_ds,cur_algo), "report", "modules_summary.tsv"),
                    sep="\t")["#_genes"]
                n_modules.append(len(n_genes_per_modules.index))
                total_num_genes.append(n_genes_per_modules.sum())
                avg_num_genes.append(n_genes_per_modules.mean())
                std_num_genes.append(n_genes_per_modules.std())
                n_enriched_terms.append(len(pd.read_csv(
                    os.path.join(constants.OUTPUT_GLOBAL_DIR, "oob", "emp_diff_modules_{}_{}_passed_oob.tsv".format(cur_ds,cur_algo)),
                    sep="\t").dropna().index))

                df_series=pd.Series({"algo": cur_algo, "dataset": cur_ds, "sig_terms": n_enriched_terms[-1],
                                     "n_genes": total_num_genes[-1], "module_size_mean": avg_num_genes[-1], "module_size_std": std_num_genes[-1], "n_modules": n_modules[-1]})
                df_series.name = "{}_{}".format(cur_ds, cur_algo)
                df_summary=df_summary.append(df_series)

            except Exception,e:
                print "no genes were found for: {}, {}".format(cur_ds, cur_algo)
                total_num_genes.append(0)
                avg_num_genes.append(0)
                std_num_genes.append(0)

    return df_summary


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="tnfa,hc,ror,shera,shezh,ers,ien,apo,cbx,ift")
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="DOMINO,jactivemodules_greedy,netbox")

    args = parser.parse_args()

    prefix = args.prefix
    datasets=args.datasets.split(",")
    algos = args.algos.split(",")

    datasets = ["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"]
    algos = ["DOMINO", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    prefix = "GE"
    ds_summary=main(datasets=datasets, algos=algos)
    ds_summary.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation", "summary_statistics_{}.tsv".format(prefix)), sep='\t')

    datasets = ["brca", "crh", "scz", "tri", "t2d", "cad", "bmd", "hgt", "amd", "af"]
    algos = ["DOMINO", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    prefix = "PASCAL_SUM"
    ds_summary=main(datasets=datasets, algos=algos)
    ds_summary.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation", "summary_statistics_{}.tsv".format(prefix)), sep='\t')


