import sys
sys.path.insert(0, '../')

import numpy as np
import os
import pandas as pd

from infra import *
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from rpy2.robjects import pandas2ri
pandas2ri.activate()

import constants
import argparse
import simplejson as json

import scipy.stats as ss

from utils.go_similarity import calc_intra_similarity


def compute_redundancy(datasets, algos, pf=10, base_folder='/home/hag007/Desktop/aggregate_report/oob',
         file_format="emp_diff_modules_{}_{}_passed_oob.tsv", sim_method='Resnik', cutoffs=[0.0,1.0,2.0,3.0,4.0,5.0], ax=None): # [2,3,4,5,6,7]
    if not os.path.exists(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "ds_2_alg_scores")):
        try:
            os.makedirs(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "ds_2_alg_scores"))
        except Exception, e:
            print "error while creating ds_2_alg_scores folder: {}".format(e)

    results={}
    for cur_algo in algos:
        results[cur_algo]=[]

    for cutoff in cutoffs:
        df_summary=pd.DataFrame()
        for cur_ds in datasets:
            df = pd.DataFrame()
            print "cur ds: {}".format(cur_ds)
            constants.update_dirs(DATASET_NAME_u=cur_ds)
            algo_go_sim_score = []
            total_num_genes = []
            algos_signals = []

            for i_algo, cur_algo in enumerate(algos):
                reduced_list = []
                print "current cur_algo: {}".format(cur_algo)
                try:

                    emp_results = pd.read_csv(
                        os.path.join(base_folder,
                                     file_format.format(cur_ds, cur_algo)), sep='\t', index_col=0)

                except:
                    print "could not find {}".format(os.path.join(base_folder,
                                                                  file_format.format(cur_ds, cur_algo)), cur_ds, cur_algo)
                    total_num_genes.append(0)
                    algos_signals.append(0)
                    algo_go_sim_score.append(1)
                    continue

                emp_results = emp_results.sort_values(by='emp_rank')
                emp_results_fdr = emp_results.dropna().loc[emp_results.dropna()["passed_oob_permutation_test"].apply(
                    lambda a: np.any(np.array(a[1:-1].split(", ")) == "True")).values, :]

                algos_signals.append(len(emp_results_fdr.index))
                all_go_terms = emp_results_fdr.index.values

                try:
                    total_num_genes.append(pd.read_csv(
                        os.path.join(constants.OUTPUT_GLOBAL_DIR, "{}_{}".format(prefix, cur_ds), cur_algo,
                                     "all_modules_general.tsv"),
                        sep="\t")["total_num_genes"][0])
                except:
                    total_num_genes.append(0)

                cache_file = os.path.join(constants.CACHE_GLOBAL_DIR,
                                          "similarity_cache_emp_{}_{}_{}.npy".format(cur_ds, cur_algo, sim_method))


                all_go_terms_r, all_go_terms_o, adj = calc_intra_similarity(all_go_terms, pf, emp_results_fdr, cache_file,
                                                                            sim_method, cutoff)

                print len(all_go_terms_o), len(all_go_terms_r)
                if len(all_go_terms_r) > -1:
                    df.loc[cur_algo, "n_reduced_terms"]=len(all_go_terms_r)
                    df_summary.loc[cur_algo, cur_ds]=len(all_go_terms_r)

            df["ranked_terms"]=df.shape[0]-ss.rankdata(df.loc[:, "n_reduced_terms"])
            for k, v in df.iterrows():
                results[k].append(v["n_reduced_terms"])

        df_summary.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "solution_richness_matrix_{}_{}.tsv".format(prefix,cutoff)),sep='\t')

    for k,v in results.iteritems():
        if len(v)>0:
            print k, np.mean(v), v
        ax.plot(cutoffs,[ np.mean([v[i_ds+i_cutoff*len(datasets)] for i_ds in np.arange(len(datasets))])  for i_cutoff in np.arange(len(cutoffs))],label=k)

    ax.set_xlabel("similarity cutoff", fontdict={"size": 22})
    ax.set_ylabel("average non-redundant terms", fontdict={"size": 22})
    ax.legend(fontsize=22)


def main(datasets, algos, pf=10, base_folder='/home/hag007/Desktop/aggregate_report/oob', file_format="emp_diff_modules_{}_{}_passed_oob.tsv", sim_method='Resnik'):


    if not os.path.exists(os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr", "ds_2_alg_scores")):
        try:
            os.makedirs(os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr", "ds_2_alg_scores"))
        except Exception, e:
            print "error while creating ds_2_alg_scores folder: {}".format(e)

    for cur_ds in datasets:
        print "cur ds: {}".format(cur_ds) 
        constants.update_dirs(DATASET_NAME_u=cur_ds)
        algo_go_sim_score = []
        total_num_genes = []
        algos_signals = []

        for i_algo, cur_algo in enumerate(algos):
            print "current cur_algo: {}".format(cur_algo)
            try:

                emp_results = pd.read_csv(
                    os.path.join(base_folder,
                                 file_format.format(cur_ds, cur_algo)), sep='\t', index_col=0)

            except:
                print "could not find {}".format(os.path.join(base_folder,
                                                              file_format.format(cur_ds, cur_algo)), cur_ds, cur_algo)
                total_num_genes.append(0)
                algos_signals.append(0)
                algo_go_sim_score.append(1)
                continue

            emp_results=emp_results.sort_values(by='emp_rank')
            emp_results_fdr=emp_results.dropna().loc[emp_results.dropna()["passed_oob_permutation_test"].apply(lambda a: np.any(np.array(a[1:-1].split(", ")) == "True")).values,:]
            

            algos_signals.append(len(emp_results_fdr.index))
            all_go_terms = emp_results_fdr.index.values

            try:
                total_num_genes.append(pd.read_csv(
                    os.path.join(constants.OUTPUT_GLOBAL_DIR, "{}_{}".format(prefix,cur_ds), cur_algo, "all_modules_general.tsv"),
                    sep="\t")["total_num_genes"][0])
            except:
                total_num_genes.append(0)

            cache_file=os.path.join(constants.CACHE_GLOBAL_DIR, "similarity_cache_emp_{}_{}_{}.npy".format(cur_ds,cur_algo, sim_method))
            all_go_terms_r, all_go_terms_o, adj = calc_intra_similarity(all_go_terms, pf, emp_results_fdr, cache_file, sim_method)

            emp_results_fdr[emp_results_fdr.index.isin(all_go_terms_r)].to_csv(os.path.join(base_folder, file_format.format(cur_ds,cur_algo)[:-4]+"_reduced.tsv"), sep='\t')
            open(cache_file, 'w+').write(json.dumps(dict(adj)))

            print "total redundancy removal: {}/{}".format(emp_results_fdr[emp_results_fdr.index.isin(all_go_terms_r)].shape[0], emp_results_fdr.shape[0])


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="TNFa_2,HC12,ROR_1,ERS_1,IEM,SHERA,SHEZH_1") # ",Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50 TNFa_2,HC12,SHERA,ROR_1,SHEZH_1,ERS_1,IEM"
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="dcem,jactivemodules_greedy,jactivemodules_sa,bionet,netbox,my_netbox_td,keypathwayminer_INES_GREEDY,hotnet2") # jactivemodules_greedy,jactivemodules_sa,bionet,netbox,my_netbox_td
    parser.add_argument('--pf', dest='pf', default=3)
    parser.add_argument('--base_folder', dest='base_folder', default=os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr","MAX"))
    parser.add_argument('--file_format', dest='file_format', default="emp_diff_modules_{}_{}_passed_oob.tsv")
    parser.add_argument('--sim_method', dest='sim_method', default="Resnik")
    args = parser.parse_args()

    prefix = args.prefix
    file_format = args.file_format
    base_folder= args.base_folder
    sim_method = args.sim_method
    datasets=["{}".format(x) for x in args.datasets.split(",")]
    algos = args.algos.split(",")
    cutoffs=[1,2,3,4,5]
    pf=int(args.pf)


    # fig, axs = plt.subplots(1,2,figsize=(20,10))
    # prefix = "GE"
    # datasets=["TNFa_2" ,"HC12","ROR_1","ERS_1","IEM","SHERA","SHEZH_1"]
    # algos=["dcem", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox", "keypathwayminer_INES_GREEDY", "hotnet2"]
    # compute_redundancy(datasets=datasets, algos=algos, pf=pf, base_folder=base_folder, file_format=file_format, sim_method=sim_method, cutoffs=cutoffs, ax=axs[0])

    prefix = "PASCAL_SUM"
    datasets=["Breast_Cancer.G50", "Crohns_Disease.G50", "Schizophrenia.G50", "Triglycerides.G50", "Type_2_Diabetes.G50"]
    algos = ["my_netbox_td"] # "dcem","jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox", "my_netbox_td"
    compute_redundancy(datasets=datasets, algos=algos, pf=pf, base_folder=base_folder, file_format=file_format, sim_method=sim_method, cutoffs=cutoffs, ax=axs[1])

    # fig.text(0.01, 0.97, "A:", weight='bold', fontsize=22)
    # fig.text(0.5, 0.97, "B:", weight='bold', fontsize=22)
    # fig.tight_layout()
    # fig.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,"figure_13.png"))
