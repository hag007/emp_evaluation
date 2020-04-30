import sys
sys.path.insert(0, '../')

import pandas as pd

from infra import *

import constants
import argparse
import simplejson as json

import scipy.stats as ss

from utils.go_similarity import calc_intra_similarity


def compute_redundancy(datasets, algos, pf=10, sim_method='Resnik', cutoffs=[1.0,2.0,3.0,4.0], base_folder=None, file_format=None):
    if not os.path.exists(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "ds_2_alg_scores")):
        try:
            os.makedirs(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "ds_2_alg_scores"))
        except Exception, e:
            print "error while creating ds_2_alg_scores folder: {}".format(e)


    for cutoff in cutoffs:
        print "cut cutoff: {}".format(cutoff)
        df_summary=pd.DataFrame()
        for cur_ds in datasets:
            df = pd.DataFrame()
            print "cur ds: {}".format(cur_ds)
            constants.update_dirs(DATASET_NAME_u=cur_ds)
            algo_go_sim_score = []
            total_num_genes = []
            algos_signals = []

            for i_algo, cur_algo in enumerate(algos):
                print "current cur_algo: {}".format(cur_algo)
                try:
                    emp_results = pd.read_csv(os.path.join(base_folder, file_format.format(cur_ds, cur_algo)), sep='\t', index_col=0)

                except:
                    print "could not find {}".format(os.path.join(base_folder, file_format.format(cur_ds, cur_algo)), cur_ds, cur_algo)
                    algos_signals.append(0)
                    algo_go_sim_score.append(1)
                    continue

                emp_results = emp_results.sort_values(by='emp_pval_max')
                emp_results_fdr = emp_results.dropna().loc[emp_results.dropna()["passed_oob_permutation_test"].apply(
                    lambda a: np.any(np.array(a[1:-1].split(", ")) == "True")).values, :]

                algos_signals.append(len(emp_results_fdr.index))
                all_go_terms = emp_results_fdr.index.values

                try:
                    total_num_genes.append(pd.read_csv(
                        os.path.join(constants.TRUE_SOLUTIONS_DIR, "{}_{}".format(cur_ds, cur_algo), "report", "modules_summary.tsv"),
                        sep="\t")["#_genes"].sum())
                except:
                    total_num_genes.append(0)

                cache_file = os.path.join(constants.CACHE_GLOBAL_DIR, "similarity_cache_emp_{}_{}_{}.npy".format(cur_ds, cur_algo, sim_method))


                all_go_terms_r, all_go_terms_o, adj = calc_intra_similarity(all_go_terms, pf, emp_results_fdr, cache_file,
                                                                            sim_method, cutoff)

                print "# original terms: ", len(all_go_terms_o), "# reduced terms: ", len(all_go_terms_r)
                if len(all_go_terms_r) > -1:
                    df.loc[cur_algo, "n_reduced_terms"]=len(all_go_terms_r)
                    df_summary.loc[cur_algo, cur_ds]=len(all_go_terms_r)

            df["ranked_terms"]=df.shape[0]-ss.rankdata(df.loc[:, "n_reduced_terms"])


        df_summary.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation", "richness_matrix_{}_{}.tsv".format(prefix,cutoff)),sep='\t')





def main(datasets, algos, pf=10, base_folder=None, file_format=None, sim_method='Resnik'):


    # if not os.path.exists(os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr", "ds_2_alg_scores")):
    #     try:
    #         os.makedirs(os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr", "ds_2_alg_scores"))
    #     except Exception, e:
    #         print "error while creating ds_2_alg_scores folder: {}".format(e)

    for cur_ds in datasets:
        print "cur ds: {}".format(cur_ds) 
        constants.update_dirs(DATASET_NAME_u=cur_ds)
        algo_go_sim_score = []
        total_num_genes = []
        algos_signals = []

        for i_algo, cur_algo in enumerate(algos):
            print "current cur_algo: {}".format(cur_algo)
            try:
                emp_results = pd.read_csv(os.path.join(base_folder, file_format.format(cur_ds, cur_algo)), sep='\t', index_col=0)

            except:
                print "could not find {}".format(os.path.join(base_folder, file_format.format(cur_ds, cur_algo)), cur_ds, cur_algo)
                total_num_genes.append(0)
                algos_signals.append(0)
                algo_go_sim_score.append(1)
                continue

            emp_results=emp_results.sort_values(by='emp_rank')
            emp_results_fdr=emp_results.dropna().loc[emp_results.dropna()["passed_oob_permutation_test"].apply(lambda a: np.any(np.array(a[1:-1].split(", ")) == "True")).values,:]

            algos_signals.append(len(emp_results_fdr.index))
            all_go_terms = emp_results_fdr.index.values

            try:
                total_num_genes.append(pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "{}_{}".format(cur_ds,cur_algo), cur_algo, "modules_summary.tsv"), sep="\t")["#_genes"].sum())
            except:
                total_num_genes.append(0)

            cache_file=os.path.join(constants.CACHE_GLOBAL_DIR, "similarity_cache_emp_{}_{}_{}.npy".format(cur_ds,cur_algo, sim_method))
            all_go_terms_r, all_go_terms_o, adj = calc_intra_similarity(all_go_terms, pf, emp_results_fdr, cache_file, sim_method)

            emp_results_fdr[emp_results_fdr.index.isin(all_go_terms_r)].to_csv(os.path.join(base_folder, file_format.format(cur_ds,cur_algo)[:-4]+"_reduced.tsv"), sep='\t')
            open(cache_file, 'w+').write(json.dumps(dict(adj)))

            print "total redundancy removal: {}/{}".format(emp_results_fdr[emp_results_fdr.index.isin(all_go_terms_r)].shape[0], emp_results_fdr.shape[0])


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="TNFa_2,HC12,ROR_1,ERS_1,IEM,SHERA,SHEZH_1,APO,CBX,IFT") # ",Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50 TNFa_2,HC12,SHERA,ROR_1,SHEZH_1,ERS_1,IEM"
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="dcem,jactivemodules_greedy,jactivemodules_sa,bionet,netbox,my_netbox_td,keypathwayminer_INES_GREEDY,hotnet2") # jactivemodules_greedy,jactivemodules_sa,bionet,netbox,my_netbox_td
    parser.add_argument('--pf', dest='pf', default=10)
    parser.add_argument('--base_folder', dest='base_folder', default=os.path.join(constants.OUTPUT_GLOBAL_DIR,"oob"))
    parser.add_argument('--file_format', dest='file_format', default="emp_diff_modules_{}_{}_passed_oob.tsv")
    parser.add_argument('--sim_method', dest='sim_method', default="Resnik")
    args = parser.parse_args()

    prefix = args.prefix
    file_format = args.file_format
    base_folder= args.base_folder
    sim_method = args.sim_method
    datasets=["{}".format(x) for x in args.datasets.split(",")]
    algos = args.algos.split(",")
    cutoffs = [1.0, 2.0, 3.0, 4.0]
    pf=int(args.pf)

    datasets=["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"]
    algos=["DOMINO", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    prefix="GE"
    compute_redundancy(datasets=datasets, algos=algos, pf=pf, base_folder=base_folder, file_format=file_format, sim_method=sim_method, cutoffs=cutoffs)

    datasets=["brca", "crh", "scz", "tri", "t2d", "cad", "bmd", "hgt", "amd", "af"]
    algos=["DOMINO", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    prefix="PASCAL_SUM"
    compute_redundancy(datasets=datasets, algos=algos, pf=pf, base_folder=base_folder, file_format=file_format, sim_method=sim_method, cutoffs=cutoffs)

