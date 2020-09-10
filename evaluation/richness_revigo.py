import sys
sys.path.insert(0, '../')

import pandas as pd

from infra import *

import constants
import argparse

from utils.go_similarity import calc_intra_similarity
from utils.daemon_multiprocessing import MyPool, func_star

def compute_redundancy_for_solution(cutoff, cur_ds, cur_algo, base_folder, file_format, sim_method, pf):
    print "current cur_algo: {}".format(cur_algo)
    try:
        emp_results = pd.read_csv(os.path.join(base_folder, file_format.format(cur_ds, cur_algo)), sep='\t', index_col=0)

    except:
        print "could not find {}".format(os.path.join(base_folder, file_format.format(cur_ds, cur_algo)), cur_ds, cur_algo)
        return cutoff, cur_ds, cur_algo, 0

    emp_results = emp_results.sort_values(by='emp_pval_max')
    emp_results_fdr = emp_results.dropna().loc[emp_results.dropna()["passed_oob_permutation_test"].apply(
        lambda a: np.any(np.array(a[1:-1].split(", ")) == "True")).values, :]


    cache_file = os.path.join(constants.CACHE_GLOBAL_DIR, "similarity_cache_emp_{}_{}_{}.npy".format(cur_ds, cur_algo, sim_method))
    all_go_terms_r, all_go_terms_o, adj = calc_intra_similarity(emp_results_fdr.index.values, pf, emp_results_fdr, cache_file,
                                                                sim_method, cutoff)
    print "# original terms: ", len(all_go_terms_o), "# reduced terms: ", len(all_go_terms_r)
    if len(all_go_terms_r) > -1:
        richness=len(all_go_terms_r)
    else:
        richness=0

    return cutoff, cur_ds, cur_algo, richness

def compute_redundancy_for_solutions(prefix, datasets, algos, pf=10, sim_method='Resnik', cutoffs=[1.0, 2.0, 3.0, 4.0], base_folder=None, file_format=None):

    output = {c: pd.DataFrame() for c in cutoffs}

    params=[]
    for cutoff in cutoffs:
        print "cut cutoff: {}".format(cutoff)
        for cur_ds in datasets:
            print "cur ds: {}".format(cur_ds)
            constants.update_dirs(DATASET_NAME_u=cur_ds)
            for i_algo, cur_algo in enumerate(algos):
                params.append([compute_redundancy_for_solution,[cutoff, cur_ds, cur_algo, base_folder, file_format, sim_method, pf]])

    p = MyPool(pf)
    results=p.map(func_star, params)
    p.close()
    for cur_res in results:
        output[cur_res[0]].loc[cur_res[2], cur_res[1]]=cur_res[3]

    for cutoff in cutoffs:
        output[cutoff].to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation", "richness_matrix_{}_{}.tsv".format(prefix,cutoff)),sep='\t')

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

    compute_redundancy_for_solutions(prefix=prefix, datasets=datasets, algos=algos, pf=pf, base_folder=base_folder, file_format=file_format, sim_method=sim_method, cutoffs=cutoffs)

