import sys
sys.path.insert(0, '../')

import pandas as pd

from fastsemsim.SemSim import *


import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as ml_colors
from matplotlib.lines import Line2D

from rpy2.robjects import pandas2ri
pandas2ri.activate()

import constants

import argparse

import simplejson as json

import seaborn as sns

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="TNFa_2,HC12,ROR_1,SHEZH_1,ERS_1,IEM") #  SHERA "Breast_Cancer.G50,Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50 TNFa_2,HC12,SHERA,ROR_1,SHEZH_1,ERS_1,IEM"
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy,jactivemodules_sa,bionet,netbox,my_netbox_td,hotnet2,keypathwayminer_INES_GREEDY") #
    # parser.add_argument('--module_indices', dest='module_indices',
    #                     default="0,1,2")  # jactivemodules_greedy,jactivemodules_sa,bionet,netbox,my_netbox_td
    parser.add_argument('--pf', dest='pf', default=3)
    parser.add_argument('--base_folder', dest='base_folder', default=os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr","MAX"))
    parser.add_argument('--sim_method', dest='sim_method', default="Resnik")
    parser.add_argument('--file_format', dest='file_format', default="emp_diff_modules_{}_{}_passed_oob.tsv")
    parser.add_argument('--cutoffs', dest='cutoffs', default="8,9,10,11")

    args = parser.parse_args()

    prefix = args.prefix
    base_folder = args.base_folder
    sim_method= args.sim_method
    file_format=args.file_format
    # semsim = SetSemSim(ontology, ac, TSS=sim_method, MSS="BMA")
    datasets=["{}".format(x) for x in args.datasets.split(",")]
    # module_indices=args.module_indices.split(",")
    algos = args.algos.split(",")
    cutoffs = np.array(args.cutoffs.split(','),dtype=float)
    pf=int(args.pf)
    print "test"
    h_scores = pd.DataFrame()
    df_full_data = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "modules_full_data_{}.tsv".format(prefix)), index_col=0, sep='\t')
    agg_dict={}
    for cur_ds in datasets:
        for cur_alg in algos:
            module_indices=[a.split("_")[-1] for a in df_full_data.loc[(df_full_data['algo']== cur_alg).values & (df_full_data['dataset']== cur_ds).values & (df_full_data['EHR']>0.2).values].index.values]
            for cur_module_index_0 in module_indices:
                set_0=df_full_data.loc["{}_{}_module_{}".format(cur_ds,cur_alg,cur_module_index_0), "tp"]
                if not type(set_0) is str:
                    set_0=[]
                else:
                    set_0=[a.split(": ")[0] for a in df_full_data.loc["{}_{}_module_{}".format(cur_ds,cur_alg,cur_module_index_0), "tp"].split("\n")]
                for cur_module_index_1 in module_indices:

                    if cur_module_index_0<cur_module_index_1:
                        cache_file=os.path.join(constants.CACHE_GLOBAL_DIR, "similarity_cache_{}_{}_{}_{}_{}.npy".format(sim_method, cur_ds,cur_alg, cur_module_index_0, cur_module_index_1))
                    else:
                        cache_file=os.path.join(constants.CACHE_GLOBAL_DIR, "similarity_cache_{}_{}_{}_{}_{}.npy".format(sim_method, cur_ds, cur_alg, cur_module_index_1, cur_module_index_0))

                    agg_dict.update(json.load(open(cache_file, 'r')))

    plt.yscale('log')
    print "maximal similarity: {}".format(agg_dict.values())
    sns.distplot(filter(lambda a: a>=0, agg_dict.values()), norm_hist=False, kde=False)
    plt.legend()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "similarity_values_dist.png"))
