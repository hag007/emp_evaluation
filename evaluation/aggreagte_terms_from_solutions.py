import sys
sys.path.insert(0, '../')


import constants
import os
import argparse
import pandas as pd
import numpy as np

import matplotlib as mpl
font = {'size': 12}
mpl.rc('font', **font)
mpl.rc('xtick', labelsize=12)    # fontsize of the tick labels
mpl.rc('ytick', labelsize=12)
mpl.rc('axes', labelsize=12)
mpl.rc('legend', fontsize=12)

import matplotlib.pyplot as plt

import seaborn as sns

def aggregate_solution(algo_sample = None, dataset_sample = None,tsv_file_name=None, filtered_go_ids_file=None, hg_th=0.05):

    try:
        output_oob = pd.read_csv(
            tsv_file_name.format(dataset_sample, algo_sample),
            sep='\t', index_col=0).dropna()
    except Exception, e:
        print e
        return 0, 0

    filtered_go_ids=open(filtered_go_ids_file,'r').read().split()
    filtered_genes=output_oob.reindex(filtered_go_ids).loc[:,["GO name", "is_emp_pval_max_significant"]].dropna()
    return filtered_genes

def main(datasets, algos, prefix, hg_th=0.05):
    df_sum=pd.DataFrame()
    filtered_go_ids_file=os.path.join(constants.GO_DIR,"filtered_go_terms.txt")
    sig_in_null_file=os.path.join(constants.GO_DIR,"sig_in_null.txt")
    sig_in_null=open(sig_in_null_file,'r').read().split()
    filtered_go_ids=open(filtered_go_ids_file,'r').read().split()
    df_ev_matrix=pd.DataFrame()
    df_non_ev_matrix=pd.DataFrame()
    df_ratio_matrix=pd.DataFrame()
    for a in filtered_go_ids:
        df_sum.loc[a,"EV"]=0
        df_sum.loc[a,"non_EV"]=0
        df_sum.loc[a,"is_sig_in_radalib"]=a in sig_in_null
        for cur_alg in algos:
            df_ev_matrix.loc[a,cur_alg]=0
            df_non_ev_matrix.loc[a,cur_alg]=0
    for cur_ds in datasets:
        for cur_alg in algos:
            filtered_genes=aggregate_solution(cur_alg, cur_ds,tsv_file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR,"oob", "emp_diff_modules_{}_{}_passed_oob.tsv"), filtered_go_ids_file=filtered_go_ids_file,  hg_th=hg_th)
            for idx, row in  filtered_genes.iterrows():
                df_sum.loc[idx,"EV"]+=int(row["is_emp_pval_max_significant"])
                df_sum.loc[idx,"non_EV"]+=int(not row["is_emp_pval_max_significant"])
                df_ev_matrix.loc[idx,cur_alg]+=int(row["is_emp_pval_max_significant"])
                df_non_ev_matrix.loc[idx,cur_alg]+=int(not row["is_emp_pval_max_significant"])

    for a in filtered_go_ids:
        if float(df_sum.loc[a,"non_EV"]+df_sum.loc[a,"EV"]) == 1:
            df_sum.loc[a,"ratio"]=df_sum.loc[a,"non_EV"]/float(df_sum.loc[a,"non_EV"]+df_sum.loc[a,"EV"])
        else:
            df_sum.loc[a,"ratio"]=np.nan

        for cur_alg in algos:
            if df_non_ev_matrix.loc[a,cur_alg]+df_ev_matrix.loc[a,cur_alg] == 1:
                df_ratio_matrix.loc[a,cur_alg]=df_non_ev_matrix.loc[a,cur_alg]/float(df_non_ev_matrix.loc[a,cur_alg]+df_ev_matrix.loc[a,cur_alg])
            else:
                df_ratio_matrix.loc[a,cur_alg]=np.nan

    df_ratio_matrix.loc[:,"is_sig_in_radalib"]=False
    df_ratio_matrix.loc[sig_in_null,"is_sig_in_radalib"]=True

    sns.distplot(df_sum.loc[df_sum.loc[:,"is_sig_in_radalib"]].loc[:,"ratio"].dropna(), kde=False)
    print(df_sum.loc[df_sum.loc[:,"is_sig_in_radalib"]].loc[:,"ratio"].mean())
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "ev_dist", "sig_ratio_dist_{}.png".format(prefix)))
    plt.clf()
    sns.distplot(df_sum.loc[~df_sum.loc[:,"is_sig_in_radalib"]].loc[:,"ratio"].dropna(), kde=False)
    print(df_sum.loc[~df_sum.loc[:,"is_sig_in_radalib"]].loc[:,"ratio"].mean())
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "ev_dist", "non_sig_ratio_dist_{}.png".format(prefix)))
    plt.clf()
    sns.distplot(df_sum.loc[df_sum.loc[:,"is_sig_in_radalib"]].loc[:,"EV"].dropna(), kde=False)
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "ev_dist","sig_ev_dist_{}.png".format(prefix)))
    plt.clf()
    sns.distplot(df_sum.loc[~df_sum.loc[:,"is_sig_in_radalib"]].loc[:,"EV"].dropna(), kde=False)
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "ev_dist","non_sig_ev_dist_{}.png".format(prefix)))
    plt.clf()
    sns.distplot(df_sum.loc[df_sum.loc[:,"is_sig_in_radalib"]].loc[:,"non_EV"].dropna(), kde=False)
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "ev_dist","sig_non_ev_dist_{}.png".format(prefix)))
    plt.clf()
    sns.distplot(df_sum.loc[~df_sum.loc[:,"is_sig_in_radalib"]].loc[:,"non_EV"].dropna(), kde=False)
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "ev_dist","non_sig_non_ev_dist_{}.png".format(prefix)))
    plt.clf()

    df_sum.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation", "aggregate_solutions_{}.tsv".format(prefix)),sep='\t')
    df_ratio_matrix.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation", "aggregate_solutions_ratio_by_algo_{}.tsv".format(prefix)),sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="TNFa_2,HC12,ROR_1,SHERA,SHEZH_1,ERS_1,IEM,APO,CBX,IFT")
    parser.add_argument('--prefix', dest='prefix', default="GE") # PASCAL_SUM   GE
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy,jactivemodules_sa,bionet,netbox,keypathwayminer_INES_GREEDY,DOMINO,DOMINO2,hotnet2")
    args = parser.parse_args()
    datasets = args.datasets.split(",")
    algos = args.algos.split(",")
    prefix = args.prefix

    datasets=['tnfa', 'hc', 'ror', 'cbx', 'shera' , 'shezh' , 'ift', 'iem' , "ers", "apo"]
    algos=["DOMINO2", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    prefix="GE"
    main(datasets,algos,prefix)

    datasets=['brca', 'crh', 'scz', 'tri', 't2d', 'af', 'cad', 'amd', 'hgt' , "bmd"]
    algos=["DOMINO2", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    prefix="PASCAL_SUM"
    main(datasets,algos,prefix)
# ["DOMINO4_per_01_1000_huri", "DOMINO4_421_018_5000_huri", "DOMINO4_45205_01944_1000_huri", "DOMINO4_per_013_1000_huri", "DOMINO4_per_01_1000_huri", "DOMINO_test6_max_07_huri", "DOMINO_test6_max_08_huri", "DOMINO_test6_max_09_huri","DOMINO_test6_max_08_01_huri","DOMINO4_max_08_015_2_huri","DOMINO2_max_08_012_2_dip"] # ,"DOMINO2_max_08_dip","DOMINO2_max_08_013_dip"] # ["DOMINO4_per_01_1000_huri", "DOMINO4_421_018_5000_huri", "DOMINO4_45205_01944_1000_huri", "DOMINO4_per_013_1000_huri", "DOMINO4_per_01_1000_huri","DOMINO_test6_max_09_huri", "DOMINO_test6_max_07_huri"] # ["DOMINO", "netbox", "DOMINO3_dip", "DOMINO3_1_dip", "DOMINO3_421_dip","DOMINO_test4_max_07_dip", "DOMINO_test2_max_08_dip", "DOMINO_test3_max_09_dip", "DOMINO_test5_max_08_string3", "DOMINO_test6_max_08_huri", "netbox_string2"] # ["netbox_string2", "DOMINO3_string2", "DOMINO3_1_string2", "DOMINO3_421_string2"] # ["DOMINO", "netbox", "DOMINO3_dip", "DOMINO3_1_dip", "DOMINO3_421_dip"] # ["DOMINO", "DOMINO3_1_dip", "DOMINO3_421_dip", "DOMINO3_dip",  "DOMINO3_421_3_dip", "netbox", "DOMINO3_421_033_dip", "DOMINO2_421_04_1000_dip", "DOMINO3_421_04_5000_dip", "DOMINO2_421_03447_1000_dip", "DOMINO2_45105_03447_1000_dip", "DOMINO3_421_008_string", "DOMINO4_421_012_huri", "DOMINO4_421_022_huri", "DOMINO4_421_018_5000_huri", "DOMINO4_45205_01944_1000_huri", "DOMINO4_per_013_1000_huri", "DOMINO2_per_035_1000_dip", "DOMINO2_per_02_1000_dip","DOMINO2_per_01_1000_dip", "DOMINO4_per_01_1000_huri", "DOMINO3_per_01_1000_string2", "DOMINO3_421_string2", "DOMINO_test_per_02_5000_string2","DOMINO_test_per_015_string3", "netbox_string2"] # ["netbox2_string", "DOMINO3_1_string"] # ["netbox_string2", "DOMINO3_string2","DOMINO3_1_string2", "DOMINO3_421_string2"] # ["DOMINO3_string2","DOMINO3_1_string2", "DOMINO3_421_string2"] # ["DOMINO3_1_dip", "DOMINO3_421_dip", "DOMINO3_dip"] # "DOMINO","DOMINO1","DOMINO2","DOMINO3","DOMINO5","DOMINO_2_3_5","DOMINO_DIP", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2", "netbox2_string", "netbox_string2","DOMINO2_string2", "DOMINO3_4", "DOMINO3_1_string","DOMINO3_1_string2", "DOMINO3_421_dip", "DOMINO3_1_dip", "DOMINO3_dip","DOMINO3_string2", "DOMINO3_421_string2"]