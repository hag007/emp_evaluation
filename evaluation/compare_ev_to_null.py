import sys
sys.path.insert(0, '../')

import constants
import os
import argparse
import pandas as pd
import numpy as np
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
from scipy.stats import hypergeom, fisher_exact

def ehr_for_solution(algo_sample = None, dataset_sample = None,tsv_file_name=None, filtered_go_ids_file=None, sig_in_null_file=None, hg_th=0.05):

    try:
        output_oob = pd.read_csv(
            tsv_file_name.format(dataset_sample, algo_sample),
            sep='\t', index_col=0).dropna()
    except Exception, e:
        print e
        return 0, 0

    filtered_go_ids=open(filtered_go_ids_file,'r').read().split()
    sig_in_null=open(sig_in_null_file,'r').read().split()
    filtered_genes=output_oob.reindex(filtered_go_ids).loc[:,["GO name", "is_emp_pval_max_significant"]].dropna()
    filtered_genes.loc[:, "sig_in_null"] = filtered_genes.index.isin(sig_in_null)
    true_ev_terms=filtered_genes.loc[filtered_genes.loc[:,"is_emp_pval_max_significant"].astype(np.bool),:]
    true_null_terms=filtered_genes.loc[filtered_genes.loc[:,"sig_in_null"].astype(np.bool),:]
    false_ev_terms=filtered_genes.loc[filtered_genes.loc[:,"is_emp_pval_max_significant"].astype(np.bool),:]
    ehr=round(np.sum(filtered_genes.loc[:,"is_emp_pval_max_significant"].astype(np.bool))/float(filtered_genes.loc[:,"is_emp_pval_max_significant"].shape[0]), 2)
    print("{}_{} (EHR={}):".format(algo_sample, dataset_sample, ehr))
    # TP
    tp=np.sum(np.logical_and(~filtered_genes.loc[:, "is_emp_pval_max_significant"].astype(np.bool), filtered_genes.loc[:, "sig_in_null"]))
    print('TP: {}'.format(tp))
    # FP
    fp=np.sum(np.logical_and(filtered_genes.loc[:, "is_emp_pval_max_significant"].astype(np.bool), filtered_genes.loc[:, "sig_in_null"]))
    print('FP: {}'.format(fp))
    # TN
    tn=np.sum(np.logical_and(filtered_genes.loc[:, "is_emp_pval_max_significant"].astype(np.bool), ~filtered_genes.loc[:, "sig_in_null"]))
    print('TN: {}'.format(tn))
    # FN
    fn=np.sum(np.logical_and(~filtered_genes.loc[:, "is_emp_pval_max_significant"].astype(np.bool), ~filtered_genes.loc[:, "sig_in_null"]))
    print('FN: {}'.format(fn))

    # OR
    tpr=round(float(tp)/max(tp+fn, 10e-5), 2)
    fpr=round(float(fp)/max(fp+tn, 10e-5), 2)
    OR=round(tpr/max(fpr, 10e-5), 2)

    # HG
    sig_score = fisher_exact([[tp, fp], [fn, tn]], alternative='greater')[1] #hypergeom.sf(tp, tp+fp+tn+fn, tp+fn, tp+fp) + hypergeom.pmf(tp, tp+fp+tn+fn, tp+fn, tp+fp)

    # precision

    precision=round(float(tp)/max(1,tp+fp), 2)
    print('precision: {}'.format(precision))
    # recall
    recall=round(float(tp)/max(1,tp+fn), 2)
    print('recall: {}'.format(recall))
    # f1
    f1=round(2*(precision*recall)/max((precision+recall),1), 2)
    print('f1: {}'.format(f1))
    return {"id": "{}_{}".format(algo_sample, dataset_sample), "dataset": dataset_sample, "algo": algo_sample, "EHR": ehr, "precision": precision, "recall": recall, "f1" : f1, "TP": tp, "FP" : fp, "TN": tn , "FN" :fn, "TPR" : tpr, "FPR" : fpr, "OR" : OR, "pval": sig_score}

def main(datasets, algos, prefix, hg_th=0.05):
    df=pd.DataFrame()
    for cur_ds in datasets:
        for cur_alg in algos:
            df=df.append(ehr_for_solution(cur_alg, cur_ds,tsv_file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR,"oob", "emp_diff_modules_{}_{}_passed_oob.tsv"), filtered_go_ids_file=os.path.join(constants.GO_DIR,"filtered_go_terms.txt"), sig_in_null_file=os.path.join(constants.GO_DIR,"sig_in_null.txt"),  hg_th=hg_th), ignore_index=True)
    df.index=df.loc[:,"id"]
    df=df.drop(["id"], axis=1)
    df.loc[:,"pval"][df.loc[:,"pval"].isnull()]=1
    fdr=fdrcorrection0(df.loc[:,"pval"])[1]
    df.loc[:,"qval"]=fdr
    df.loc[:,["TP", "FP", "FN", "TN", "EHR", "f1", "precision", "recall", "TPR", "FPR", "OR", "pval", "qval"]].to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "ev_vs_null_{}.tsv".format(prefix)),sep='\t')

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