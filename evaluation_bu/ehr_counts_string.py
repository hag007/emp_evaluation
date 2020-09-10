import sys
sys.path.insert(0, '../')

import constants
import os
import argparse
import pandas as pd
import numpy as np
from statsmodels.sandbox.stats.multicomp import fdrcorrection0

def ehr_for_solution(algo_sample = None, dataset_sample = None,tsv_file_name=None, filtered_go_ids_file=None, hg_th=0.05):

    try:
        output_md = pd.read_csv(
            tsv_file_name.format(dataset_sample, algo_sample),
            sep='\t', index_col=0).dropna()
    except Exception, e:
        print e
        return 0, 0

    filtered_go_ids=open(filtered_go_ids_file,'r').read().split()
    filtered_genes=output_md.reindex(filtered_go_ids).loc[:,["GO name","hg_pval_max", "emp_pval_max", "passed_oob_permutation_test"]]

    print "total n_genes with pval:{}/{}".format(np.size(filtered_genes["hg_pval_max"].values), len(filtered_go_ids))

    sorted_genes_hg = filtered_genes.sort_values(by=['hg_pval_max'], ascending=False)
    sig_genes_hg_pval=np.append(sorted_genes_hg["hg_pval_max"].values,np.zeros(len(filtered_go_ids)-np.size(sorted_genes_hg["hg_pval_max"].values)))
    sig_genes_hg_pval = [10**(-x) for x in sig_genes_hg_pval]
    fdr_results = fdrcorrection0(sig_genes_hg_pval, alpha=hg_th, method='indep', is_sorted=False)
    n_hg_true = len([cur for cur in fdr_results[0] if cur])
    sig_hg_genes= sorted_genes_hg.iloc[:n_hg_true, :] 
    if len(sig_hg_genes) == 0:
        HG_CUTOFF = 0
        print "HG cutoff: {}, n=0".format(HG_CUTOFF)

    else:
        HG_CUTOFF = 10**(-sig_hg_genes.iloc[- 1]["hg_pval_max"])
        print "HG cutoff: {}, n={}".format(HG_CUTOFF, len(sig_hg_genes.index))
    print(len(sig_hg_genes.index))
    sorted_genes_emp = filtered_genes.sort_values(by=['emp_pval_max'])[:len(sig_hg_genes.index)] 
    sorted_genes_emp.loc[sorted_genes_emp['emp_pval_max']==0,'emp_pval_max']=1.0/5000
    sig_genes_emp_pval = sorted_genes_emp["emp_pval_max"].values
    fdr_results = fdrcorrection0(sig_genes_emp_pval, alpha=0.05, method='indep', is_sorted=False)
    n_emp_true = sum(fdr_results[0])
    sig_emp_genes = sorted_genes_emp.iloc[:n_emp_true, :]
    EMP_CUTOFF = sig_emp_genes.iloc[- 1]["emp_pval_max"] if n_emp_true > 0 else 0
    print "EMP cutoff: {}, n={}".format(EMP_CUTOFF, len(sig_emp_genes.index))

    return sig_hg_genes, sig_emp_genes


def main(datasets, algos, prefix):

    hg_th=0.05
    df_matrix = pd.DataFrame()
    df_count_matrix = pd.DataFrame()
    for cur_ds in datasets:
        df_ds=pd.DataFrame()
        for cur_alg in algos:
            sig_hg_genes, sig_emp_genes=ehr_for_solution(cur_alg, cur_ds,tsv_file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR,"oob", "emp_diff_modules_{}_{}_passed_oob.tsv"), filtered_go_ids_file=os.path.join(constants.GO_DIR,"filtered_go_terms.txt"), hg_th=hg_th) # _gwas

            df_ds.loc["{}_{}".format(cur_ds,cur_alg), "algo"]=cur_alg
            df_ds.loc["{}_{}".format(cur_ds,cur_alg), "dataset"]=cur_ds
            if type(sig_emp_genes)==int or type(sig_hg_genes)==int:
                df_ds.loc["{}_{}".format(cur_ds, cur_alg), "n_emp"] = 0
                df_ds.loc["{}_{}".format(cur_ds, cur_alg), "n_hg"] = 0
            else:
                df_ds.loc["{}_{}".format(cur_ds,cur_alg),"n_emp"]=np.sum(sig_emp_genes["passed_oob_permutation_test"].apply(lambda x: np.any(np.array(x[1:-1].split(', '),dtype=np.bool))).values)
                df_ds.loc["{}_{}".format(cur_ds,cur_alg), "n_hg"] = len(sig_hg_genes.index)
            df_ds.loc["{}_{}".format(cur_ds,cur_alg), "ratio"] = round(float(df_ds.loc["{}_{}".format(cur_ds,cur_alg), "n_emp"])/max(df_ds.loc["{}_{}".format(cur_ds,cur_alg), "n_hg"],1),3)
            df_count_matrix.loc[cur_alg, cur_ds] = df_ds.loc["{}_{}".format(cur_ds, cur_alg), "n_emp"]
            df_matrix.loc[cur_alg, cur_ds] = df_ds.loc["{}_{}".format(cur_ds, cur_alg), "ratio"]

    df_matrix.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation", "ehr_matrix_{}.tsv".format(prefix)), sep='\t')
    df_count_matrix.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation", "count_matrix_{}.tsv".format(prefix)), sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="TNFa_2,HC12,ROR_1,SHERA,SHEZH_1,ERS_1,IEM,APO,CBX,IFT") # "TNFa_2,HC12,ROR_1,SHERA,SHEZH_1,ERS_1,IEM,APO,CBX,IFT" Breast_Cancer.G50,Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50
    parser.add_argument('--prefix', dest='prefix', default="GE") # PASCAL_SUM   GE
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy,jactivemodules_sa,bionet,netbox,keypathwayminer_INES_GREEDY,dcem2") # ,dcem2,dcem3,dcem4,my_netbox_td,hotnet2

    args = parser.parse_args()
    datasets = args.datasets.split(",")
    algos = args.algos.split(",")
    prefix = args.prefix


    datasets=["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"]
    algos=['DOMINO4', 'netbox2_string'] # ["DOMINO", "DOMINO2", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    prefix="GE"
    main(datasets,algos,prefix)

    datasets=["brca", "crh", "scz", "tri", "t2d", "cad", "bmd", "hgt", "amd", "af"]
    algos=['DOMINO4', 'netbox2_string'] # ["DOMINO", "DOMINO2", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    prefix="PASCAL_SUM"
    main(datasets,algos,prefix)
