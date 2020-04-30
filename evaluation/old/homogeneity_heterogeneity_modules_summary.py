import matplotlib
matplotlib.use('Agg')

import sys
sys.path.insert(0, '../')

import argparse
import seaborn as sns
sns.set(color_codes=True)
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from infra import *
import pandas as pd
import numpy as np
from statsmodels.sandbox.stats.multicomp import fdrcorrection0


def calc_emp_pval(cur_rv, cur_dist):
    pos = np.size(cur_dist) - np.searchsorted(np.sort(cur_dist), cur_rv, side='left')
    return pos / float(np.size(cur_dist))


def calc_modules_ehr(algo_sample = None, dataset_sample = None, prefix=None, terms_file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "MAX", "emp_diff_{}_{}_oob.tsv"),
         modules_file_name=os.path.join("/media/hag007/Data/bnet/output/{}_{}/{}/modules_summary.tsv"), emp_ratio_th=0.5):
    statistics = {}
    full_data = pd.DataFrame()

    output_terms = pd.read_csv(
        terms_file_name.format(dataset_sample, algo_sample),
        sep='\t', index_col=0).dropna()

    output_modules = pd.read_csv(
        modules_file_name.format(prefix,dataset_sample, algo_sample),
        sep='\t', index_col=0).dropna()



    output_terms = output_terms.rename(columns={"filtered_pval": "hg_pval_max"})
    filtered_genes=output_terms.loc[np.logical_and.reduce([output_terms["n_genes"].values > 5, output_terms["n_genes"].values < 500]), ["GO name","hg_pval_max", "emp_pval_max", "passed_oob_permutation_test"]]


    print "total n_genes with pval:{}/{}".format(np.size(filtered_genes["hg_pval_max"].values), 7435)

    sorted_genes_hg = filtered_genes.sort_values(by=['hg_pval_max'], ascending=False)
    sig_genes_hg_pval=np.append(sorted_genes_hg["hg_pval_max"].values,np.zeros(7435-np.size(sorted_genes_hg["hg_pval_max"].values)))
    sig_genes_hg_pval = [10**(-x) for x in sig_genes_hg_pval]
    fdr_results = fdrcorrection0(sig_genes_hg_pval, alpha=0.05, method='indep', is_sorted=False)
    n_hg_true = len([cur for cur in fdr_results[0] if cur])
    sig_hg_genes= sorted_genes_hg.iloc[:n_hg_true, :] if n_hg_true > 0 else 0
    if type(sig_hg_genes)==int:
        HG_CUTOFF = 0
        statistics["hg_cutoff"] = 0
        statistics["hg_corrected"] = 0
        print "HG cutoff: {}, n={}".format(HG_CUTOFF, 0)
    else:
        HG_CUTOFF = sig_hg_genes.iloc[- 1]["hg_pval_max"]
        statistics["hg_cutoff"]=round(HG_CUTOFF,2)
        statistics["hg_corrected"]=len(sig_hg_genes.index)
        print "HG cutoff: {}, n={}".format(HG_CUTOFF, len(sig_hg_genes.index))

    sorted_genes_emp = filtered_genes.sort_values(by=['emp_pval_max'])
    sorted_genes_emp.loc[sorted_genes_emp['emp_pval_max']==0,'emp_pval_max']=1.0/5000
    sig_genes_emp_pval = sorted_genes_emp["emp_pval_max"].values
    fdr_results = fdrcorrection0(sig_genes_emp_pval, alpha=0.05, method='indep', is_sorted=False)
    n_emp_true = len([cur for cur in fdr_results[0] if cur])
    sig_emp_genes = sorted_genes_emp.iloc[:n_emp_true, :]
    EMP_CUTOFF = sig_emp_genes.iloc[- 1]["emp_pval_max"] if n_emp_true > 0 else 0
    statistics["emp_cutoff"]=EMP_CUTOFF
    statistics["emp_corrected"]=len(sig_emp_genes.index)
    print "EMP cutoff: {}, n={}".format(EMP_CUTOFF, len(sig_emp_genes.index))

    tps={}
    fps = {}
    if output_terms.shape[0]==0:
        n_modules=0
    else:
        n_modules=output_terms.iloc[0]["hg_pval"].count(',')+1
    statistics["n_modules"]=n_modules
    for a in range(n_modules):
        tps[a]=[]
        fps[a]=[]

    for go_id , cur in output_terms.iterrows():
        hgs=list(np.array(cur["hg_pval"][1:-1].split(", "),dtype=np.float32))
        emps=list(np.array(cur["emp_pval"][1:-1].split(", "),dtype=np.float32))
        for i, v in enumerate(zip(hgs,emps)):
            hg,emp=v
            if hg >= HG_CUTOFF and emp <= EMP_CUTOFF:
                tps[i].append("{}: {}".format(go_id, cur["GO name"]))
            elif hg >= HG_CUTOFF and emp > EMP_CUTOFF:
                fps[i].append("{}: {}".format(go_id, cur["GO name"]))

    real_modules_counter=0
    for a in range(n_modules):
        n_hg_terms=len(tps[a]) + len(fps[a])
        # print "module_{}: {} ({}/{})".format(a, float(len(tps[a]))/n_hg_terms, len(tps[a]), n_hg_terms)
        statistics["module_{}_emp_ratio".format(a)]=round(float(len(tps[a]))/max(len(tps[a])+len(fps[a]),1),2)
        statistics["module_{}_tp".format(a)] = len(tps[a])
        statistics["module_{}_fp".format(a)] = len(fps[a])
        statistics["module_{}_total".format(a)] = len(tps[a])+len(fps[a])
        print dataset_sample, algo_sample
        statistics["module_{}_size".format(a)] = output_modules.loc[a, '#_genes']

        if statistics["module_{}_emp_ratio".format(a)] >emp_ratio_th:
            real_modules_counter+=1
            # for cur in tps[a]:
            #     print cur
        full_data.loc["{}_{}_module_{}".format(dataset_sample, algo_sample, a), "EHR"] = round(float(len(tps[a])) / max(len(tps[a]) + len(fps[a]), 1),2)
        full_data.loc["{}_{}_module_{}".format(dataset_sample, algo_sample, a),"tp"]="\n".join(tps[a])
        full_data.loc["{}_{}_module_{}".format(dataset_sample, algo_sample, a),"fp"]="\n".join(fps[a])
        full_data.sort_index(inplace=True)
    statistics["real_modules_ratio"]=round(float(real_modules_counter)/max(n_modules,1),2)
    # print "real_modules_counter: {}. ratio: {}".format(real_modules_counter, statistics["real_modules_ratio"])


    for cur_key in np.sort(statistics.keys()):
        print "{}: {}".format(cur_key, statistics[cur_key])

    file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR, "full_modules_report_{}_{}.tsv".format(algo_sample, dataset_sample))
    full_data.to_csv(file_name,sep='\t')

    with file(file_name,'a+') as f:
        f.write("\nsolution statistics: \n")
        for k in np.sort(statistics.keys()):
            if not k.startswith("modules_"):
                f.write("{}\t{}\n".format(k,statistics[k]))

        for k in np.sort(statistics.keys()):
            if k.startswith("modules_"):
                f.write("{}\t{}\n".format(k,statistics[k]))


    return tps, fps, sig_hg_genes, sig_emp_genes, statistics, full_data



def plot_modules_ehr_summary(prefix, datasets, algos, base_folder, terms_file_name_format):


    terms_limit = 0
    df_rank_matrix = pd.DataFrame()
    df_matrix = pd.DataFrame()
    df_count_matrix = pd.DataFrame()
    df_all = pd.DataFrame()
    df_statistics = pd.DataFrame()
    df_full_data = pd.DataFrame()
    for cur_ds in datasets:
        df_ds=pd.DataFrame()
        for cur_alg in algos:
            terms_file_name = os.path.join(base_folder, terms_file_name_format.format(cur_ds, cur_alg))
            modules_file_name = os.path.join(constants.OUTPUT_GLOBAL_DIR, "{}_{}".format(prefix, cur_ds), cur_alg, "modules_summary.tsv")

            tps, fps, sig_hg_genes, sig_emp_genes, statistics, full_data =calc_modules_ehr(cur_alg, cur_ds, prefix,
                                                                   terms_file_name=terms_file_name,
                                                                   modules_file_name=modules_file_name)
            statistics["algo"]=cur_alg
            statistics["dataset"]=cur_ds
            statistics["id"]="{}_{}".format(cur_alg, cur_ds)

            full_data["algo"] = cur_alg
            full_data["dataset"] = cur_ds

            df_statistics=df_statistics.append(statistics, ignore_index=True)
            df_full_data=df_full_data.append(full_data)


    df_statistics=df_statistics.set_index("id")
    df_statistics.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"modules_statistics_{}.tsv".format(prefix)), sep='\t')
    df_full_data.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "modules_full_data_{}.tsv".format(prefix) ), sep='\t')

    df_statistics = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "modules_statistics_{}.tsv".format(prefix)), sep='\t',
                                index_col=0)
    df_full_data = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "modules_full_data_{}.tsv".format(prefix)), sep='\t',
                                index_col=0)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="Breast_Cancer.G50,Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50,Coronary_Artery_Disease.G50,Bone_Mineral_Density.G50,Height1.G50,Alzheimer.G50,Age_Related_Macular_Degeneration.G50,Atrial_Fibrillation.G50")  # TNFa_2,HC12,SHERA,SHEZH_1,ROR_1,ERS_1,IEM,APO,CBX,IFT ## Breast_Cancer.G50,Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50,Coronary_Artery_Disease.G50,Bone_Mineral_Density.G50,Height1.G50,Alzheimer.G50,Age_Related_Macular_Degeneration.G50,Atrial_Fibrillation.G50
    parser.add_argument('--prefix', dest='prefix', default="PASCAL_SUM") # GE   PASCAL_SUM
    parser.add_argument('--base_folder', dest='base_folder', default=os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr/MAX"))
    parser.add_argument('--terms_file_name_format', dest='terms_file_name_format', default="emp_diff_modules_{}_{}_passed_oob.tsv")
    parser.add_argument('--algos', dest='algos',
                        default="dcem,jactivemodules_greedy,jactivemodules_sa,bionet,netbox,keypathwayminer_INES_GREEDY") # hotnet2 ")  # ,keypathwayminer_INES_GREEDY,hotnet2,my_netbox_td,hotnet2,keypathwayminer_INES_GREEDY

    args = parser.parse_args()

    datasets = args.datasets.split(",")
    algos = args.algos.split(",")
    prefix = args.prefix
    base_folder = args.base_folder
    terms_file_name_format = args.terms_file_name_format

    plot_modules_ehr_summary(prefix, datasets, algos, base_folder, terms_file_name_format)
