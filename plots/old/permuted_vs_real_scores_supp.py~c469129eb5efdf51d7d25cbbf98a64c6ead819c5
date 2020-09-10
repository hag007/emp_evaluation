import sys
sys.path.insert(0, '../')

import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from infra import *
import pandas as pd
import multiprocessing
from utils.daemon_multiprocessing import func_star
from statsmodels.sandbox.stats.multicomp import fdrcorrection0

def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


def false_positive_example(algo="jactivemodules_greedy",dataset="TNFa_2", results=None):
    print "running combination: {}, {}".format(dataset,algo)
    output=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr","MAX","emp_diff_modules_{}_{}.tsv".format(dataset, algo)), sep='\t', error_bad_lines=False)
    output_md = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "MAX", "emp_diff_modules_{}_{}_md.tsv".format(dataset, algo)), sep='\t', error_bad_lines=False)

    if not "GO id" in output_md.columns:
        df_summary.loc[dataset, algo]="0/0=0"
        return 

    output_md_ids=output_md.loc[np.logical_and.reduce([output_md["n_genes"].values > 5, output_md["n_genes"].values < 500])]['GO id']
    output=output.loc[output["GO id"].isin(output_md_ids.values)]
    output_md=output_md[output_md["GO id"].isin(output_md_ids)]

    n_of_bg=1
    for cur_bg in np.arange(n_of_bg):

        output["cur_bg_hg_pval"]=output.apply(lambda a: float(a["dist_n_samples"][1:-1].split(", ")[cur_bg]) if not pd.isnull(a["GO name"]) and a["dist_n_samples"].startswith("[") else np.nan ,axis=1)
        bg_pvals=output[output["cur_bg_hg_pval"] != -1].dropna()


    bg_sig_ids=bg_pvals[fdrcorrection0([10**-a for a in  bg_pvals["cur_bg_hg_pval"]], alpha=0.05, method='indep', is_sorted=False)[0]]["GO id"]

    real_pvals=output_md[output_md["hg_pval_max"].apply(lambda a: isFloat(a))].dropna()
    real_pvals["hg_pval_max"]=real_pvals["hg_pval_max"].values.astype(np.float)

    values=real_pvals["hg_pval_max"][real_pvals["hg_pval_max"]!=-1]
    # if is_corrected:
    #     values=-np.log10(fdrcorrection0(10**-values, alpha=0.05, method='indep', is_sorted=False)[1])


    real_sig_ids = real_pvals[fdrcorrection0([10**-a for a in real_pvals["hg_pval_max"]], alpha=0.05, method='indep', is_sorted=False)[0]]["GO id"]
    # out=venn2([set(real_sig_ids), set(bg_sig_ids)], set_labels=('Original\ndataset', 'Permuted\ndataset'), ax=axs[1])
    n_intersection = len(set(real_sig_ids).intersection(set(bg_sig_ids)))
    n_real = len(real_sig_ids)
    # _summary.loc[dataset,algo]=round(float(n_intersection)/max(n_real,1),2)
    print "done combination: {}, {}".format(dataset,algo)
    return results.append((dataset,algo,round(float(n_intersection)/max(n_real,1),2)))


if __name__ == "__main__":
       
    for datasets,prefix in [(["Breast_Cancer.G50"], "PASCAL_SUM"), (["ERS_1"], "GE")]:
        algos=["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox", "keypathwayminer_INES_GREEDY","hotnet2"] 
        # if prefix=="GE": continue
        params=[]                
        p=multiprocessing.Pool(30)
        manager = multiprocessing.Manager()
        results = manager.list()
        i=0
        df_summary=pd.DataFrame()
        for dataset in datasets:
            for algo in algos:
                params.append([false_positive_example, [algo, dataset, results]])
                # false_positive_example(algo=algo, dataset=dataset,df_summary=df_summary) # "Breast_Cancer.G50"
                i += 1
        p.map(func_star, params)
        for d,a,v in results:
             df_summary.loc[d,a]=v
        df_summary.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_11_{}_matrix.tsv".format(prefix)), sep='\t', index_label="algo")

