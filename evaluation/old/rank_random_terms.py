import sys
sys.path.insert(0, '../')
import matplotlib
matplotlib.use('Agg')

import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from infra import *
import pandas as pd
import multiprocessing
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
from scipy.stats import rankdata

def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


def false_positive_example(algo="jactivemodules_greedy",dataset="TNFa_2", results=None):
    print "running combination: {}, {}".format(dataset,algo)
    output=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"bg","emp_diff_modules_{}_{}.tsv".format(dataset, algo)), sep='\t', error_bad_lines=False)
    output_md = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"md", "emp_diff_modules_{}_{}_md.tsv".format(dataset, algo)), sep='\t', error_bad_lines=False)

    if not "GO id" in output_md.columns:
        df_summary.loc[dataset, algo]="0/0=0"
        return 

    output_md_ids=output_md.loc[np.logical_and.reduce([output_md["n_genes"].values > 5, output_md["n_genes"].values < 500])]['GO id']
    output_md_names=output_md.loc[np.logical_and.reduce([output_md["n_genes"].values > 5, output_md["n_genes"].values < 500])]['GO name']
    go_rank=pd.DataFrame(index=output_md_ids)
    go_rank.loc[:,"GO name"]=output_md_names.values
    go_rank.loc[:,"n_sig"]=0
    output=output.loc[output["GO id"].isin(output_md_ids.values)]
    output_md=output_md[output_md["GO id"].isin(output_md_ids)]

    real_pvals=output_md[output_md["hg_pval_max"].apply(lambda a: isFloat(a))].dropna()
    real_pvals["hg_pval_max"]=real_pvals["hg_pval_max"].values.astype(np.float)
    fdr=fdrcorrection0([10**-a for a in real_pvals["hg_pval_max"]], alpha=0.05, method='indep', is_sorted=False)
    real_sig_ids = real_pvals[fdr[0]]["GO id"]
    go_rank.loc[real_sig_ids,"real_sig"]=fdr[1]
    go_rank.loc[:,"real_sig_rank"]=rankdata(-go_rank.loc[:,"real_sig"])
    n_of_bg=5000
    n_intersections=[]
    n_sig_permuteds=[]
    n_sig_reals=[]
    params=[]
    p=multiprocessing.Pool(10)
    for cur_bg in np.arange(n_of_bg):
        params.append([output,cur_bg])

    bg_sig_ids_arr = p.map(get_bg_sig_ids, params)
    for bg_sig_ids in bg_sig_ids_arr:
        for a in bg_sig_ids:
            go_rank.loc[a,"n_sig"]=go_rank.loc[a,"n_sig"]+1

    go_rank.loc[:,"n_sig_rank"]=rankdata(-go_rank.loc[:,"n_sig"])


    return go_rank


def get_bg_sig_ids(args):
    output, cur_bg = args
    output["cur_bg_hg_pval"] = output.apply(
        lambda a: float(a["dist_n_samples"][1:-1].split(", ")[cur_bg]) if not pd.isnull(a["GO name"]) and a[
            "dist_n_samples"].startswith("[") else np.nan, axis=1)
    bg_pvals = output[output["cur_bg_hg_pval"] != -1].dropna()
    bg_sig_ids = bg_pvals[
        fdrcorrection0([10 ** -a for a in bg_pvals["cur_bg_hg_pval"]], alpha=0.05, method='indep', is_sorted=False)[0]]["GO id"]
    return bg_sig_ids


def calc_score(d,a,rs,ps,ins, df_summary):
    jaccards=[round(i/float(r+p+i+10e-7),2) for r,p,i in zip(rs,ps,ins)]
    if np.sum(rs)==0 and np.sum(ps)==0:
        df_summary.loc[constants.ALGOS_ACRONYM[a],constants.DATASETS_ACRONYM[d]]=np.nan
    else:
        df_summary.loc[constants.ALGOS_ACRONYM[a],constants.DATASETS_ACRONYM[d]]=1-np.mean(jaccards)

    # print constants.ALGOS_ACRONYM[a],constants.DATASETS_ACRONYM[d],rs,ps,ins,np.mean(jaccards)


if __name__ == "__main__":

    fig, ax = plt.subplots(2,1,figsize=(20,10))
    fig.subplots_adjust(left=0.2, right=0.8, wspace=2)
    # for datasets, prefix, ax_index in [(["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"], "GE", 0) , (["brca", "crh", "scz", "tri", "t2d", "cad", "cmd", "hgt", "amd", "af"], "PASCAL_SUM", 1)]:
    #     algos=["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox", "keypathwayminer_INES_GREEDY","hotnet2"]
    for datasets, prefix, ax_index in [(["scz"], "GE", 0) ]:
        algos=["jactivemodules_greedy"] # ["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox", "keypathwayminer_INES_GREEDY","hotnet2"]
        ax[ax_index].set_facecolor('#ffffff')
        ax[ax_index].set_title(prefix if prefix=="GE" else "GWAS", fontdict={"size":30})
        params=[]
        p=multiprocessing.Pool(80)
        manager = multiprocessing.Manager()
        results = manager.list()
        i=0
        df_summary=pd.DataFrame()
        for dataset in datasets:
            for algo in algos:
                go_rank=false_positive_example(algo, dataset, results)
                go_rank.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"go_ranks_{}_{}.tsv".format(dataset,algo)), sep='\t')

