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
from utils.daemon_multiprocessing import func_star
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
from scipy.spatial.distance import jaccard
import seaborn as sns

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
    output=output.loc[output["GO id"].isin(output_md_ids.values)]
    output_md=output_md[output_md["GO id"].isin(output_md_ids)]

    real_pvals=output_md[output_md["hg_pval_max"].apply(lambda a: isFloat(a))].dropna()
    real_pvals["hg_pval_max"]=real_pvals["hg_pval_max"].values.astype(np.float)
    real_sig_ids = real_pvals[fdrcorrection0([10**-a for a in real_pvals["hg_pval_max"]], alpha=0.05, method='indep', is_sorted=False)[0]]["GO id"]

    n_of_bg=100
    n_intersections=[]
    n_sig_permuteds=[]
    n_sig_reals=[]
    for cur_bg in np.arange(n_of_bg):

        output["cur_bg_hg_pval"]=output.apply(lambda a: float(a["dist_n_samples"][1:-1].split(", ")[cur_bg]) if not pd.isnull(a["GO name"]) and a["dist_n_samples"].startswith("[") else np.nan ,axis=1)
        bg_pvals=output[output["cur_bg_hg_pval"] != -1].dropna()
        bg_sig_ids=bg_pvals[fdrcorrection0([10**-a for a in  bg_pvals["cur_bg_hg_pval"]], alpha=0.05, method='indep', is_sorted=False)[0]]["GO id"]

        n_intersections.append(len(set(real_sig_ids).intersection(set(bg_sig_ids))))
        n_sig_permuteds.append(len(bg_sig_ids))
        n_sig_reals.append(len(real_sig_ids))


    print "done combination: {}, {}".format(dataset,algo)
    return results.append((dataset,algo,n_sig_reals,n_sig_permuteds,n_intersections))


def calc_score(d,a,rs,ps,ins, df_summary):
    jaccards=[round(i/float(r+p+i+10e-7),2) for r,p,i in zip(rs,ps,ins)]
    if np.sum(rs)==0 and np.sum(ps)==0:
        df_summary.loc[constants.ALGOS_ACRONYM[a],d]=np.nan
    else:
        df_summary.loc[constants.ALGOS_ACRONYM[a],d]=1-np.mean(jaccards)

    # print constants.ALGOS_ACRONYM[a],constants.DATASETS_ACRONYM[d],rs,ps,ins,np.mean(jaccards)


def plot_jaccards(df_summary, algos, ax):

    df_summary=df_summary.loc[algos,:]

    ### heatmap ###
    # df_heatmap=pd.DataFrame(df_summary)
    # for a, jacs in df_summary.iterrows():
    #     df_heatmap.loc[a]=jacs.sort_values(ascending=False).values
    #
    # df_heatmap.columns=["" for a in df_heatmap.columns]
    # g=sns.heatmap(df_heatmap, cmap="Reds", ax=ax, mask=df_heatmap.isnull())
    # g.set_facecolor("#EEEEEE")

    ### heat dots ###
    cmap = matplotlib.cm.get_cmap('Reds')
    jac_min = np.nanmin(df_summary)
    jac_max = np.nanmax(df_summary)
    df_summary[np.isnan(df_summary)] = -1
    for a, jacs in df_summary.iterrows():
        ax.scatter(algos.index(a), 20,
                                 s=13200, c=(0, 0, 0, 1))
        for b, cur_jac in enumerate(jacs.sort_values(ascending=True)):
            ax.scatter(algos.index(a), 20,
                                 s=10000 * (1.1 - float(b) / len(jacs)) ** 2.5, c=(
                    cmap((cur_jac - jac_min) / (jac_max - jac_min)) if cur_jac != -1 else (0.8, 0.8, 0.8, 1)))
    ax.set_xticklabels([""] + algos, rotation=20)
    ax.set_yticklabels(["" for a in ax.get_yticks()], rotation=20)


if __name__ == "__main__":

    fig, ax = plt.subplots(2,1,figsize=(20,10))
    fig.subplots_adjust(left=0.2, right=0.8, wspace=2)

    for datasets, prefix, ax_index in [(["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"], "GE", 0) , (["brca", "crh", "scz", "tri", "t2d", "cad", "bmd", "hgt", "amd", "af"], "PASCAL_SUM", 1)]:
        algos=["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox", "keypathwayminer_INES_GREEDY","hotnet2"]
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
                params.append([false_positive_example, [algo, dataset, results]])
                # false_positive_example(algo=algo, dataset=dataset,df_summary=df_summary) # "Breast_Cancer.G50"
                # df_summary.loc[constants.ALGOS_ACRONYM[algo],constants.DATASETS_ACRONYM[dataset]]=np.random.rand()
                i += 1
        p.map(func_star, params)
        p.close()
        for d,a,rs,ps,ins in results:
            calc_score(d,a,rs,ps,ins, df_summary)
        df_summary.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation", "figure_11_{}_matrix.tsv".format(prefix)), sep='\t')
        plot_jaccards(df_summary, [constants.ALGOS_ACRONYM[a] for a in  algos], ax[ax_index])

    plt.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots","figure_11_heat.png".format(prefix)))


    # df_summary.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_11_{}_matrix.tsv".format(prefix)), sep='\t', index_label="algo")


