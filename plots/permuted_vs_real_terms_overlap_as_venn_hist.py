from matplotlib import style

import sys
sys.path.insert(0, '../')
import matplotlib
matplotlib.use('Agg')

import seaborn as sns
sns.set(color_codes=True)
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from infra import *
import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt

from statsmodels.sandbox.stats.multicomp import fdrcorrection0

from matplotlib_venn import venn2

fontsize=35

font = {'size'   : fontsize}
mpl.rc('xtick', labelsize=fontsize)    # fontsize of the tick labels
mpl.rc('ytick', labelsize=fontsize)


def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


def false_positive_example(algo="jactivemodules_greedy",dataset="TNFa_2", is_corrected=False, axs=None):

    output=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"bg","emp_diff_modules_{}_{}.tsv".format(dataset, algo)), sep='\t', error_bad_lines=False)
    output.index=output.loc[:,"GO id"]
    output_md = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "md", "emp_diff_modules_{}_{}_md.tsv".format(dataset, algo)), sep='\t', error_bad_lines=False)
    try:
        output_md.index=output_md["GO id"]
        output_md_ids=output_md.loc[np.logical_and.reduce([output_md["n_genes"].values > 5, output_md["n_genes"].values < 500])]['GO id']
        output_md=output_md[output_md["GO id"].isin(output_md_ids)]
        output=output.loc[output["GO id"].isin(output_md_ids.values)]

    except:
        pass




    axs[0].set_facecolor('#ffffff')
    axs[1].set_facecolor('#ffffff')

    x_label = "-log10(pval)"
    if is_corrected:
        x_label = "-log10(qval)"

    n_of_bg=1
    for cur_bg in np.arange(n_of_bg):

        output["cur_bg_hg_pval"]=output.apply(lambda a: float(a["dist_n_samples"][1:-1].split(", ")[cur_bg+11] if a["dist_n_samples"].count(',') > 0 else 0) if a["dist_n_samples"].startswith("[") else np.nan ,axis=1)

        bg_pvals=output[output["cur_bg_hg_pval"] != -1] # .dropna()

        values=bg_pvals["cur_bg_hg_pval"].values
        if is_corrected:
            values = -np.log10(fdrcorrection0(10**-values, alpha=0.05, method='indep', is_sorted=False)[1])

        sns.distplot([a for a in values if a > -np.log10(0.05)],norm_hist=False, kde=False, label="Permuted dataset", ax=axs[0], bins=50)

    bg_sig_ids=bg_pvals[fdrcorrection0([10**-a for a in  bg_pvals["cur_bg_hg_pval"]], alpha=0.05, method='indep', is_sorted=False)[0]]["GO id"]

    # plt.clf()
    real_pvals=output_md[output_md["hg_pval_max"].apply(lambda a: isFloat(a))].dropna()
    real_pvals["hg_pval_max"]=real_pvals["hg_pval_max"].values.astype(np.float)

    values=real_pvals["hg_pval_max"][real_pvals["hg_pval_max"]!=-1]
    if is_corrected:
        values=-np.log10(fdrcorrection0(10**-values, alpha=0.05, method='indep', is_sorted=False)[1])

    sns.distplot(values, norm_hist=False, kde=False, label="Original dataset", ax=axs[0], bins=50)
    axs[0].set_yscale('log')

    axs[0].set_ylabel("# GO terms", fontsize=fontsize)
    axs[0].set_xlabel(x_label, fontsize=fontsize)
    axs[0].set_title(constants.ALGOS_ACRONYM[algo], fontsize=fontsize)
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "pval_dist_real_{}_{}.png".format(dataset, constants.ALGOS_ACRONYM[algo])))

    try:
        real_sig_ids = real_pvals[fdrcorrection0([10**-a for a in real_pvals["hg_pval_max"]], alpha=0.05, method='indep', is_sorted=False)[0]]["GO id"]
    except:
        real_sig_ids =[]
    out=venn2([set(real_sig_ids), set(bg_sig_ids)], set_labels=('Original\ndataset', 'Permuted\ndataset'), ax=axs[1])
    for text in out.set_labels:
        text.set_fontsize(fontsize)
    for text in out.subset_labels:
        if text is not None:
            text.set_fontsize(fontsize)
    # axs[1].set_title("overlap of enriched GO terms\nbetween original and permuted datasets", fontsize=fontsize)
    axs[1].set_title("Enriched GO terms", fontsize=fontsize)



        # cur_bg_pvals = []
        # agg_counter=0
        # for k,v in bg_dict:
        #     cur_bg_pvals.append(v[cur_bg])
        #     print "agg_counter: {}/{}".format(agg_counter, len(bg_dict))
        #     agg_counter+=1





if __name__ == "__main__":

    for is_corrected in [True]:
        suffix = ("qval" if is_corrected else "pval")
        datasets=["cbx"]# ["TNFa_2", "Schizophrenia.G50"]  #, "Schizophrenia.G50"
        algos= ["jactivemodules_greedy","jactivemodules_sa","netbox","bionet","hotnet2", "keypathwayminer_INES_GREEDY"] #
        i=0
        fig,axs=plt.subplots(3, 4, figsize=(56, 12*len(algos)/2))
        for dataset in datasets:
            for i_a, algo in enumerate(algos):
                print "current solution: {}-{}".format(dataset,algo)
                false_positive_example(algo=algo, dataset=dataset,is_corrected=is_corrected, axs=axs[i_a%3,np.arange(2) if i_a<3 else np.arange(2,4) ] ) # "Breast_Cancer.G50"
                i += 1

        fig.legend(["Permuted dataset", "Original dataset"], prop={'size': fontsize}, facecolor="#FFFFFF", loc=(0,-0.5))
        # fig.text(0.01,0.98, "A:", weight='bold',fontsize=fontsize)
        # fig.text(0.55, 0.98, "B:", weight='bold',fontsize=fontsize)
        # fig.text(0.01,0.5, "C:", weight='bold',fontsize=fontsize)
        # fig.text(0.55, 0.5, "D:", weight='bold',fontsize=fontsize)
        # fig.text(0.01, 0.5, "E:", weight='bold', fontsize=fontsize)
        # fig.text(0.55, 0.5, "F:", weight='bold', fontsize=fontsize)
        # fig.text(0.01, 0.25, "G:", weight='bold', fontsize=fontsize)
        # fig.text(0.55, 0.25, "H:", weight='bold', fontsize=fontsize)
        # fig.tight_layout()
        fig.subplots_adjust(wspace=0.34, hspace=0.22)
        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "figure_11_venn_hist_{}.png".format(suffix)))
