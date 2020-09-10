import sys
sys.path.insert(0, '../')
import matplotlib
matplotlib.use('Agg')

import pandas as pd

from fastsemsim.SemSim import *

import matplotlib
matplotlib.use("Agg")
import seaborn as sns

import constants

import argparse

import matplotlib.pyplot as plt


def main(algos=None, datasets=None, prefix="", cutoffs=[1.0, 2.0, 3.0, 4.0, 5.0], ax=None, axs_violin=None, title=""):

    ax.set_facecolor('#ffffff')
    ax.grid(color='gray')

    df_summary_agg = pd.DataFrame()
    for cutoff in cutoffs:
        df_summary=pd.read_csv(
            os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation", "richness_matrix_{}_{}.tsv".format(prefix, cutoff)),
            sep='\t', index_col=0).loc[algos,datasets]

        for a in set(df_summary.index).intersection(constants.ALGOS):
            for b in df_summary.columns:
                df_summary_agg=df_summary_agg.append({"algo": a, "dataset" : b, "cutoff": cutoff, "value" : df_summary.loc[a,b]}, ignore_index=True)


    if axs_violin is not None:
            for i_cutoff, cutoff in enumerate(cutoffs):
                my_order = df_summary_agg[df_summary_agg["cutoff"]==cutoff].groupby(by=["algo"])["value"].median().sort_values().index
                g=sns.violinplot(x="algo", y="value", data=df_summary_agg[df_summary_agg["cutoff"]==cutoff], ax=axs_violin[i_cutoff], order=my_order, palette={a: constants.COLORDICT[a] for a in my_order})
                g.set_xticklabels(g.get_xticklabels(), rotation=45)
    results = {}
    for cutoff in cutoffs:
        df_summary=pd.read_csv(
            os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation", "richness_matrix_{}_{}.tsv".format(prefix, cutoff)),
            sep='\t', index_col=0).loc[algos,datasets]

        for k, v in df_summary.iterrows():
            if k not in results:
                results[k] = []
            results[k].append(v.values)


    i=0
    df_mean=pd.DataFrame()
    df_median=pd.DataFrame()
    df_std = pd.DataFrame()
    for k,v in sorted(list(results.iteritems()),key=lambda a: a[0]):
        if len(v)>0:
            print k, v
        y_median=[ np.nanmedian(cur_measurement) for cur_measurement in v]
        y_mean = [np.nanmean(cur_measurement) for cur_measurement in v]
        y_stds= [np.std(cur_measurement) for cur_measurement in v]
        df_mean = df_mean.append(pd.Series(y_mean,index=[int(a) for a in cutoffs], name=k))
        df_median = df_median.append(pd.Series(y_median,index=[int(a) for a in cutoffs], name=k))
        df_std = df_std.append(pd.Series(y_stds,index=[int(a) for a in cutoffs], name=k))

        ax.plot(cutoffs,y_median,label=constants.ALGOS_ACRONYM[k], c=constants.COLORDICT[k])

        i += 1

    df_std.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"biological_richness_std_{}.tsv".format(prefix)), sep='\t')
    df_mean.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"biological_richness_average_{}.tsv".format(prefix)), sep='\t')
    df_median.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"biological_richness_median_{}.tsv".format(prefix)), sep='\t')

    ax.set_xlabel("similarity cutoff", fontdict={"size": 22})
    ax.set_ylabel("median # of non-redundant EV terms", fontdict={"size": 22})
    ax.set_title(title, fontdict={"size":22})
    patches = [Line2D([0], [0], marker='o', color='gray', label=constants.ALGOS_ACRONYM[a], markersize=12,
                      markerfacecolor=constants.COLORDICT[a], alpha=0.7) for a in algos]
    ax.legend(fontsize=17, loc=(0,1.1), ncol=2, facecolor='#ffffff')



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="TNFa_2,HC12,ROR_1,ERS_1,IEM,SHERA,SHEZH_1") # ",Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50 TNFa_2,HC12,SHERA,ROR_1,SHEZH_1,ERS_1,IEM"
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="dcem,jactivemodules_greedy,jactivemodules_sa,bionet,netbox,my_netbox_td,keypathwayminer_INES_GREEDY,hotnet2") # jactivemodules_greedy,jactivemodules_sa,bionet,netbox,my_netbox_td
    parser.add_argument('--pf', dest='pf', default=3)
    parser.add_argument('--base_folder', dest='base_folder', default=os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr","MAX"))
    parser.add_argument('--file_format', dest='file_format', default="emp_diff_modules_{}_{}_passed_oob.tsv")
    parser.add_argument('--sim_method', dest='sim_method', default="Resnik")
    args = parser.parse_args()

    prefix = args.prefix
    file_format = args.file_format
    base_folder= args.base_folder
    sim_method = args.sim_method
    datasets=["{}".format(x) for x in args.datasets.split(",")]
    algos = args.algos.split(",")
    cutoffs=[1.0,2.0,3.0,4.0]
    pf=int(args.pf)


    fig, axs = plt.subplots(1,2,figsize=(20,10))
    fig_violin, axs_violin = plt.subplots(2,len(cutoffs),figsize=(4*len(cutoffs)*2,10))

    main(algos=algos, datasets=datasets, prefix=prefix, cutoffs=cutoffs, ax=axs[0], axs_violin=axs_violin[0], title="GE")
