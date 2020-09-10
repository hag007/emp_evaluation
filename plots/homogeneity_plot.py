import sys
sys.path.insert(0, '../')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from fastsemsim.SemSim import *
import constants

from matplotlib.lines import Line2D

import seaborn as sns


def variability_plot(algos, datasets, base_folder, homogeneity_file_format, suffix, cutoffs, ax, axs_violin, title=""):

    ax.set_facecolor('#ffffff')
    ax.grid(color='gray')

    homogeneities = pd.DataFrame()

    for i_cutoff, cutoff in enumerate(cutoffs):
        df_homogeneity = pd.read_csv(os.path.join(base_folder, homogeneity_file_format.format(suffix, cutoff)),
                                     sep='\t', index_col=0)
        df_homogeneity=df_homogeneity.loc[algos,datasets]
        df_homogeneity[df_homogeneity.isna()] = 0

        homogeneities = pd.concat([homogeneities, df_homogeneity.mean(axis=1)], axis=1)

        if axs_violin is not None:
            df_summary_agg=pd.DataFrame()
            for algo in set(df_homogeneity.index).intersection(constants.ALGOS):
                for dataset in df_homogeneity.columns:
                    df_summary_agg=df_summary_agg.append({"algo":algo, "dataset": dataset, "value": df_homogeneity.loc[algo,dataset]},ignore_index=True)

            df_summary_agg=df_summary_agg.dropna(axis=0)
            my_order = df_summary_agg.groupby(by=["algo"])["value"].median().sort_values().index

            g = sns.violinplot(x="algo", y="value", data=df_summary_agg,
                               ax=axs_violin[i_cutoff], order=my_order,
                               palette={a: constants.COLORDICT[a] for a in my_order})
            g.set_xticklabels(g.get_xticklabels(), rotation=45)

    homogeneities = homogeneities.loc[algos]


    homogeneities.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"homogeneity_summary_{}.tsv".format(title)), sep='\t')
    for cur in set(homogeneities.index).intersection(constants.ALGOS):
        ax.plot(cutoffs, homogeneities.loc[cur], c=constants.COLORDICT[cur])

    ax.set_xlabel("similarity cutoff", fontsize=22)
    ax.set_ylabel("mean intra-module homogeneity", fontsize=22)
    ax.set_title(title, fontsize=22)

    patches = [Line2D([0], [0], marker='o', color='gray', label=constants.ALGOS_ACRONYM[a], markersize=12,
                      markerfacecolor=constants.COLORDICT[a], alpha=0.7) for a in algos]
    ax.legend(handles=patches, loc=(0,1.1), ncol=3, fontsize=22, facecolor='#ffffff')


if __name__=="__main__":

    base_folder=os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation")
    homogeneity_file_format = "homogeneity_avg_matrix_{}_{}.tsv"
    cutoffs = [1.0, 2.0, 3.0, 4.0]

    fig,axs=plt.subplots(1,2, figsize=(22,10))
    fig_violin, axs_violin = plt.subplots(2, len(cutoffs), figsize=(4 * len(cutoffs) * 2, 10))

    datasets = ["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"]
    algos = ['DOMINO4', 'netbox2_string']
    prefix = "GE"
    variability_plot(algos, datasets, base_folder, homogeneity_file_format, prefix, cutoffs, axs[0], axs_violin[0], title="GE")

    datasets = ["brca", "crh", "scz", "tri", "t2d", "cad", "bmd", "hgt", "amd", "af"]
    algos = ['DOMINO4', 'netbox2_string']
    prefix = "PASCAL_SUM"
    variability_plot(algos, datasets, base_folder, homogeneity_file_format, prefix, cutoffs, axs[1], axs_violin[1], title="GWAS")

    fig.text(0.01,0.97, "A:", weight='bold',fontsize=22)
    fig.text(0.5, 0.97, "B:", weight='bold',fontsize=22)
    fig.tight_layout()
    fig.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "figure_15.png"))

