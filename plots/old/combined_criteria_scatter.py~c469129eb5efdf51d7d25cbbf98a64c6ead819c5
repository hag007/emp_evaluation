import pandas as pd

from fastsemsim.SemSim import *

import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
import constants

from scipy.cluster import hierarchy
import scipy.cluster.hierarchy as hcl
from scipy.spatial.distance import squareform


import utils.go
import utils.go_hierarcies
import math
import random
import matplotlib.cm as cm

from matplotlib.lines import Line2D
import matplotlib.colors as ml_colors

from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size

import seaborn as sns



ENABLE_GO_GRAPH = False
IS_GO_GRAPH_ONLY=False
GO_PCA=False
TERMS_SIMILARITY_TO_NUM_OF_TERMS=True
RATIO_TO_GO_TERM = True
QVAL_TH = 0.01
SIM_TH= 0.4


def plot_agg_scatter(df_measurements_x, df_measurements_y, ax, title, x_label, y_label, x_agg_metric, y_agg_metric):



    # if ax is None:
    #     fig, ax = plt.subplots(figsize=(13, 13))

    ax.set_facecolor('#ffffff')
    ax.grid(color='gray')
    df_measurements_x=df_measurements_x.dropna(how='all')
    df_measurements_y=df_measurements_y.dropna(how='all')
    datasets = df_measurements_x.columns
    algos = set(df_measurements_x.dropna().index).intersection(df_measurements_y.index)
    for cur_algo in algos:
        x= getattr(df_measurements_x.loc[cur_algo, :], x_agg_metric)()
        y= getattr(df_measurements_y.loc[cur_algo, :], y_agg_metric)()
        sc = ax.scatter(x, y, s=400, c=constants.COLORDICT[cur_algo])

    patches = [Line2D([0], [0], marker='o', markersize=12, color='gray', label=constants.ALGOS_ACRONYM[a],
                      markerfacecolor=constants.COLORDICT[a]) for i, a in zip(list(range(len(algos))), algos)]
    ax.legend(handles=patches, loc=(0.0,1.1), prop={'size': 20},ncol=2, facecolor='#ffffff')
    # ax.margins(0.03, 0.03)
    ax.set_xlabel(x_label, fontdict={"size":22})
    ax.set_ylabel(y_label, fontdict={"size": 22})
    ax.set_title(title, fontdict={"size":22})
    plt.subplots_adjust(left=0.25, right=0.99, top=0.99, bottom=0.05)




def plot_fig(file_name_x_ge, file_name_y_ge, file_name_x_gwas, file_name_y_gwas, title, x_label, y_label, fig_suffix, x_agg_metric="mean", y_agg_metric="mean", df_total=pd.DataFrame()):

    fig, axs = plt.subplots(1,2, figsize=(25, 12))

    plot_single_ax(axs[0], file_name_x_ge, file_name_y_ge, title="{}, GE".format(title), x_label=x_label, y_label=y_label, x_agg_metric=x_agg_metric, y_agg_metric=y_agg_metric, df_total=df_total)
    plot_single_ax(axs[1], file_name_x_gwas, file_name_y_gwas, title="{}, GWAS".format(title), x_label=x_label, y_label=y_label, x_agg_metric=x_agg_metric, y_agg_metric=y_agg_metric, df_total=df_total)

    fig.text(0.01, 0.97, "A:", weight='bold', fontsize=22)
    fig.text(0.5, 0.97, "B:", weight='bold', fontsize=22)
    fig.tight_layout()
    fig.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_agg_{}.png".format(fig_suffix)))


def plot_given_fig(file_name_x_ge, file_name_y_ge, file_name_x_gwas, file_name_y_gwas, title, x_label, y_label, fig_suffix, axs, x_agg_metric="mean", y_agg_metric="mean", df_total=pd.DataFrame()):


    plot_single_ax(axs[0], file_name_x_ge, file_name_y_ge, title="{}, GE".format(title), x_label=x_label, y_label=y_label, x_agg_metric=x_agg_metric, y_agg_metric=y_agg_metric, df_total=df_total)
    plot_single_ax(axs[1], file_name_x_gwas, file_name_y_gwas, title="{}, GWAS".format(title), x_label=x_label, y_label=y_label, x_agg_metric=x_agg_metric, y_agg_metric=y_agg_metric, df_total=df_total)


def plot_single_ax(ax, file_name_x, file_name_y, title, x_label, y_label, x_agg_metric, y_agg_metric, df_total):
    df_measurements_x = pd.read_csv(file_name_x, sep='\t',
                                    index_col=0)  # .loc[np.sort(constants.ALGOS_ACRONYM.keys()),:]
    df_measurements_y = pd.read_csv(file_name_y, sep='\t',
                                    index_col=0)  # .loc[np.sort(constants.ALGOS_ACRONYM.keys()),:

    df_measurements_x[df_measurements_x.isna()]=0
    df_measurements_y[df_measurements_y.isna()]=0

    if title=="GWAS":
        df_measurements_x=df_measurements_x.drop(labels=['Alzheimer.G50'],axis=1,errors='ignore')
        df_measurements_y=df_measurements_y.drop(labels=['Alzheimer.G50'], axis=1,errors='ignore')

    df_measurements_x =df_measurements_x.loc[set(constants.ALGOS).intersection(df_measurements_x.index).intersection(df_measurements_y.index)].dropna()
    df_measurements_y =df_measurements_y.loc[set(constants.ALGOS).intersection(df_measurements_x.index).intersection(df_measurements_y.index)]

    df_measurements_x=df_measurements_x.loc[np.sort(df_measurements_x.index), np.sort(df_measurements_x.columns)]
    df_measurements_y=df_measurements_y.loc[np.sort(df_measurements_y.index), np.sort(df_measurements_y.columns)]
    df_total.loc[:, "{}_{}_{}".format(x_label, x_agg_metric, title.split(',')[1].strip())] = getattr(df_measurements_x.T, x_agg_metric)()
    df_total.loc[:, "{}_{}_{}".format(x_label, "std", title.split(',')[1].strip())] = df_measurements_x.T.std()
    df_total.loc[:, "{}_{}_{}".format(y_label, y_agg_metric, title.split(',')[1].strip())] = getattr(df_measurements_y.T, y_agg_metric)()
    df_total.loc[:, "{}_{}_{}".format(y_label, "std", title.split(',')[1].strip())] = df_measurements_y.T.std()
    plot_agg_scatter(df_measurements_x, df_measurements_y, ax=ax, title=title, x_label=x_label, y_label=y_label, x_agg_metric=x_agg_metric, y_agg_metric=y_agg_metric)


def main():
    df_total=pd.DataFrame(index=constants.ALGOS)
    main_folder= os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation")


    file_name_x_ge = "ehr_matrix_GE.tsv"
    file_name_y_ge = "../mehr_cache_files/summary_mEHR_mean_10_GE.tsv"
    file_name_x_gwas = "ehr_matrix_PASCAL_SUM.tsv"
    file_name_y_gwas = "../mehr_cache_files/summary_mEHR_mean_10_PASCAL_SUM.tsv" # ""recovery_results_PASCAL_SUM_100_0.3.tsv"
    x_label = "EHR"
    y_label = "mEHR"
    title = "{} vs {}".format(x_label,y_label)
    fig_suffix=1
    cutoff =3.0
    plot_fig(os.path.join(main_folder, file_name_x_ge), os.path.join(main_folder, file_name_y_ge),
             os.path.join(main_folder, file_name_x_gwas), os.path.join(main_folder, file_name_y_gwas),
             title, x_label, y_label, fig_suffix, df_total=df_total)


    file_name_x_ge = "ehr_matrix_GE.tsv"
    file_name_y_ge = "robustness_GE_100_0.2_matrix_f1.tsv"
    file_name_x_gwas = "ehr_matrix_PASCAL_SUM.tsv"
    file_name_y_gwas = "robustness_PASCAL_SUM_100_0.2_matrix_f1.tsv" # ""recovery_results_PASCAL_SUM_100_0.3.tsv"
    x_label = "EHR"
    y_label = "Robustness (F1)"
    title = "{} vs {}".format(x_label,y_label)
    fig_suffix=2
    fig, axs = plt.subplots(2, 2, figsize=(25, 25))

    plot_given_fig(os.path.join(main_folder, file_name_x_ge), os.path.join(main_folder, file_name_y_ge),
             os.path.join(main_folder, file_name_x_gwas), os.path.join(main_folder, file_name_y_gwas),
             title, x_label, y_label, fig_suffix, axs[0], df_total=df_total)


    # file_name_x_ge = "ratio_matrix.tsv"
    # file_name_y_ge = "count_matrix.tsv"
    # file_name_x_gwas = "ratio_matrix.tsv"
    # file_name_y_gwas = "count_matrix.tsv" # ""recovery_results_PASCAL_SUM_100_0.3.tsv"
    # x_label = "EHR"
    # y_label = "term count"
    # title = "{} vs {}".format(x_label,y_label)
    # fig_suffix=2
    # fig, axs = plt.subplots(2, 2, figsize=(25, 25))
    #
    # plot_given_fig(os.path.join(main_ge_folder, file_name_x_ge), os.path.join(main_ge_folder, file_name_y_ge),
    #          os.path.join(main_gwas_folder, file_name_x_gwas), os.path.join(main_gwas_folder, file_name_y_gwas),
    #          title, x_label, y_label, fig_suffix, axs[0], df_total=df_total)


    file_name_y_ge = "robustness_auc_GE_100_0.2.tsv"
    file_name_y_gwas = "robustness_auc_PASCAL_SUM_100_0.2.tsv" # ""recovery_results_PASCAL_SUM_100_0.3.tsv"
    y_label = "Robustness (AUPR)"
    title = "{} vs {}".format(x_label,y_label)

    plot_given_fig(os.path.join(main_folder, file_name_x_ge), os.path.join(main_folder, file_name_y_ge),
                   os.path.join(main_folder, file_name_x_gwas), os.path.join(main_folder, file_name_y_gwas),
                   title, x_label, y_label, fig_suffix, axs[1], df_total=df_total)


    fig.text(0.01, 0.97, "A:", weight='bold', fontsize=22)
    fig.text(0.5, 0.97, "B:", weight='bold', fontsize=22)
    fig.text(0.01, 0.5, "C:", weight='bold', fontsize=22)
    fig.text(0.5, 0.5, "D:", weight='bold', fontsize=22)
    fig.tight_layout()
    fig.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_agg_{}.png".format(fig_suffix)))


    file_name_x_ge = "ehr_matrix_GE.tsv"
    file_name_y_ge = "richness_matrix_GE_{}.tsv".format(cutoff)
    file_name_x_gwas = "ehr_matrix_PASCAL_SUM.tsv"
    file_name_y_gwas = "richness_matrix_PASCAL_SUM_{}.tsv".format(cutoff) # ""recovery_results_PASCAL_SUM_100_0.3.tsv"
    x_label = "EHR"
    y_label = "Biological Richness"
    title = "{} vs {}".format(x_label,y_label)
    fig_suffix=3

    plot_fig(os.path.join(main_folder, file_name_x_ge), os.path.join(main_folder, file_name_y_ge),
             os.path.join(main_folder, file_name_x_gwas), os.path.join(main_folder, file_name_y_gwas),
             title, x_label, y_label, fig_suffix, x_agg_metric="mean", y_agg_metric="median", df_total=df_total)

    file_name_x_ge = "richness_matrix_GE_{}.tsv".format(cutoff)
    file_name_y_ge = "homogeneity_avg_matrix_GE_{}.tsv".format(cutoff)
    file_name_x_gwas = "richness_matrix_PASCAL_SUM_{}.tsv".format(cutoff)
    file_name_y_gwas = "homogeneity_avg_matrix_PASCAL_SUM_{}.tsv".format(cutoff) # ""recovery_results_PASCAL_SUM_100_0.3.tsv"
    x_label = "Biological Richness"
    y_label = "Intra-Module Homomgeneity"
    title = "{} vs {}".format(x_label,y_label)
    fig_suffix=4

    plot_fig(os.path.join(main_folder, file_name_x_ge), os.path.join(main_folder, file_name_y_ge),
             os.path.join(main_folder, file_name_x_gwas), os.path.join(main_folder, file_name_y_gwas),
             title, x_label, y_label, fig_suffix, x_agg_metric="median", y_agg_metric="mean", df_total=df_total)


    file_name_x_ge = "../mehr_cache_files/summary_mEHR_mean_10_GE.tsv"
    file_name_y_ge = "homogeneity_avg_matrix_GE_{}.tsv".format(cutoff)
    file_name_x_gwas = "../mehr_cache_files/summary_mEHR_mean_10_PASCAL_SUM.tsv"
    file_name_y_gwas = "homogeneity_avg_matrix_PASCAL_SUM_{}.tsv".format(cutoff) # ""recovery_results_PASCAL_SUM_100_0.3.tsv"
    x_label = "mEHR"
    y_label = "Intra-Module Homomgeneity"
    title = "{} vs {}".format(x_label,y_label)
    fig_suffix=5

    plot_fig(os.path.join(main_folder, file_name_x_ge), os.path.join(main_folder, file_name_y_ge),
             os.path.join(main_folder, file_name_x_gwas), os.path.join(main_folder, file_name_y_gwas),
             title, x_label, y_label, fig_suffix, df_total=df_total)

    omic="GE"
    algos=["jactivemodules_greedy", "jactivemodules_sa", "netbox", "bionet", "keypathwayminer_INES_GREEDY", "DOMINO"]
    agg_type=["mean", "median"]
    df_total_tmp=df_total.loc[algos,[a for a in df_total.columns if omic in a and any([b in a for b in agg_type])]]
    df_total_tmp=df_total_tmp.rename(columns={a:a.split("_")[0] for a in df_total_tmp}, index={a: constants.ALGOS_ACRONYM[a] for a in algos})
    df_total_tmp=df_total_tmp.apply(lambda a: a.apply(lambda b: '%.2E' % b))
    df_total_tmp.to_csv(os.path.join(main_folder, "criteria_summary_{}_{}_{}.tsv".format(agg_type[0],omic, cutoff)), sep='\t')
    agg_type=["std"]
    df_total_tmp=df_total.loc[algos,[a for a in df_total.columns if omic in a and any([b in a for b in agg_type])]]
    df_total_tmp=df_total_tmp.rename(columns={a:a.split("_")[0] for a in df_total_tmp}, index={a: constants.ALGOS_ACRONYM[a] for a in algos})
    df_total_tmp=df_total_tmp.apply(lambda a: a.apply(lambda b: '%.2E' % b))
    df_total_tmp.to_csv(os.path.join(main_folder, "criteria_summary_{}_{}_{}.tsv".format(agg_type[0],omic, cutoff)), sep='\t')


    omic="GWAS"
    algos=["jactivemodules_greedy", "jactivemodules_sa", "netbox", "bionet", "DOMINO"]
    agg_type=["mean", "median"]
    df_total_tmp=df_total.loc[algos,[a for a in df_total.columns if omic in a and any([b in a for b in agg_type])]]
    df_total_tmp=df_total_tmp.rename(columns={a:a.split("_")[0] for a in df_total_tmp}, index={a: constants.ALGOS_ACRONYM[a] for a in algos})
    df_total_tmp=df_total_tmp.apply(lambda a: a.apply(lambda b: '%.2E' % b))
    df_total_tmp.to_csv(os.path.join(main_folder, "criteria_summary_{}_{}_{}.tsv".format(agg_type[0],omic, cutoff)), sep='\t')
    agg_type=["std"]
    df_total_tmp=df_total.loc[algos,[a for a in df_total.columns if omic in a and any([b in a for b in agg_type])]]
    df_total_tmp=df_total_tmp.rename(columns={a:a.split("_")[0] for a in df_total_tmp}, index={a: constants.ALGOS_ACRONYM[a] for a in algos})
    df_total_tmp=df_total_tmp.apply(lambda a: a.apply(lambda b: '%.2E' % b))
    df_total_tmp.to_csv(os.path.join(main_folder, "criteria_summary_{}_{}_{}.tsv".format(agg_type[0],omic, cutoff)), sep='\t')


    df_total.to_csv(os.path.join(main_folder, "aggergated_criteria_summary_{}.tsv".format(cutoff)), sep='\t')

    df_total_ge=df_total.loc[:,[a for a in df_total.columns if "GE" in a]]
    df_total_ge_agg=df_total_ge.loc[:,[a for a in df_total_ge.columns if "mean" in a or "median" in a]]
    df_total_ge_agg=df_total_ge_agg.rename(columns={a: a.split("_")[0] for a in df_total_ge_agg.columns})
    df_total_ge_agg.to_csv(os.path.join(main_folder, "aggergated_criteria_summary_{}_{}_{}.tsv".format("ge","agg", cutoff)), sep='\t')

    df_total_ge_std=df_total_ge.loc[:,[a for a in df_total_ge.columns if "std" in a]]
    df_total_ge_std=df_total_ge_std.rename(columns={a: a.split("_")[0] for a in df_total_ge_std.columns})
    df_total_ge_std.to_csv(os.path.join(main_folder, "aggergated_criteria_summary_{}_{}_{}.tsv".format("ge","std", cutoff)), sep='\t')

    df_total_gwas=df_total.loc[:,[a for a in df_total.columns if "GWAS" in a]]
    df_total_gwas_agg=df_total_gwas.loc[:,[a for a in df_total_gwas.columns if "mean" in a or "median" in a]]
    df_total_gwas_agg=df_total_gwas_agg.rename(columns={a: a.split("_")[0] for a in df_total_gwas_agg.columns})
    df_total_gwas_agg.to_csv(os.path.join(main_folder, "aggergated_criteria_summary_{}_{}_{}.tsv".format("gwas","agg", cutoff)), sep='\t')
    df_total_gwas_std=df_total_gwas.loc[:,[a for a in df_total_gwas.columns if "std" in a]]
    df_total_gwas_std=df_total_gwas_std.rename(columns={a: a.split("_")[0] for a in df_total_gwas_std.columns})
    df_total_gwas_std.to_csv(os.path.join(main_folder, "aggergated_criteria_summary_{}_{}_{}.tsv".format("gwas","std", cutoff)), sep='\t')


if __name__=="__main__":
    main()