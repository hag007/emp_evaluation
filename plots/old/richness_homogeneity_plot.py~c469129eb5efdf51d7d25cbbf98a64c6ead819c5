import pandas as pd

from fastsemsim.SemSim import *

import matplotlib.pyplot as plt

from rpy2.robjects import pandas2ri
pandas2ri.activate()
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





def plot_richness_homogeneity(df_richness, df_homogeneity, grid_type, ax=None, title="", algos=constants.ALGOS_ACRONYM.keys()):
    df_richness=df_richness.loc[np.sort(df_richness.index), np.sort(df_richness.columns)].loc[algos]
    df_homogeneity=df_homogeneity.loc[np.sort(df_homogeneity.index), np.sort(df_homogeneity.columns)].loc[algos]

    if ax is None:
        fig, ax = plt.subplots(figsize=(13, 13))

    ax.set_facecolor('#ffffff')
    ax.grid(color='gray')

    i_x = list(df_richness.columns)
    i_y = list(df_richness.index)
    size = np.array([df_richness.loc[b, a] for a in i_x for b in i_y])
    x = np.array([df_homogeneity.loc[b, a] for a in i_x for b in i_y])
    y = np.array([df_richness.loc[b, a] for a in i_x for b in i_y])
    labels=np.array([constants.COLORDICT[b] for a in i_x for b in i_y])
    sc = ax.scatter(x, y, s=200, c=labels, cmap='jet') # size/float(max(size))*2000+20

    # sorted_list = sorted([[x[i], y[i]] for i in range(len(x))], reverse=True)
    # pareto_front = [sorted_list[0]]
    # for pair in sorted_list[1:]:
    #     if True:
    #         if pair[1] >= pareto_front[-1][1]:
    #             pareto_front.append(pair)
    #     else:
    #         if pair[1] <= pareto_front[-1][1]:
    #             pareto_front.append(pair)
    #
    # pf_X = [pair[0] for pair in pareto_front]
    # pf_Y = [pair[1] for pair in pareto_front]
    # ax.plot(pf_X, pf_Y)

    patches= [Line2D([0], [0], marker='o', color='gray', label=constants.ALGOS_ACRONYM[a], markersize=12,
                        markerfacecolor=constants.COLORDICT[a], alpha=0.7) for a in algos]
    ax.legend(handles=patches, loc=(0,1.1), ncol=2, prop={'size': 20}, facecolor='#ffffff')
    # ax.margins(0.03, 0.03)
    ax.set_xlabel("homogeneity", fontdict={"size":22})
    ax.set_ylabel("richness", fontdict={"size": 22})
    ax.set_title(title, fontdict={"size": 22})
    plt.subplots_adjust(left=0.25, right=0.99, top=0.99, bottom=0.05)


def plot_richness_homogeneity_agg(df_richness, df_homogeneity, grid_type, ax=None, title="", algos=constants.ALGOS_ACRONYM.keys()):
    df_richness=df_richness.loc[np.sort(df_richness.index), np.sort(df_richness.columns)].loc[algos]
    df_homogeneity=df_homogeneity.loc[np.sort(df_homogeneity.index), np.sort(df_homogeneity.columns)].loc[algos]

    if ax is None:
        fig, ax = plt.subplots(figsize=(13, 13))

    ax.set_facecolor('#ffffff')
    ax.grid(color='gray')

    for cur_algo in df_richness.index:
        y=df_richness.loc[cur_algo,:].median()
        x=df_homogeneity.loc[cur_algo,:].mean()
        sc = ax.scatter(x, y, s=200, c=constants.COLORDICT[cur_algo])

    patches = [Line2D([0], [0], marker='o', color='gray', label=constants.ALGOS_ACRONYM[a], markersize=12,
                      markerfacecolor=constants.COLORDICT[a], alpha=0.7) for a in algos]
    ax.legend(handles=patches, loc=(0,1.1), ncol=2, prop={'size': 20}, facecolor='#ffffff')
    # ax.margins(0.03, 0.03)
    ax.set_xlabel("homogeneity", fontdict={"size":22})
    ax.set_ylabel("richness", fontdict={"size": 22})
    ax.set_title(title, fontdict={"size": 22})
    plt.subplots_adjust(left=0.25, right=0.99, top=0.99, bottom=0.05)



def main():

    cutoffs=[1.0,2.0,3.0,4.0]

    for cutoff in cutoffs:


        fig_1, axs_1 = plt.subplots(1,2, figsize=(23, 12))
        fig_2, axs_2 = plt.subplots(1, 2, figsize=(23, 12))

        suffix="GE_{}".format(cutoff)
        algos = ["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox", "keypathwayminer_INES_GREEDY"]
        df_richness=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "solution_richness_matrix_{}.tsv".format(suffix)), sep='\t', index_col=0).loc[np.sort(algos),:]
        df_homogeneity = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "homogeneity_avg_matrix_{}.tsv".format(suffix)), sep='\t', index_col=0).loc[np.sort(algos),:]
        plot_richness_homogeneity(df_richness, df_homogeneity, "size", ax=axs_1[0], title="GE", algos=algos)
        plot_richness_homogeneity_agg(df_richness, df_homogeneity, "size", ax=axs_2[0], title="GE", algos=algos)


        suffix="PASCAL_SUM_{}".format(cutoff)
        algos = ["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox"]
        df_richness=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "solution_richness_matrix_{}.tsv".format(suffix)), sep='\t', index_col=0).loc[np.sort(algos),:]
        df_homogeneity = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "homogeneity_avg_matrix_{}.tsv".format(suffix)), sep='\t', index_col=0).loc[np.sort(algos),:]
        plot_richness_homogeneity(df_richness, df_homogeneity, "size", ax=axs_1[1], title="GWAS", algos=algos)
        plot_richness_homogeneity_agg(df_richness, df_homogeneity, "size", ax=axs_2[1], title="GWAS", algos=algos)

        fig_1.text(0.01,0.97, "A:", weight='bold',fontsize=22)
        fig_1.text(0.5, 0.97, "B:", weight='bold',fontsize=22)
        fig_1.tight_layout()
        fig_1.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_16_{}.png".format(cutoff)))

        fig_2.text(0.01,0.97, "A:", weight='bold',fontsize=22)
        fig_2.text(0.5, 0.97, "B:", weight='bold',fontsize=22)
        fig_2.tight_layout()
        fig_2.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_17_{}.png".format(cutoff)))


if __name__=="__main__":
    main()