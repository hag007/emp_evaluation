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

ENABLE_GO_GRAPH = False
IS_GO_GRAPH_ONLY=False
GO_PCA=False
TERMS_SIMILARITY_TO_NUM_OF_TERMS=True
RATIO_TO_GO_TERM = True
QVAL_TH = 0.01
SIM_TH= 0.4

import seaborn as sns


def plot_grid(df_all_hg_pval, df_all_hg_size, grid_type):
    df_all_hg_pval=df_all_hg_pval.loc[np.sort(df_all_hg_pval.index), np.sort(df_all_hg_pval.columns)]
    df_all_hg_size=df_all_hg_size.loc[np.sort(df_all_hg_size.index), np.sort(df_all_hg_size.columns)]


    fig, ax = plt.subplots(figsize=(15, 15))
    ax.set_facecolor('#fffde3')
    ax.grid(color='gray')
    x = np.arange(len(df_all_hg_pval.columns))
    y = np.arange(len(df_all_hg_pval.index))
    xx, yy = zip(*[(a, b) for a in x for b in y]) # if df_all_hg_pval.iloc[b, a] != 0
    c = [df_all_hg_pval.iloc[b, a] for a in x for b in y] #  if df_all_hg_pval.iloc[b, a] != 0
    s = [df_all_hg_size.iloc[b, a] for a in x for b in y ] # if df_all_hg_pval.iloc[b, a] != 0
    s = np.array(s)
    # im = plt.imshow(np.array(c).reshape(len(c), 1), cmap='bwr')
    # im.remove()
    sc = ax.scatter(xx, yy, s / float(max(s)) * 3000, c=c, cmap='bwr', vmin=np.percentile(c, 10),
                    vmax=np.percentile(c, 90))
    for s_i, a in enumerate(s):
        ax.annotate( "{}\n({})".format(int(a), round(c[s_i],2)), (xx[s_i], yy[s_i]), color='green', size=20)
    ax.legend(loc='upper left')
    ax.margins(0.03, 0.03)
    # ax.locator_params(axis='x', nbins=len(df_all_hg_pval.columns))
    # ax.locator_params(axis='y', nbins=len(df_all_hg_pval.index))
    ax.set_xlabel("datasets", fontdict={"size" : 25})
    plt.subplots_adjust(left=0.25, right=0.99, top=0.99, bottom=0.05)
    plt.xticks(np.arange(len(df_all_hg_pval.columns)), tuple([x[3:] if x.startswith("GE_") else x for x in list(df_all_hg_pval.columns.values)]), rotation='vertical', size=25)
    ax.set_ylabel("algos", fontdict={"size" : 25})
    plt.yticks(np.arange(len(df_all_hg_pval.index)), tuple([x[3:] if x.startswith("GE_") else x for x in list(df_all_hg_pval.index.values)]), size=25)
    ax_ = plt.gca()
    aspect = 20
    pad_fraction = 0.5
    divider = make_axes_locatable(ax_)
    width = axes_size.AxesY(ax_, aspect=1. / aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=0.3, pad=0.4)
    plt.colorbar(mappable=sc, cax=cax)
    plt.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "datasets_algo_grid_{}.png".format(grid_type)))
    plt.clf()


def plot_scatter(df_all_hg_pval, df_all_hg_size, grid_type):
    df_all_hg_pval=df_all_hg_pval.loc[np.sort(df_all_hg_pval.index), np.sort(df_all_hg_pval.columns)]
    df_all_hg_size=df_all_hg_size.loc[np.sort(df_all_hg_size.index), np.sort(df_all_hg_size.columns)]
    algos = list(df_all_hg_pval.index)

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.grid(color='gray')
    i_x = np.arange(len(df_all_hg_pval.columns))
    i_y = np.arange(len(df_all_hg_pval.index))
    x = np.array([df_all_hg_pval.iloc[b, a] for a in i_x for b in i_y])
    y = np.array([df_all_hg_size.iloc[b, a] for a in i_x for b in i_y ])
    labels=np.array([b for a in i_x for b in i_y])
    sc = ax.scatter(x, y, c=labels/float(len(algos)-1), cmap='jet')

    colormap = cm.jet

    colorlist = [ml_colors.rgb2hex(colormap(i)) for i in
                 np.array(list(range(len(algos)))) / float(len(algos) - 1)]
    patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                      markerfacecolor=c) for i, a, c in zip(list(range(len(algos))), algos, colorlist)]
    ax.legend(handles=patches, loc='upper left')
    ax.set_xlabel("X label")
    ax.set_ylabel("Y label")
    plt.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "datasets_algo_scatter_{}.png".format(grid_type)))
    plt.clf()

def plot_scatter3d(df_size, df_ratio, df_variability, grid_type):
    df_size=df_size.loc[np.sort(df_size.index), np.sort(df_size.columns)]
    df_ratio=df_ratio.loc[np.sort(df_ratio.index), np.sort(df_ratio.columns)]
    df_variability = df_variability.loc[np.sort(df_variability.index), np.sort(df_variability.columns)]
    algos = list(df_size.index)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.grid(color='gray')
    i_x = np.arange(len(df_size.columns))
    i_y = np.arange(len(df_size.index))
    x = np.array([df_size.iloc[b, a] for a in i_x for b in i_y])
    y = np.array([df_ratio.iloc[b, a] for a in i_x for b in i_y])
    z = np.array([df_variability.iloc[b, a] for a in i_x for b in i_y])
    labels=np.array([b for a in i_x for b in i_y])
    sc = ax.scatter(x, y, z,  c=labels/float(len(algos)-1), cmap='jet')

    colormap = cm.jet

    colorlist = [ml_colors.rgb2hex(colormap(i)) for i in
                 np.array(list(range(len(algos)))) / float(len(algos) - 1)]
    patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                      markerfacecolor=c) for i, a, c in zip(list(range(len(algos))), algos, colorlist)]
    ax.legend(handles=patches, loc='upper left')
    # ax.legend()
    ax.margins(0.03, 0.03)
    # ax.locator_params(axis='x', nbins=len(df_all_hg_pval.columns))
    # ax.locator_params(axis='y', nbins=len(df_all_hg_pval.index))
    ax.set_xlabel("X label")
    plt.subplots_adjust(left=0.25, right=0.99, top=0.99, bottom=0.05)
    ax.set_ylabel("Y label")
    ax_ = plt.gca()
    aspect = 20
    pad_fraction = 0.5
    divider = make_axes_locatable(ax_)
    width = axes_size.AxesY(ax_, aspect=1. / aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=0.3, pad=0.4)
    plt.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "datasets_algo_scatter_{}.png".format(grid_type)))
    plt.clf()


def plot_scatter_with_size(df_size, df_ratio, df_variability, grid_type, ax=None):
    df_size=df_size.loc[np.sort(df_size.index), np.sort(df_size.columns)]
    df_ratio=df_ratio.loc[np.sort(df_ratio.index), np.sort(df_ratio.columns)]
    df_variability=df_variability.loc[np.sort(df_variability.index), np.sort(df_variability.columns)]
    algos = list(df_size.index)
    sns.set_palette("husl", len(algos))

    if ax is None:
        fig, ax = plt.subplots(figsize=(13, 13))

    ax.set_facecolor('#fffde3')
    ax.grid(color='gray')

    i_x = np.arange(len(df_size.columns))
    i_y = np.arange(len(df_size.index))
    size = np.array([df_size.iloc[b, a] for a in i_x for b in i_y])
    x = np.array([df_ratio.iloc[b, a] for a in i_x for b in i_y])
    y = np.array([df_variability.iloc[b, a] for a in i_x for b in i_y])
    labels=np.array([sns.color_palette()[b] for a in i_x for b in i_y])
    sc = ax.scatter(x, y, s=size/float(max(size))*2000+20, c=labels, cmap='jet', alpha=0.7)

    sorted_list = sorted([[x[i], y[i]] for i in range(len(x))], reverse=True)
    pareto_front = [sorted_list[0]]
    for pair in sorted_list[1:]:
        if True:
            if pair[1] >= pareto_front[-1][1]:
                pareto_front.append(pair)
        else:
            if pair[1] <= pareto_front[-1][1]:
                pareto_front.append(pair)

    pf_X = [pair[0] for pair in pareto_front]
    pf_Y = [pair[1] for pair in pareto_front]
    ax.plot(pf_X, pf_Y)


    colormap = cm.jet

    colorlist = [sns.color_palette()[i] for i in
                 np.array(list(range(len(algos))))]
    patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                      markerfacecolor=c) for i, a, c in zip(list(range(len(algos))), algos, colorlist)]
    ax.legend(handles=patches, loc='lower left', prop={'size': 12})
    # ax.legend()
    ax.margins(0.03, 0.03)
    # ax.locator_params(axis='x', nbins=len(df_all_hg_pval.columns))
    # ax.locator_params(axis='y', nbins=len(df_all_hg_pval.index))
    ax.set_xlabel("EHR", fontdict={"size":12})
    plt.subplots_adjust(left=0.25, right=0.99, top=0.99, bottom=0.05)
    ax.set_ylabel("heterogeneity", fontdict={"size":12})
    ax.set_title("Pareto Fontier")
    plt.xticks(size=12)
    plt.yticks(size=12)
    plt.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "datasets_algo_scatter_{}.png".format(grid_type)))
    # plt.clf()


def main():
    main_path = "/home/hag007/Desktop/aggregate_report/visual"
    df_counts=pd.read_csv(os.path.join( main_path, "true_positive_counts.tsv"), sep='\t', index_col=0)
    df_ratio=pd.read_csv(os.path.join(main_path, "true_positive_ratio.tsv"), sep='\t', index_col=0)
    plot_grid(df_ratio, df_counts, "true_positive")
    plot_scatter(df_ratio, df_counts, "true_positive")
    df_counts = pd.read_csv(os.path.join(main_path, "empirical_terms_counts.tsv"), sep='\t', index_col=0)
    df_variability = pd.read_csv(os.path.join(main_path, "empirical_terms_variability.tsv"), sep='\t', index_col=0)
    df_variability[df_counts == 0]=0
    df_variability += np.abs(np.min(df_variability.values))+1
    df_variability[df_counts==0]=0
    plot_grid(df_variability, df_counts, "empirical_terms")
    plot_scatter(df_variability, df_counts, "empirical_terms")

    plot_scatter3d(df_ratio, df_counts, df_variability, "3d")
    plot_scatter_with_size(df_counts, df_ratio, df_variability, "size")

if __name__=="__main__":
    main()