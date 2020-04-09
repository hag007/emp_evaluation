import sys
sys.path.insert(0, '../')

import pandas as pd

from fastsemsim.SemSim import *
from fastsemsim.Ontology import ontologies
from fastsemsim.Ontology import AnnotationCorpus
from fastsemsim.SemSim.SetSemSim import SetSemSim
import matplotlib
matplotlib.use("Agg")

import constants

import multiprocessing

from utils.daemon_multiprocessing import func_star

import argparse

import math
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt

import scipy.spatial.distance as ssd





METHOD='weighted'

def solutions_hierarchical_clustering(color_strategy=None, color_strategy_title="none", ax=None, distance_file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR, "cancer_type_go_distance.tsv")):


    df_summary=pd.read_csv(distance_file_name, sep='\t', index_col=0)
    # df_summary-=1

    df_summary+=np.abs(np.nanmin(df_summary.values))+1
    df_summary = df_summary.where(~np.isnan(df_summary.values), other=1)
    for cur in df_summary.columns:
        df_summary.loc[cur,cur]=0

    distArray = ssd.squareform(df_summary.values )

    linked = linkage(distArray, method=METHOD, metric='euclidean')


    if ax is None:
        plt.figure(figsize=(10, 7))
        ax = plt.gca()

    lbl=filter(lambda f: f is not None, [a[:-len(b)]+ constants.ALGOS_ACRONYM[b] if b in a else None for a in df_summary.index.values for b in constants.ALGOS_ACRONYM.keys()])

    dendrogram(linked,
               color_threshold=0,
               orientation='top',
               labels=lbl,
               distance_sort='descending',
               show_leaf_counts=True,
               leaf_rotation=90,
               leaf_font_size=20,
               ax=ax)

    ax.set_title('colored by {}'.format(color_strategy_title))
    ax.set_ylabel("distance")
    ax.set_xlabel("solutions")


    # Apply the right color to each label
    if color_strategy is not None:
        my_color=["red", "green", "blue", "black", "brown", "orange", "purple"]
        # colors = [[i_y for i_y, y in enumerate(algos) if y in x][0] for x in list(df_summary.index)]

        xlbls = ax.get_xmajorticklabels()
        num = -1
        for lbl in xlbls:
            val = my_color[[i_y for i_y, y in enumerate(color_strategy) if y in lbl._text][0]]
            lbl.set_color(val)

    plt.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "hierarchical_go_clustering.png"))



def algos_hierarchical_clustering(datasets, algos, distance_file_name, ax=None):


    df_summary=pd.read_csv(distance_file_name, sep='\t', index_col=0)
    # df_summary-=1

    df_summary+=np.abs(np.nanmin(df_summary.values))+1
    df_summary = df_summary.where(~np.isnan(df_summary.values), other=1)
    for cur in df_summary.columns:
        df_summary.loc[cur,cur]=0

    distArray = ssd.squareform(df_summary.values )


    df_summary3d = np.empty((len(algos), len(algos), len(datasets)),dtype=np.object)
    for cur_col in df_summary.columns:
        for cur_index in df_summary.index:
            cur_ds_col = [a for a in datasets if a in cur_col][0]
            cur_algo_col = [a for a in algos if a in cur_col][0]
            cur_ds_index = [a for a in datasets if a in cur_index][0]
            cur_algo_index = [a for a in algos if a in cur_index][0]
            if cur_ds_col != cur_ds_index: continue

            df_summary3d[algos.index(cur_algo_index), algos.index(cur_algo_col)
            , datasets.index(cur_ds_index)]=df_summary.loc[cur_index, cur_col]

    df_summary3d[df_summary3d == None] = 5
    for cur in range(len(algos)):
        for cur_z in range(len(datasets)):
            df_summary3d[cur, cur,cur_z]=0

    df_summary_avg=np.nanmean(df_summary3d,axis=2)
    df_summary_max = np.nanmax(df_summary3d, axis=2)
    df_summary_min = np.nanmin(df_summary3d, axis=2)



    np_similarity = df_summary_avg # df_summary.values
    labels=[constants.ALGOS_ACRONYM[a] for a in algos]

    for i, cur_ds in enumerate(datasets):
        # plot_hierarchical_algos(df_summary3d[:,:,i], labels, cur_ds)
        pd.DataFrame(data=df_summary3d[:,:,i], index=algos, columns=algos).to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "similarity_matrix_{}.tsv".format(datasets[i])),sep='\t')

    plot_hierarchical_algos(df_summary_avg, labels, "average", ax=ax)
    # plot_hierarchical_algos(df_summary_min, labels, "min", ax=ax)
    # plot_hierarchical_algos(df_summary_max, labels, "max", ax=ax)


def plot_hierarchical_algos(df_summary_avg, labels, agg_strategy, ax=None):

    distArray = ssd.squareform(df_summary_avg)
    linked = linkage(distArray, method=METHOD, metric='euclidean')
    if ax==None:
        plt.figure(figsize=(10, 7))

    dendrogram(linked,
               orientation='top',
               labels=labels,
               color_threshold=0,
               leaf_rotation=90,
               distance_sort='descending',
               show_leaf_counts=True ,ax=ax)

    ax.set_title("algos hierarchical clustering\n(aggregated by {})".format(agg_strategy))
    ax.set_ylabel("distance")
    ax.set_xlabel("algos")
    # plt.tight_layout()
    # plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "hierarchical_go_clustering_{}.png".format(agg_strategy)))

    plt.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_2.png"))






if __name__ == "__main__":

    font = {'size': 22}

    matplotlib.rc('font', **font)

    # algos = ["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox", "hotnet2", "keypathwayminer_INES_GREEDY"]
    datasets = ["TNFa_2", "HC12", "SHERA", "ROR_1", "SHEZH_1", "ERS_1", "IEM"]
    distance_file_name = os.path.join(constants.OUTPUT_GLOBAL_DIR, "cancer_type_go_distance_original_1.tsv")

    figure, axs=plt.subplots(3,1,figsize=(20,30))


    plt.figtext(0.01, 0.98, "A:", weight='bold', size=20)
    plt.figtext(0.01, 0.63, "B:", weight='bold', size=20)
    plt.figtext(0.01, 0.28, "C:", weight='bold', size=20)



    solutions_hierarchical_clustering(color_strategy=datasets, color_strategy_title='datasets', ax=axs[0], distance_file_name=distance_file_name)
    solutions_hierarchical_clustering(color_strategy=constants.ALGOS_ACRONYM.values(), color_strategy_title='algos', ax=axs[1], distance_file_name=distance_file_name)

    algos_hierarchical_clustering(algos=constants.ALGOS_ACRONYM.keys(), datasets=datasets, distance_file_name=distance_file_name, ax=axs[2])



