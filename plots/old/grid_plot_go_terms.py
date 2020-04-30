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


def plot_grid(path_format, datasets, algos):

    for dataset in datasets:
        df_grid=pd.DataFrame()
        for algo in algos:

            path=path_format.format(dataset=dataset, algo=algo)
            df_terms=pd.read_csv(path, sep='\t')
            # df_terms = df_terms[df_terms["emp_pval"]==0]
            df_terms = df_terms.iloc[:5,:]

            for cur_i, cur_row in df_terms.iterrows():
                df_grid=df_grid.append({"algo": algo, 'go_name': cur_row['GO name'], 'hg_pval': cur_row['hg_pval']}, ignore_index=True)

        go_list = np.unique(df_grid[["go_name"]])
        fig, ax = plt.subplots(figsize=(40, 40))
        ax.set_facecolor('#fffde3')
        ax.grid(color='gray')

        xs=[algos.index(a) for a in df_grid["algo"]]
        ys=[list(go_list).index(a) for a in df_grid["go_name"]]
        cs=[a for a in df_grid["hg_pval"]]
        sc = ax.scatter(xs, ys, 300, c=cs, cmap='bwr', vmin=np.percentile(cs, 10),
                        vmax=np.percentile(cs, 90))
        for x, y ,c  in zip(xs, ys, cs):
            ax.annotate( "{}".format(round(c,2)), (x, y), color='green', size=20)
        ax.legend(loc='upper left')
        ax.margins(0.03, 0.03)

        ax.set_xlabel("algos", fontdict={"size" : 35})
        plt.subplots_adjust(left=0.25, right=0.99, top=0.99, bottom=0.05)
        plt.xticks(np.arange(len(algos)), tuple(algos), rotation='vertical', size=35)
        ax.set_ylabel("algos", fontdict={"size" : 25})
        plt.yticks(np.arange(len(go_list)), tuple(go_list), size=35)
        ax_ = plt.gca()
        aspect = 20
        pad_fraction = 0.5
        divider = make_axes_locatable(ax_)
        width = axes_size.AxesY(ax_, aspect=1. / aspect)
        pad = axes_size.Fraction(pad_fraction, width)
        cax = divider.append_axes("right", size=0.3, pad=0.4)
        plt.colorbar(mappable=sc, cax=cax)
        cax.tick_params(labelsize=35)
        plt.tight_layout()
        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "grid_go_terms_{}.png".format(dataset)))
        plt.clf()


def main():
    datasets=["TNFa_2", "HC12", "SHERA", "ROR_1", "SHEZH_1", "ERS_1", "IEM"]
    algos = ["jactivemodules_greedy","jactivemodules_sa","netbox","hotnet2","bionet","keypathwayminer_INES_GREEDY"]
    path_format = "/home/hag007/Desktop/aggregate_report/oob/emp_diff_{dataset}_{algo}_passed_oob.tsv"

    plot_grid(path_format, datasets, algos)


if __name__=="__main__":
    main()