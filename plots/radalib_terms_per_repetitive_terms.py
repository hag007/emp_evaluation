import sys
sys.path.insert(0, '../')


import constants
import os
import argparse
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def main(datasets, algos, prefix, ax1=None, ax2=None, hg_th=0.05):

    filtered_go_ids_file=os.path.join(constants.GO_DIR,"filtered_go_terms.txt")
    sig_in_null_file=os.path.join(constants.GO_DIR,"sig_in_null.txt")


    df_ratio_matrix=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation", "aggregate_solutions_ratio_by_algo_{}.tsv".format(prefix)),sep='\t', index_col=0)
    df_count_matrix=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation", "aggregate_solutions_count_by_algo_{}.tsv".format(prefix)), sep='\t', index_col=0)

    df_curve_ratio=pd.DataFrame()
    df_curve_count = pd.DataFrame()

    for count in np.arange(11):
        for alg in algos:
            df_count_filtered=df_count_matrix.loc[(df_count_matrix.loc[:, alg] == count).values, "is_sig_in_radalib"]
            df_curve_count.loc[alg, count] = df_count_filtered.sum() / float(df_count_filtered.shape[0])
            df_ratio_filtered = df_ratio_matrix.loc[(df_count_matrix.loc[:, alg] == count).values, alg]
            df_curve_ratio.loc[alg, count] = df_ratio_filtered.mean()

    ax1.set_facecolor('#ffffff')
    ax1.grid(color='gray')
    ax2.set_facecolor('#ffffff')
    ax2.grid(color='gray')

    for alg in algos:
        ax1.plot(np.arange(11),df_curve_count.loc[alg,:],label=constants.ALGOS_ACRONYM[alg], c=constants.COLORDICT[alg])
        ax2.plot(np.arange(11), df_curve_ratio.loc[alg, :], label=constants.ALGOS_ACRONYM[alg],
                c=constants.COLORDICT[alg])


    patches = [Line2D([0], [0], marker='o', color='gray', label=constants.ALGOS_ACRONYM[a], markersize=12,
                      markerfacecolor=constants.COLORDICT[a], alpha=0.7) for a in algos]

    ax1.set_xlabel("# of repetitions")
    ax1.set_ylabel("N-terms proportion")
    ax2.set_xlabel("# of repetitions")
    ax2.set_ylabel("non-EV ratio")


    ax1.legend(handles=patches, loc=(0,1.1), ncol=3, fontsize=22, facecolor='#ffffff')
    ax2.legend(handles=patches, loc=(0, 1.1), ncol=3, fontsize=22, facecolor='#ffffff')



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="TNFa_2,HC12,ROR_1,SHERA,SHEZH_1,ERS_1,IEM,APO,CBX,IFT")
    parser.add_argument('--prefix', dest='prefix', default="GE") # PASCAL_SUM   GE
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy,jactivemodules_sa,bionet,netbox,keypathwayminer_INES_GREEDY,DOMINO,DOMINO2,hotnet2")
    args = parser.parse_args()
    datasets = args.datasets.split(",")
    algos = args.algos.split(",")
    prefix = args.prefix

    fig,axs=plt.subplots(2,2, figsize=(22,22))

    datasets=['tnfa', 'hc', 'ror', 'cbx', 'shera' , 'shezh' , 'ift', 'iem' , "ers", "apo"]
    algos=["DOMINO2", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    prefix="GE"
    main(datasets,algos,prefix, axs[0][0], axs[1][0])

    datasets=['brca', 'crh', 'scz', 'tri', 't2d', 'af', 'cad', 'amd', 'hgt' , "bmd"]
    algos=["DOMINO2", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    prefix="PASCAL_SUM"
    main(datasets,algos,prefix, axs[0][1], axs[1][1])

    # fig.text(0.01,0.97, "A:", weight='bold',fontsize=22)
    # fig.text(0.5, 0.97, "B:", weight='bold',fontsize=22)
    fig.tight_layout()
    fig.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "radalib_curve.png"))
