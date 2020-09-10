
import matplotlib.pyplot as plt
import pandas as pd
from fastsemsim.SemSim import *

import constants

from matplotlib.lines import Line2D

import seaborn as sns
from scipy.stats import mannwhitneyu


ENABLE_GO_GRAPH = False
IS_GO_GRAPH_ONLY=False
GO_PCA=False
TERMS_SIMILARITY_TO_NUM_OF_TERMS=True
RATIO_TO_GO_TERM = True
QVAL_TH = 0.01
SIM_TH= 0.4

def bar_plot(df_measurements, y_label, filter_zeros=False, ax=None, algos=None, title=""):
    plt.minorticks_off()

    ax[0].set_facecolor('#ffffff')
    ax[0].grid(color='gray')
    ax[0].xaxis.grid(False)
    ax0=ax[0]

    ax[1].set_facecolor('#ffffff')
    ax[1].grid(color='gray')
    ax[1].xaxis.grid(False)
    ax1=ax[1]

    df_measurements=df_measurements.loc[:, algos].T
    sig_in_null_file=os.path.join(constants.GO_DIR,"sig_in_null.txt")
    sig_in_null=open(sig_in_null_file,'r').read().split()
    df_n_terms=df_measurements.loc[:,sig_in_null]
    df_non_n_terms=df_measurements.loc[:,~df_measurements.columns.isin(sig_in_null)]

    # patches += [Line2D([0], [0], marker='D', color='gray', markersize=0)]
    patches = [Line2D([0], [0], marker='s', color='gray', label='sig. terms in Radalib', markersize=12, markerfacecolor='palegreen', alpha=0.7),  Line2D([0], [0], marker='s', color='gray', label='non sig. terms in Radalib', markersize=12, markerfacecolor='orange', alpha=0.7)]
    i=0
    my_order = [constants.ALGOS_ACRONYM[a] for a in algos] # df_new.groupby(by=["algo"])[y_label].mean().sort_values().index

    for df_n_terms_elem, df_non_sig_in_null_elem in zip(df_n_terms.iterrows(), df_non_n_terms.iterrows()):

        index_sig, n_terms = df_n_terms_elem
        index_non_sig, row_non_n_terms = df_non_sig_in_null_elem
        n_terms=n_terms.dropna()
        row_non_n_terms=row_non_n_terms.dropna()

        mn_sig = np.nanmean(n_terms.values) if not np.isnan(np.nanmean(n_terms.values)) else 0
        std_sig = np.nanstd(n_terms.values)

        mn_non_sig =  np.nanmean(row_non_n_terms.values) if not np.isnan(np.nanmean(row_non_n_terms.values)) else 0
        std_non_sig = np.nanstd(row_non_n_terms.values)

        x=list(my_order).index(constants.ALGOS_ACRONYM[index_sig])
        max_y=max(mn_sig, mn_non_sig)+0.01

        ax0.bar(np.array([x])-0.05, [mn_sig], color='palegreen', width=0.2)
        ax0.bar(np.array([x])+0.25, [mn_non_sig], color='orange', width=0.2)

        ax0.plot([x, x, x+0.2, x+0.2], [max_y, max_y+0.05, max_y+0.05, max_y], lw=1.5, c='black')
        # ax.text(x-0.15, 0, "n={}".format(len(n_terms.values)), ha='center', va='bottom')
        # ax.text(x+0.35, 0, "n={}".format(len(row_non_n_terms.values)), ha='center', va='bottom')
        try:
            mwu=mannwhitneyu(n_terms.values, row_non_n_terms.values, alternative='greater')[1]
        except Exception:
            mwu=np.nan

        ax0.text((x+x+0.2)*.5, max_y+0.05, "{:.2e}".format(mwu), ha='center', va='bottom')

        ax1.bar(np.array([x])-0.05, len(n_terms.values), color='palegreen', width=0.2)
        ax1.bar(np.array([x])+0.25, len(row_non_n_terms.values), color='orange', width=0.2)


        ax0.set_xlabel(ax0.get_xlabel(), fontdict={"size" : 22})
        ax0.set_ylabel("non-EV ratio", fontdict={"size": 22})
        ax0.set_title("non-EV ratio, " + title, fontdict={"size":22})
        ax0.legend(handles=patches, loc=(0.0,1.0), fontsize=20, facecolor='#ffffff')
        ax0.set_xticklabels([""]+my_order, rotation=45)

        ax1.set_xlabel(ax1.get_xlabel(), fontdict={"size" : 22})
        ax1.set_ylabel("# of terms", fontdict={"size": 22})
        ax1.set_title("# of terms, " + title, fontdict={"size":22})
        ax1.legend(handles=patches, loc=(0.0,1.0), fontsize=20, facecolor='#ffffff')
        ax1.set_xticklabels([""]+my_order, rotation=45)

    # plt.tight_layout()
    # plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "datasets_algo_dot_{}.png".format(title)))



def dot_plot(df_measurements, y_label, filter_zeros=False, ax=None, algos=None, title=""):
    plt.minorticks_off()
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))

    ax.set_facecolor('#ffffff')
    ax.grid(color='gray')

    df_measurements=df_measurements.loc[:, algos].T
    df_new = pd.DataFrame()
    for i, cur_row in df_measurements.iterrows():
        for cur_entry in cur_row:
            df_new=df_new.append({"Alg": constants.ALGOS_ACRONYM[i], y_label : cur_entry}, ignore_index=True)


    df_new=df_new.dropna(axis=0)
    my_order = [constants.ALGOS_ACRONYM[a] for a in algos] # df_new.groupby(by=["algo"])[y_label].mean().sort_values().index
    ax = sns.violinplot(x='Alg', y=y_label, data=df_new, ax=ax, order=my_order, palette={constants.ALGOS_ACRONYM[a]: constants.COLORDICT[a] for a in algos}, jitter=True, size=10) # {a: colorlist[algo_acrs.index(a)] for a in my_order}  # jitter=True, ,  size=10
    ax.set_xlabel("")
    ax.set_title("{}".format(title), fontdict={"size":22})

    for i, a in enumerate(algos):
        ax.text(i, 1.1, "n={}".format(df_measurements.loc[a].dropna().shape[0]), ha='center', va='bottom')

    ax.set_ylim([-0.1,1.2])
    ax.set_xlabel(ax.get_xlabel(), fontdict={"size" : 22})
    ax.set_ylabel(ax.get_ylabel(), fontdict={"size": 22})
    ax.set_title(title, fontdict={"size":22})
    # ax.legend(handles=patches, loc=(0.0,1.0), fontsize=20, facecolor='#ffffff')
    ax.set_xticklabels([a for a in ax.get_xticklabels()], rotation=45)


def main():


    limits=np.arange(4,11)
    for limit in limits:

        fig_0, axs_0 = plt.subplots(1, 2, figsize=(30, 10))
        fig_1, axs_1 = plt.subplots(2, 2, figsize=(40, 20))

        print("cur limit: {}".format(limit))

        prefix="GE"
        algos=["jactivemodules_sa"] # ["DOMINO2" , "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
        main_path=os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation")
        df_measurements_ratio=pd.read_csv(os.path.join(main_path, "aggregate_solutions_ratio_by_algo_{}_{}.tsv".format(limit, prefix)), sep='\t', index_col=0) # .loc[np.sort(constants.ALGOS_ACRONYM.keys()),:]
        dot_plot(df_measurements_ratio, "non-EV ratio", ax=axs_0[0], algos=algos, title="non-EV ratio, {}".format(prefix))

        # prefix = "PASCAL_SUM"
        # algos=["jactivemodules_sa"] # ["DOMINO2", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
        # df_measurements_ratio=pd.read_csv(os.path.join(main_path, "aggregate_solutions_ratio_by_algo_{}_{}.tsv".format(limit, prefix)), sep='\t', index_col=0) # .loc[np.sort(constants.ALGOS_ACRONYM.keys()),:]
        # dot_plot(df_measurements_ratio, "non-EV ratio", ax=axs_0[1], algos=algos, title="non-EV ratio, {}".format("GWAS"))

        prefix="GE"
        algos=["jactivemodules_sa"] # ["DOMINO2" , "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
        main_path=os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation")
        df_measurements_ratio=pd.read_csv(os.path.join(main_path, "aggregate_solutions_ratio_by_algo_{}_{}.tsv".format(limit, prefix)), sep='\t', index_col=0) # .loc[np.sort(constants.ALGOS_ACRONYM.keys()),:]
        bar_plot(df_measurements_ratio, "non-EV ratio", ax=axs_1[:,0], algos=algos, title=", {}".format(prefix))

        # prefix = "PASCAL_SUM"
        # algos=["jactivemodules_sa"]  #["DOMINO2", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
        # df_measurements_ratio=pd.read_csv(os.path.join(main_path, "aggregate_solutions_ratio_by_algo_{}_{}.tsv".format(limit, prefix)), sep='\t', index_col=0) # .loc[np.sort(constants.ALGOS_ACRONYM.keys()),:]
        # bar_plot(df_measurements_ratio, "non-EV ratio", ax=axs_1[:,1], algos=algos, title=", {}".format("GWAS"))


        # fig_4.text(0.01,0.98, "A", weight='bold',fontsize=22)
        # fig_4.text(0.51, 0.98, "B", weight='bold',fontsize=22)
        # fig_4.text(0.01, 0.49, "C", weight='bold',fontsize=22)
        # fig_4.text(0.51, 0.49, "D", weight='bold',fontsize=22)
        fig_0.tight_layout()
        fig_0.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "figure_20_{}.png".format(limit)))

        fig_1.tight_layout()
        fig_1.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "figure_21_{}.png".format(limit)))

        # fig_5.text(0.01, 0.97, "A", weight='bold', fontsize=22)
        # fig_5.text(0.5, 0.97, "B", weight='bold', fontsize=22)
        # fig_5.tight_layout()
        # fig_5.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "figure_5.png"))

        plt.clf()

if __name__=="__main__":
    main()