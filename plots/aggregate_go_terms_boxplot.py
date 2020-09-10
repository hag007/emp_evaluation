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



    df_measurements=df_measurements.loc[:, algos].T
    sig_in_null_file=os.path.join(constants.GO_DIR,"sig_in_null.txt")
    sig_in_null=open(sig_in_null_file,'r').read().split()
    df_sig_in_null=df_measurements.loc[:,sig_in_null]
    df_non_sig_in_null=df_measurements.loc[:,~df_measurements.columns.isin(sig_in_null)]

    # patches += [Line2D([0], [0], marker='D', color='gray', markersize=0)]
    patches = [Line2D([0], [0], marker='s', color='gray', label='sig. terms in Radalib', markersize=12, markerfacecolor='green', alpha=0.7),  Line2D([0], [0], marker='s', color='gray', label='non sig. terms in Radalib', markersize=12, markerfacecolor='orange', alpha=0.7)]
    i=0
    my_order = [constants.ALGOS_ACRONYM[a] for a in algos] # df_new.groupby(by=["algo"])[y_label].mean().sort_values().index

    for df_sig_in_null_elem, df_non_sig_in_null_elem in zip(df_sig_in_null.iterrows(), df_non_sig_in_null.iterrows()):

        index_sig, row_sig = df_sig_in_null_elem
        index_non_sig, row_non_sig = df_non_sig_in_null_elem
        row_sig=row_sig.dropna()
        row_non_sig=row_non_sig.dropna()

        mn_sig = np.nanmean(row_sig.values) if not np.isnan(np.nanmean(row_sig.values)) else 0
        std_sig = np.nanstd(row_sig.values)

        mn_non_sig =  np.nanmean(row_non_sig.values) if not np.isnan(np.nanmean(row_non_sig.values)) else 0
        std_non_sig = np.nanstd(row_non_sig.values)

        x=list(my_order).index(constants.ALGOS_ACRONYM[index_sig])
        max_y=max(mn_sig, mn_non_sig)+0.1

        ax.bar(np.array([x])-0.05, [mn_sig], color='green', width=0.2)
        # err1 = ax.errorbar([list(my_order).index(constants.ALGOS_ACRONYM[index_sig]+"_{}".format(df_measurements.loc[index_sig].dropna().shape[0]))], [mn_sig], yerr=np.float32(std_sig), lolims=True, capsize=5, ls='None', color='black', elinewidth=0.5, label='_nolegend_')
        # err1[1][0].set_marker('_')
        # err1[1][0].set_markersize(10)

        ax.bar(np.array([x])+0.25, [mn_non_sig], color='orange', width=0.2)
        # err1 = ax.errorbar(np.array([list(my_order).index(constants.ALGOS_ACRONYM[index_non_sig]+"_{}".format(df_measurements.loc[index_non_sig].dropna().shape[0]))])+0.2, [mn_non_sig], yerr=np.float32(std_non_sig), lolims=True, capsize=5, ls='None', color='black', elinewidth=0.5, label='_nolegend_')
        # err1[1][0].set_marker('_')
        # err1[1][0].set_markersize(10)

        ax.plot([x, x, x+0.2, x+0.2], [max_y, max_y+0.05, max_y+0.05, max_y], lw=1.5, c='black')
        ax.text(x-0.15, mn_sig, "n={}".format(len(row_sig.values)), ha='center', va='bottom')
        ax.text(x+0.35, mn_non_sig, "n={}".format(len(row_non_sig.values)), ha='center', va='bottom')
        ax.text((x+x+0.2)*.5, max_y+0.05, "{:.2e}".format(mannwhitneyu(row_sig.values, row_non_sig.values)[1]), ha='center', va='bottom')


    ax.set_xlabel(ax.get_xlabel(), fontdict={"size" : 22})
    ax.set_ylabel(ax.get_ylabel(), fontdict={"size": 22})
    ax.set_title(title, fontdict={"size":22})
    ax.legend(handles=patches, loc=(0.0,1.0), fontsize=20, facecolor='#ffffff')
    ax.set_xticklabels([""]+my_order, rotation=45)

    # plt.tight_layout()
    # plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "datasets_algo_dot_{}.png".format(title)))



def dot_plot(df_measurements, y_label, filter_zeros=False, ax=None, algos=None, title=""):
    plt.minorticks_off()
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))

    ax.set_facecolor('#ffffff')
    ax.grid(color='gray')

    if "EV terms" in title:
        ax_2=ax.twinx()
        ax_2.set_yticks([])
        ax_2.set_ylim([-0.005,1])
        ax_2.grid(False)


    df_measurements=df_measurements.loc[:, algos].T
    df_new = pd.DataFrame()
    for i, cur_row in df_measurements.iterrows():
        for cur_entry in cur_row:
            df_new=df_new.append({"Alg": constants.ALGOS_ACRONYM[i], y_label : cur_entry}, ignore_index=True)


    df_new=df_new.dropna(axis=0)
    my_order = [constants.ALGOS_ACRONYM[a] for a in algos] # df_new.groupby(by=["algo"])[y_label].mean().sort_values().index
    ax = sns.violinplot(x='Alg', y=y_label, data=df_new, ax=ax, order=my_order, palette={constants.ALGOS_ACRONYM[a]: constants.COLORDICT[a] for a in algos}, jitter=True, size=10) # {a: colorlist[algo_acrs.index(a)] for a in my_order}  # jitter=True, ,  size=10
    if "EV terms" in title:
        df_new_copy=pd.DataFrame(df_new)
        df_new_copy.loc[:,'# EV terms'][df_new_copy.loc[:,'# EV terms']!=0]=np.nan
        ax_2 = sns.violinplot(x='Alg', y=y_label, data=df_new_copy, ax=ax_2, order=my_order, palette={constants.ALGOS_ACRONYM[a]: constants.COLORDICT[a] for a in algos}, jitter=True, size=10) # {a: colorlist[algo_acrs.index(a)] for a in my_order} #
    ax.set_xlabel("")
    ax.set_title("{}".format(title), fontdict={"size":22})

    for i, a in enumerate(algos):
        ax.text(i, 1.1, "n={}".format(df_measurements.loc[a].dropna().shape[0]), ha='center', va='bottom')

    # patches += [Line2D([0], [0], marker='D', color='gray', markersize=0)]
    patches = [Line2D([0], [0], marker='s', color='gray', label='mean', markersize=12, markerfacecolor='blue', alpha=0.7),  Line2D([0], [0], marker='P', color='gray', label='median', markersize=12, markerfacecolor='red', alpha=0.7)]
    i=0
    for index, row in df_measurements.iterrows():

        # if constants.ALGOS_ACRONYM[index] not in my_order:
        #     continue

        if filter_zeros:
            mn = np.nanmean(row.values[row.values != 0])
            mdn = np.median(row.values[row.values != 0])
        else:
            mn = np.nanmean(row.values)
            std = np.nanstd(row.values)
            mdn = np.median(row.values)


        # ax.bar([list(my_order).index(constants.ALGOS_ACRONYM[index])], [mn], color='blue')#, marker='D', color='red', edgecolors='b', s=[250], alpha=0.7 ,zorder=4)
        # err1 = ax.errorbar([list(my_order).index(constants.ALGOS_ACRONYM[index])-0.1], [mn], yerr=np.float32(std), lolims=True, capsize=5, ls='None', color='black', elinewidth=0.5, label='_nolegend_')
        # err1[1][0].set_marker('_')
        # err1[1][0].set_markersize(10)

        # ax.scatter(list(my_order).index(constants.ALGOS_ACRONYM[index]), [mdn], marker='P', color='red', edgecolors='b', s=[200], alpha=0.7, zorder=5)
        i+=1


    # # else:
    # #     ax.set_yticklabels([a.get_text() if float(np.abs(a.get_position()[1])) <= 1 else "" for a in ax.get_yticklabels()])

    ax.set_ylim([-0.1,1.2])
    ax.set_xlabel(ax.get_xlabel(), fontdict={"size" : 22})
    ax.set_ylabel(ax.get_ylabel(), fontdict={"size": 22})
    ax.set_title(title, fontdict={"size":22})
    # ax.legend(handles=patches, loc=(0.0,1.0), fontsize=20, facecolor='#ffffff')
    ax.set_xticklabels([a for a in ax.get_xticklabels()], rotation=45)

    # plt.tight_layout()
    # plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "datasets_algo_dot_{}.png".format(title)))



def main():

    fig_4, axs_4 = plt.subplots(2,2, figsize=(40, 20))



    prefix="GE"
    algos=["DOMINO2" , "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    main_path=os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation")
    df_measurements_ratio=pd.read_csv(os.path.join(main_path, "aggregate_solutions_ratio_by_algo_{}.tsv".format(prefix)), sep='\t', index_col=0) # .loc[np.sort(constants.ALGOS_ACRONYM.keys()),:]
    dot_plot(df_measurements_ratio, "non-EV ratio", ax=axs_4[0][0], algos=algos, title="non-EV EHR ratio, {}".format(prefix))

    prefix = "PASCAL_SUM"
    algos=["DOMINO2", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    df_measurements_ratio=pd.read_csv(os.path.join(main_path, "aggregate_solutions_ratio_by_algo_{}.tsv".format(prefix)), sep='\t', index_col=0) # .loc[np.sort(constants.ALGOS_ACRONYM.keys()),:]
    dot_plot(df_measurements_ratio, "non-EV ratio", ax=axs_4[0][1], algos=algos, title="non-EV EHR ratio, {}".format("GWAS"))

    prefix="GE"
    algos=["DOMINO2" , "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    main_path=os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation")
    df_measurements_ratio=pd.read_csv(os.path.join(main_path, "aggregate_solutions_ratio_by_algo_{}.tsv".format(prefix)), sep='\t', index_col=0) # .loc[np.sort(constants.ALGOS_ACRONYM.keys()),:]
    bar_plot(df_measurements_ratio, "non-EV ratio", ax=axs_4[1][0], algos=algos, title="non-EV EHR ratio, {}".format(prefix))

    prefix = "PASCAL_SUM"
    algos=["DOMINO2", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    df_measurements_ratio=pd.read_csv(os.path.join(main_path, "aggregate_solutions_ratio_by_algo_{}.tsv".format(prefix)), sep='\t', index_col=0) # .loc[np.sort(constants.ALGOS_ACRONYM.keys()),:]
    bar_plot(df_measurements_ratio, "non-EV ratio", ax=axs_4[1][1], algos=algos, title="non-EV EHR ratio, {}".format("GWAS"))


    # fig_4.text(0.01,0.98, "A", weight='bold',fontsize=22)
    # fig_4.text(0.51, 0.98, "B", weight='bold',fontsize=22)
    # fig_4.text(0.01, 0.49, "C", weight='bold',fontsize=22)
    # fig_4.text(0.51, 0.49, "D", weight='bold',fontsize=22)
    fig_4.tight_layout()
    fig_4.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "figure_20.png"))

    # fig_5.text(0.01, 0.97, "A", weight='bold', fontsize=22)
    # fig_5.text(0.5, 0.97, "B", weight='bold', fontsize=22)
    # fig_5.tight_layout()
    # fig_5.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "figure_5.png"))

if __name__=="__main__":
    main()