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



ENABLE_GO_GRAPH = False
IS_GO_GRAPH_ONLY=False
GO_PCA=False
TERMS_SIMILARITY_TO_NUM_OF_TERMS=True
RATIO_TO_GO_TERM = True
QVAL_TH = 0.01
SIM_TH= 0.4


def scatter_plot(df_size, df_ratio, ax=None, title="", algos=constants.ALGOS_ACRONYM.keys()):
    df_size = df_size.loc[set(algos).intersection(df_size.index).intersection(df_ratio.index)]
    df_ratio = df_ratio.loc[set(algos).intersection(df_size.index).intersection(df_ratio.index)]
    df_size[df_size.isna()]=0
    df_ratio[df_ratio.isna()]=0
    df_size=df_size.loc[np.sort(df_size.index), np.sort(df_size.columns)]
    df_ratio=df_ratio.loc[np.sort(df_ratio.index), np.sort(df_ratio.columns)]

    if ax is None:
        fig, ax = plt.subplots(figsize=(13, 13))

    ax.set_facecolor('#ffffff')
    ax.grid(color='gray')

    i_x = np.arange(len(df_size.columns))
    i_y = np.arange(len(df_size.index))
    size = np.array([df_size.iloc[b, a] for a in i_x for b in i_y])
    x = np.array([df_ratio.iloc[b, a] for a in i_x for b in i_y])
    y = np.array([df_size.iloc[b, a] for a in i_x for b in i_y])
    labels=np.array([constants.COLORDICT[b] for a in df_size.columns for b in df_size.index])
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


    patches = [Line2D([0], [0], marker='o', markersize=12, color='gray', label=constants.ALGOS_ACRONYM[a],
                      markerfacecolor=constants.COLORDICT[a]) for i, a in zip(list(range(len(algos))), algos)]
    ax.legend(handles=patches, loc=(0.0,1.1), prop={'size': 20},ncol=2, facecolor='#ffffff')
    # ax.margins(0.03, 0.03)
    ax.set_xlabel("# of non-redundant GO terms", fontdict={"size":22})
    ax.set_ylabel("EHR", fontdict={"size": 22})
    ax.set_title(title, fontdict={"size":22})
    plt.subplots_adjust(left=0.25, right=0.99, top=0.99, bottom=0.05)


# def plot_agg_scatter(df_size, df_ratio, ax=None, title=""):
#     df_size = df_size.loc[set(constants.ALGOS).intersection(df_size.index).intersection(df_ratio.index)].dropna()
#     df_ratio = df_ratio.loc[set(constants.ALGOS).intersection(df_size.index).intersection(df_ratio.index)].dropna()
#
#     df_size=df_size.loc[np.sort(df_size.index), np.sort(df_size.columns)]
#     df_ratio=df_ratio.loc[np.sort(df_ratio.index), np.sort(df_ratio.columns)]
#     algos = list(df_size.index)
#     sns.set_palette("husl", len(algos))
#
#     if ax is None:
#         fig, ax = plt.subplots(figsize=(13, 13))
#
#     ax.set_facecolor('#ffffff')
#     ax.grid(color='gray')
#
#     datasets = df_size.columns
#     algos = df_size.index
#     for cur_algo in algos:
#         y=df_size.loc[cur_algo,:].mean()
#         x=df_ratio.loc[cur_algo,:].median()
#         sc = ax.scatter(x, y, s=200, c=constants.COLORDICT[cur_algo])
#
#     patches = [Line2D([0], [0], marker='o', markersize=12, color='gray', label=constants.ALGOS_ACRONYM[a],
#                       markerfacecolor=constants.COLORDICT[a]) for i, a in zip(list(range(len(algos))), algos)]
#     ax.legend(handles=patches, loc=(0.0,1.1), prop={'size': 20},ncol=2, facecolor='#ffffff')
#     # ax.margins(0.03, 0.03)
#     ax.set_xlabel("non-redundant GO terms", fontdict={"size":22})
#     ax.set_ylabel("EHR", fontdict={"size": 22})
#     ax.set_title(title, fontdict={"size":22})
#     plt.subplots_adjust(left=0.25, right=0.99, top=0.99, bottom=0.05)




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


    df_measurements=df_measurements.loc[algos].dropna()
    df_new = pd.DataFrame()
    for i, cur_row in df_measurements.iterrows():
        for cur_entry in cur_row:
            df_new=df_new.append({"Alg": constants.ALGOS_ACRONYM[i], y_label : cur_entry}, ignore_index=True)


    df_new=df_new.dropna(axis=0)
    my_order = [constants.ALGOS_ACRONYM[a] for a in algos] # df_new.groupby(by=["algo"])[y_label].mean().sort_values().index
    ax = sns.stripplot(x='Alg', y=y_label, data=df_new, jitter=True, size=10, ax=ax, order=my_order, palette={a: "#AAAAAA" for a in my_order}) # {a: colorlist[algo_acrs.index(a)] for a in my_order}
    if "EV terms" in title:
        df_new_copy=pd.DataFrame(df_new)
        df_new_copy.loc[:,'# EV terms'][df_new_copy.loc[:,'# EV terms']!=0]=np.nan
        ax_2 = sns.stripplot(x='Alg', y=y_label, data=df_new_copy, jitter=True, size=10, ax=ax_2, order=my_order, palette={a: "#AAAAAA" for a in my_order}) # {a: colorlist[algo_acrs.index(a)] for a in my_order}
    ax.set_xlabel("")
    ax.set_title("{}".format(title), fontdict={"size":22})

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


        ax.bar([list(my_order).index(constants.ALGOS_ACRONYM[index])], [mn], color='blue')#, marker='D', color='red', edgecolors='b', s=[250], alpha=0.7 ,zorder=4)
        err1 = ax.errorbar([list(my_order).index(constants.ALGOS_ACRONYM[index])-0.1], [mn], yerr=np.float32(std), lolims=True, capsize=5, ls='None', color='black', elinewidth=0.5, label='_nolegend_')
        err1[1][0].set_marker('_')
        err1[1][0].set_markersize(10)

        ax.scatter(list(my_order).index(constants.ALGOS_ACRONYM[index]), [mdn], marker='P', color='red', edgecolors='b', s=[200], alpha=0.7, zorder=5)
        i+=1

    if "EV terms" in title:
        ax.set_yscale('log', subsy=[])
        ax_2.set_ylabel("")
        ax.set_yticklabels([])

    # # else:
    # #     ax.set_yticklabels([a.get_text() if float(np.abs(a.get_position()[1])) <= 1 else "" for a in ax.get_yticklabels()])

    ax.set_xlabel(ax.get_xlabel(), fontdict={"size" : 22})
    ax.set_ylabel(ax.get_ylabel(), fontdict={"size": 22})
    ax.set_title(title, fontdict={"size":22})
    ax.legend(handles=patches, loc=(0.0,1.0), fontsize=20, facecolor='#ffffff')
    ax.set_xticklabels([a for a in ax.get_xticklabels()], rotation=45)

    # plt.tight_layout()
    # plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "datasets_algo_dot_{}.png".format(title)))



def main():

    fig_4, axs_4 = plt.subplots(2,2, figsize=(20, 20))
    fig_5, axs_5 = plt.subplots(1,2, figsize=(25, 12))

    prefix="GE"
    algos=["DOMINO4", "netbox2_string"] # ["jactivemodules_greedy", "netbox", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY" ,"hotnet2", "DOMINO", "DOMINO2"]
    main_path =os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation")
    df_measurements_ratio=pd.read_csv(os.path.join(main_path, "ehr_matrix_{}.tsv".format(prefix)), sep='\t', index_col=0) # .loc[np.sort(constants.ALGOS_ACRONYM.keys()),:]

    dot_plot(df_measurements_ratio, "EHR", ax=axs_4[0][0], algos=algos, title="EHR, {}".format(prefix))

    df_measurements_counts = pd.read_csv(os.path.join(main_path, "count_matrix_{}.tsv".format(prefix)), sep='\t', index_col=0) # .loc[np.sort(constants.ALGOS_ACRONYM.keys()),:] # ratio_matrix.tsv
    dot_plot(df_measurements_counts, "# EV terms", ax=axs_4[1][0], algos=algos,  title="EV terms, {}".format(prefix))

    df_measurements_richness=pd.read_csv(os.path.join(main_path, "richness_matrix_{}_2.0.tsv".format(prefix)), sep='\t', index_col=0) # .loc[np.sort(constants.ALGOS_ACRONYM.keys()),:]

    scatter_plot(df_measurements_ratio, df_measurements_richness, ax=axs_5[0], title="EHR-Richness, {}".format(prefix), algos=algos)

    prefix = "PASCAL_SUM"
    algos=["DOMINO4", "netbox2_string"] # ["jactivemodules_greedy", "netbox", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY","hotnet2", "DOMINO", "DOMINO2"]
    df_measurements_ratio=pd.read_csv(os.path.join(main_path, "ehr_matrix_{}.tsv".format(prefix)), sep='\t', index_col=0) # .loc[np.sort(constants.ALGOS_ACRONYM.keys()),:]
    dot_plot(df_measurements_ratio, "EHR", ax=axs_4[0][1], algos=algos, title="EHR, {}".format("GWAS"))

    df_measurements_counts = pd.read_csv(os.path.join(main_path, "count_matrix_{}.tsv".format(prefix)), sep='\t', index_col=0) # .loc[np.sort(constants.ALGOS_ACRONYM.keys()),:] # ratio_matrix.tsv
    dot_plot(df_measurements_counts, "# EV terms", ax=axs_4[1][1], algos=algos, title="EV terms, {}".format("GWAS"))

    df_measurements_richness=pd.read_csv(os.path.join(main_path, "richness_matrix_{}_2.0.tsv".format(prefix)), sep='\t', index_col=0) # .loc[np.sort(constants.ALGOS_ACRONYM.keys()),:]

    scatter_plot(df_measurements_ratio, df_measurements_richness, ax=axs_5[1], title="EHR-Richness, {}".format(prefix), algos=algos)


    # fig_4.text(0.01,0.98, "A", weight='bold',fontsize=22)
    # fig_4.text(0.51, 0.98, "B", weight='bold',fontsize=22)
    # fig_4.text(0.01, 0.49, "C", weight='bold',fontsize=22)
    # fig_4.text(0.51, 0.49, "D", weight='bold',fontsize=22)
    fig_4.tight_layout()
    fig_4.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "figure_4.png"))

    # fig_5.text(0.01, 0.97, "A", weight='bold', fontsize=22)
    # fig_5.text(0.5, 0.97, "B", weight='bold', fontsize=22)
    fig_5.tight_layout()
    fig_5.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "figure_5.png"))

if __name__=="__main__":
    main()
