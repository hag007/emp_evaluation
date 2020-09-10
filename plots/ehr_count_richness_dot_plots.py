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

    if ax is None:
        fig, ax = plt.subplots(figsize=(13, 13))

    ax.set_facecolor('#ffffff')
    ax.grid(color='gray')
    patches = [Line2D([0], [0], marker='o', markersize=12, color='gray', label=constants.ALGOS_ACRONYM[a],
                      markerfacecolor=constants.COLORDICT[a]) for i, a in zip(list(range(len(algos))), algos)]
    ax.legend(handles=patches, loc=(0.0,1.1), prop={'size': 20},ncol=2, facecolor='#ffffff')
    ax.set_xlabel("# of non-redundant GO terms", fontdict={"size":22})
    ax.set_ylabel("EHR", fontdict={"size": 22})
    ax.set_title(title, fontdict={"size":22})
    plt.subplots_adjust(left=0.25, right=0.99, top=0.99, bottom=0.05)


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
    my_order = [constants.ALGOS_ACRONYM[a] for a in algos]
    ax = sns.stripplot(x='Alg', y=y_label, data=df_new, jitter=True, size=10, ax=ax, order=my_order, palette={a: "#AAAAAA" for a in my_order}) # {a: colorlist[algo_acrs.index(a)] for a in my_order}
    if "EV terms" in title:
        df_new_copy=pd.DataFrame(df_new)
        df_new_copy.loc[:,'# EV terms'][df_new_copy.loc[:,'# EV terms']!=0]=np.nan
        ax_2 = sns.stripplot(x='Alg', y=y_label, data=df_new_copy, jitter=True, size=10, ax=ax_2, order=my_order, palette={a: "#AAAAAA" for a in my_order}) # {a: colorlist[algo_acrs.index(a)] for a in my_order}
    ax.set_xlabel("")
    ax.set_title("{}".format(title), fontdict={"size":22})

    patches = [Line2D([0], [0], marker='s', color='gray', label='mean', markersize=12, markerfacecolor='blue', alpha=0.7),  Line2D([0], [0], marker='P', color='gray', label='median', markersize=12, markerfacecolor='red', alpha=0.7)]
    i=0
    for index, row in df_measurements.iterrows():

        if filter_zeros:
            mn = np.nanmean(row.values[row.values != 0])
            mdn = np.median(row.values[row.values != 0])
        else:
            mn = np.nanmean(row.values)
            std = np.nanstd(row.values)
            mdn = np.median(row.values)


        ax.bar([list(my_order).index(constants.ALGOS_ACRONYM[index])], [mn], color='blue')
        err1 = ax.errorbar([list(my_order).index(constants.ALGOS_ACRONYM[index])-0.1], [mn], yerr=np.float32(std), lolims=True, capsize=5, ls='None', color='black', elinewidth=0.5, label='_nolegend_')
        err1[1][0].set_marker('_')
        err1[1][0].set_markersize(10)

        ax.scatter(list(my_order).index(constants.ALGOS_ACRONYM[index]), [mdn], marker='P', color='red', edgecolors='b', s=[200], alpha=0.7, zorder=5)
        i+=1

    if "EV terms" in title:
        ax.set_yscale('log', subsy=[])
        ax_2.set_ylabel("")
        ax.set_yticklabels([])

    ax.set_xlabel(ax.get_xlabel(), fontdict={"size" : 22})
    ax.set_ylabel(ax.get_ylabel(), fontdict={"size": 22})
    ax.set_title(title, fontdict={"size":22})
    ax.legend(handles=patches, loc=(0.0,1.0), fontsize=20, facecolor='#ffffff')
    ax.set_xticklabels([a for a in ax.get_xticklabels()], rotation=45)


def main(algos):

    fig_4, axs_4 = plt.subplots(2,2, figsize=(30, 30))
    main_path = os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation")


    prefix="GE"
    df_measurements_ratio=pd.read_csv(os.path.join(main_path, "ehr_matrix_{}.tsv".format(prefix)), sep='\t', index_col=0)
    dot_plot(df_measurements_ratio, "EHR", ax=axs_4[0][0], algos=algos, title="EHR, {}".format(prefix))

    df_measurements_counts = pd.read_csv(os.path.join(main_path, "count_matrix_{}.tsv".format(prefix)), sep='\t', index_col=0)
    dot_plot(df_measurements_counts, "# EV terms", ax=axs_4[1][0], algos=algos,  title="EV terms, {}".format(prefix))

    prefix = "PASCAL_SUM"
    df_measurements_ratio = pd.read_csv(os.path.join(main_path, "ehr_matrix_{}.tsv".format(prefix)), sep='\t', index_col=0)
    dot_plot(df_measurements_ratio, "EHR", ax=axs_4[0][1], algos=algos, title="EHR, {}".format(prefix))

    df_measurements_counts = pd.read_csv(os.path.join(main_path, "count_matrix_{}.tsv".format(prefix)), sep='\t', index_col=0)
    dot_plot(df_measurements_counts, "# EV terms", ax=axs_4[1][1], algos=algos, title="EV terms, {}".format(prefix))

    fig_4.tight_layout()
    fig_4.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "figure_4.png"))

