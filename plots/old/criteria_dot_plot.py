import pandas as pd

from fastsemsim.SemSim import *

import matplotlib.pyplot as plt

from rpy2.robjects import pandas2ri
pandas2ri.activate()
import constants

import matplotlib.cm as cm

from matplotlib.lines import Line2D

import seaborn as sns


from plots.old.grid_plot_criteria import plot_scatter_with_size

ENABLE_GO_GRAPH = False
IS_GO_GRAPH_ONLY=False
GO_PCA=False
TERMS_SIMILARITY_TO_NUM_OF_TERMS=True
RATIO_TO_GO_TERM = True
QVAL_TH = 0.01
SIM_TH= 0.4


algos_acronym={"jactivemodules_greedy":"jAM_greedy",
               "jactivemodules_sa": "jAM_SA",
               "netbox": "netbox",
               "keypathwayminer_INES_GREEDY": "KPM",
               "hotnet2": "hotnet2",
               "bionet": "bionet",
               "my_netbox_td": "netbox_td"
               }




def dot_grid(df_measurements, grid_type, y_label, filter_zeros=False, ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))

    ax.set_facecolor('#fffde3')
    ax.grid(color='gray')

    colormap = cm.jet

    algos=df_measurements.index.values
    sns.set_palette("husl", algos.shape[0])

    df_new = pd.DataFrame()
    for i, cur_row in df_measurements.iterrows():
        for cur_entry in cur_row:
            df_new=df_new.append({"algo": algos_acronym[i], y_label : cur_entry}, ignore_index=True)

    ax = sns.stripplot(x='algo', y=y_label, data=df_new, jitter=True, size=10, ax=ax)
    ax.set_title("{} criterion".format(y_label), fontdict={"size":22})
    colorlist = [sns.color_palette()[i] for i in
                 np.array(list(range(len(algos))))] #  / float(len(algos) - 1)]
    patches = [Line2D([0], [0], marker='o', color='gray', label=algos_acronym[a], markersize=12,
                      markerfacecolor=c, alpha=0.7) for i, a, c in zip(list(range(len(algos))), algos, colorlist)]

    patches += [Line2D([0], [0], marker='D', color='gray', markersize=0)]
    patches += [Line2D([0], [0], marker='D', color='gray', label='mean', markersize=12, markerfacecolor='red', alpha=0.7)]
    patches += [Line2D([0], [0], marker='P', color='gray', label='median', markersize=12, markerfacecolor='red', alpha=0.7)]
    i=0
    for index, row in df_measurements.iterrows():

        if filter_zeros:
            mn = row.values[row.values != 0].mean()
            mdn = np.median(row.values[row.values != 0])
        else:
            mn = row.values.mean()
            mdn = np.median(row.values)


        ax.scatter([i], [mn], marker='D', color='red', edgecolors='b', s=[200], alpha=0.7 ,zorder=4)


        ax.scatter([i], [mdn], marker='P', color='red', edgecolors='b', s=[200], alpha=0.7, zorder=4)
        i+=1

    ax.set_xlabel("algos", fontdict={"size" : 22})
    ax.set_ylabel(y_label, fontdict={"size": 22})
    ax.legend(handles=patches, loc='upper left', fontsize=22)
    # ax.set_xticks(np.arange(len(algos)),
    #            tuple([algos_acronym[x[3:]] if x.startswith("GE_") else algos_acronym[x] for x in algos]))
    ax.set_xticklabels(tuple([algos_acronym[x[3:]] if x.startswith("GE_") else algos_acronym[x] for x in algos]), rotation=45)

    plt.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "datasets_algo_dot_{}.png".format(grid_type)))



def main():

    fig, axs = plt.subplots(2,2, figsize=(20, 20))

    main_path = "/home/hag007/Desktop/aggregate_report/visual"
    df_measurements_counts=pd.read_csv(os.path.join(main_path, "empirical_terms_counts.tsv"), sep='\t', index_col=0).loc[np.sort(algos_acronym.keys()),:]
    # df_measurements_counts=df_measurements_counts.drop(labels=["SHERA"], axis=1)
    dot_grid(df_measurements_counts, "counts", "Terms Count", ax=axs[0][1])


    df_measurements_ratio = pd.read_csv(os.path.join(main_path, "true_positive_ratio.tsv"), sep='\t', index_col=0).loc[np.sort(algos_acronym.keys()),:]
    dot_grid(df_measurements_ratio, "ratio", "EHR", ax=axs[0][0])


    df_measurements_heterogeneity = pd.read_csv(os.path.join(main_path, "empirical_terms_variability.tsv"), sep='\t', index_col=0).loc[np.sort(algos_acronym.keys()),:]
    i_zeros = df_measurements_heterogeneity==0
    df_measurements_heterogeneity +=np.abs(np.min(df_measurements_heterogeneity.values))+1
    # df_measurements=1.0/df_measurements
    df_measurements_heterogeneity[i_zeros]=0
    dot_grid(df_measurements_heterogeneity, "variability", "Heterogeneity", filter_zeros=False, ax=axs[1][0])

    df_counts = pd.read_csv(os.path.join(main_path, "empirical_terms_counts.tsv"), sep='\t', index_col=0)
    plot_scatter_with_size(df_measurements_counts, df_measurements_ratio, df_measurements_heterogeneity, "size", ax=axs[1][1])

    plt.figtext(0.01,0.99, "A:", weight='bold')
    plt.figtext(0.5, 0.99, "B:", weight='bold')
    plt.figtext(0.01, 0.5, "C:", weight='bold')
    plt.figtext(0.5, 0.5, "D:", weight='bold')

    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_4.png"))

if __name__=="__main__":
    main()