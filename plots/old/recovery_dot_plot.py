import pandas as pd

from fastsemsim.SemSim import *

import matplotlib.pyplot as plt

from rpy2.robjects import pandas2ri
pandas2ri.activate()
import constants

import matplotlib.cm as cm

from matplotlib.lines import Line2D

import seaborn as sns

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
               "my_netbox_td": "my_netbox_td",
               "bionet": "bionet",
               "hotnet2": "hotnet2",
               "keypathwayminer_INES_GREEDY": "KPM",
               "dcem": "dcem"
}

# "hotnet2": "hotnet2",
# "keypathwayminer_INES_GREEDY": "KPM",

# "my_netbox_td": "netbox_td", "keypathwayminer_INES_GREEDY": "KPM",  "hotnet2": "hotnet2",




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
    ax.set_title("{} criterion".format(y_label))
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
            mn = row.values[row.values != 0 & ~np.isnan(row.values)].mean()
            mdn = np.median(row.values[row.values != 0 & ~np.isnan(row.values)])
        else:
            mn = row.values[~np.isnan(row.values)].mean()
            mdn = np.median(row.values[~np.isnan(row.values)])

        ax.scatter([i], [mn], marker='D', color='red', edgecolors='b', s=[200], alpha=0.7 ,zorder=4)


        ax.scatter([i], [mdn], marker='P', color='red', edgecolors='b', s=[200], alpha=0.7, zorder=4)
        i+=1

    ax.set_xlabel("algos", fontdict={"size" : 12})
    ax.set_ylabel(y_label, fontdict={"size": 12})
    ax.legend(handles=patches, loc='upper left', fontsize=12)
    # ax.set_xticks(np.arange(len(algos)),
    #            tuple([algos_acronym[x[3:]] if x.startswith("GE_") else algos_acronym[x] for x in algos]))
    ax.set_xticklabels(tuple([algos_acronym[x[3:]] if x.startswith("GE_") else algos_acronym[x] for x in algos]), rotation=45)

    plt.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "datasets_algo_dot_{}.png".format(grid_type)))



def main():

    fig, axs = plt.subplots(2,2, figsize=(20, 20))

    suffix='PASCAL_SUM_100_0.1'

    if suffix.startswith("GE"):
        omic_type=''
    elif suffix.startswith("PASCAL"):
        omic_type='_gwas'

    main_path = os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr/MAX")
    df_measurements=pd.read_csv(os.path.join(main_path, "recovery_results_{}_matrix_p.tsv".format(suffix)), sep='\t', index_col=0).loc[np.sort(algos_acronym.keys()),:]

    df_zeros = pd.read_csv(os.path.join('/home/hag007/Desktop/aggregate{}_report/venn/'.format(omic_type), "count_matrix.tsv"), sep='\t', index_col=0).loc[np.sort(algos_acronym.keys()), df_measurements.columns.values]
    df_measurements[np.logical_or(df_zeros==0, np.isnan(df_zeros)).values]=np.nan

    dot_grid(df_measurements, "precision", "Precision", ax=axs[0][0])

    df_measurements=pd.read_csv(os.path.join(main_path, "recovery_results_{}_matrix_r.tsv".format(suffix)), sep='\t', index_col=0).loc[np.sort(algos_acronym.keys()),:]

    df_zeros = pd.read_csv(os.path.join('/home/hag007/Desktop/aggregate{}_report/venn/'.format(omic_type), "count_matrix.tsv"),
                sep='\t', index_col=0).loc[np.sort(algos_acronym.keys()), df_measurements.columns.values]
    df_measurements[np.logical_or(df_zeros == 0, np.isnan(df_zeros)).values] = np.nan

    dot_grid(df_measurements, "recall", "Recall", ax=axs[0][1])

    df_measurements=pd.read_csv(os.path.join(main_path, "recovery_results_{}_matrix_f1.tsv".format(suffix)), sep='\t', index_col=0).loc[np.sort(algos_acronym.keys()),:]

    df_zeros = pd.read_csv(os.path.join('/home/hag007/Desktop/aggregate{}_report/venn/'.format(omic_type), "count_matrix.tsv"),
                sep='\t', index_col=0).loc[np.sort(algos_acronym.keys()), df_measurements.columns.values]
    df_measurements[np.logical_or(df_zeros == 0, np.isnan(df_zeros)).values] = np.nan

    dot_grid(df_measurements, "f1", "F1", ax=axs[1][0])

    plt.figtext(0.01,0.99, "A:", weight='bold')
    plt.figtext(0.5, 0.99, "B:", weight='bold')
    plt.figtext(0.01, 0.5, "C:", weight='bold')
    # plt.figtext(0.5, 0.5, "D:", weight='bold')

    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_9_{}.png".format(suffix)))

if __name__=="__main__":
    main()