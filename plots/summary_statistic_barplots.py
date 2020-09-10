import sys
sys.path.insert(0, "../")

import matplotlib
matplotlib.use("Agg")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import constants
import os

ENABLE_GO_GRAPH = False
IS_GO_GRAPH_ONLY=False
GO_PCA=False
TERMS_SIMILARITY_TO_NUM_OF_TERMS=True
RATIO_TO_GO_TERM = True
QVAL_TH = 0.01
SIM_TH= 0.4

n_col=3

def align_yaxis(ax1, v1, ax2, v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny+dy, maxy+dy)


def barplot(df_summary_statistic, prefix, algos, dataset, ax, has_label=False):
    df_summary_statistic=df_summary_statistic[df_summary_statistic['dataset'] == dataset]
    algos=list(df_summary_statistic.loc[:,'algo'])
    objects = tuple(algos)
    y_pos = np.arange(len(objects))
    avg_module_size = df_summary_statistic.loc[:,'module_size_mean']
    std_module_size = df_summary_statistic.loc[:,'module_size_std']
    n_modules = df_summary_statistic.loc[:,'n_modules']
    ax.set_facecolor('#ffffff')
    ax.grid(False)
    ax2 = ax.twinx()
    ax2.grid(False)

    label1='_nolegend_'
    label2 = '_nolegend_'
    if has_label:
        label1="average module size"
        label2="# of modules"

    ax.bar(y_pos-0.1, avg_module_size, width=0.2, align='center', alpha=0.5, ecolor='black', color='red', label=label1 )
    err1 = ax.errorbar(y_pos-0.1, avg_module_size, yerr=np.float32(std_module_size), lolims=True, capsize=5, ls='None', color='black', elinewidth=0.5, label='_nolegend_')
    err1[1][0].set_marker('_')
    err1[1][0].set_markersize(10)
    ax2.bar(y_pos+0.1, n_modules, width=0.2, align='center', alpha=0.5, color='blue', label=label2)

    ax.set_xticks(y_pos)
    ax.set_xticklabels([constants.ALGOS_ACRONYM[a] for a in objects], rotation='45', fontsize=16)
    # ax.set_xlabel('alg', fontsize=16)
    ax.set_ylabel('average module size', fontsize=16)
    ax2.set_ylabel('# of modules', fontsize=16)
    ax.set_title('dataset: {}'.format(dataset), fontsize=16)

    align_yaxis(ax, 0, ax2, 0)
    raw_factor=1.2
    ax1_factor=max((max(avg_module_size) / ax.get_ylim()[1]) * raw_factor, raw_factor)

    ax2_factor=max((max(n_modules) / ax2.get_ylim()[1]) * raw_factor, raw_factor)

    ax_factor=max(ax1_factor, ax2_factor)

    print ax_factor

    ax2.get_ylim()[1]

    ylim=ax2.get_ylim()
    print ylim
    y_lim_factor = 1.2
    ylim=(ylim[0] * y_lim_factor * max(n_modules)/(ylim[1] * ax_factor), y_lim_factor * max(n_modules))
    ax2.set_ylim(ylim)
    print ylim

    y_lim_factor = 1.2
    ylim = ax.get_ylim()
    print ylim
    ylim_1 = (ylim[0] * y_lim_factor * (np.nanmax(avg_module_size) + np.nanmax(std_module_size)) /(ylim[1] * ax_factor),  y_lim_factor* (np.nanmax(avg_module_size) + np.nanmax(std_module_size)) )
    ax.set_ylim(ylim_1)
    print ylim_1


    for algo, bp in zip(algos, avg_module_size):
        y_pos_n_genes=ax.get_ylim()[1]-(ax.get_ylim()[1]-ax.get_ylim()[0])*0.1 #
        txt="{}".format(int(df_summary_statistic.loc[(df_summary_statistic.loc[:,"algo"]==algo) & (df_summary_statistic.dataset==dataset), "n_genes"][0])) # # genes:
        ax.text(algos.index(algo)-0.3, y_pos_n_genes, txt, color='green', fontsize=16)



def barplots(df_summary_statistic, prefix, axs):

    algos =  np.sort(np.unique(df_summary_statistic.index.values))
    datasets = np.sort(np.unique(df_summary_statistic['dataset'].values))

    for i, dataset in enumerate(datasets):
        has_label = i==0
        try:
            barplot(df_summary_statistic, prefix, algos, dataset, axs[i/n_col, np.mod(i,n_col)], has_label) #
        except Exception, e:
            print e
            pass

    for j in np.arange(i+1,axs.size): #
        axs[j / n_col, np.mod(j, n_col)].set_axis_off()



def main():

    fig, axs = plt.subplots(int(np.ceil(10.0/n_col)), n_col, figsize=(20, 10.0*int(np.ceil(7.0/n_col))))
    # fig, axs = plt.subplots(1,1, figsize=(20, 20))

    main_path = constants.OUTPUT_GLOBAL_DIR
    prefix="GE"
    df_statistic=pd.read_csv(os.path.join(main_path, "evaluation", "summary_statistics_{}.tsv".format(prefix)), sep='\t', index_col=0)
    barplots(df_statistic, prefix, axs=axs)
    plt.figlegend(loc=(0.84,0.05),
                        facecolor='#ffffff', ncol=1, prop={'size': 14})
    fig.tight_layout()
    plt.subplots_adjust(left=0.07, right=0.95, top=0.95, bottom=0.07)
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "figure_18_{}.png".format(prefix)))

    fig, axs = plt.subplots(int(np.ceil(10.0 / 3)), 3, figsize=(20, 10.0*int(np.ceil(7.0/n_col))))
    prefix="PASCAL_SUM"
    df_statistic=pd.read_csv(os.path.join(main_path, "evaluation", "summary_statistics_{}.tsv".format(prefix)), sep='\t', index_col=0)
    barplots(df_statistic, prefix, axs=axs)
    plt.figlegend(loc=(0.84,0.05),
                        facecolor='#ffffff', ncol=1, prop={'size': 14})
    fig.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.1)
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "figure_18_{}.png".format(prefix)))



if __name__=="__main__":
    main()
