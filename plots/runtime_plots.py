import sys
sys.path.insert(0,'../')
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

def dot_plot(df_summary, y_label, ax=None, algos=None, title="", df_avg=None):
    plt.minorticks_off()
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))

    ax.set_facecolor('#ffffff')
    ax.grid(color='gray')
    inf_runtimes={}
    for cur_alg in algos:
        df_algo=df_summary[df_summary.loc[:, "algo"]==cur_alg]
        inf_runtimes[cur_alg]=df_algo.loc[:, "runtime"][df_algo.loc[:, "runtime"]==9999].shape[0]

    df_summary=df_summary[df_summary.loc[:, "runtime"]!=9999]
    df_summary.loc[:, "Alg"]=df_summary.loc[:, "algo"].apply(lambda a: constants.ALGOS_ACRONYM[a])
    df_summary=df_summary.rename(columns={"runtime": "Seconds"})
    my_order = [constants.ALGOS_ACRONYM[a] for a in algos]
    ax = sns.stripplot(x='Alg', y=y_label, data=df_summary, jitter=True, size=10, ax=ax, order=my_order, palette={a: "#AAAAAA" for a in my_order}) # {a: colorlist[algo_acrs.index(a)] for a in my_order}

    patches = [Line2D([0], [0], marker='s', color='gray', label='mean', markersize=12, markerfacecolor='blue', alpha=0.7),  Line2D([0], [0], marker='P', color='gray', label='median', markersize=12, markerfacecolor='red', alpha=0.7)]
    i=0
    for i_alg, cur_alg in enumerate(algos):
    
        vals=df_summary[df_summary.loc[:, "algo"]==cur_alg].loc[:,"Seconds"].values.reshape(-1)
        mn = np.nanmean(vals)
        std = np.nanstd(vals)
        mdn = np.median(vals)
        mx = np.nanmax(list(vals)+[0])
        df_avg.loc[cur_alg,"mean"]=mn
        df_avg.loc[cur_alg,"# long runs"]=inf_runtimes[cur_alg]    
        ax.bar([list(my_order).index(constants.ALGOS_ACRONYM[cur_alg])], [mn], color='blue')
        err1 = ax.errorbar([list(my_order).index(constants.ALGOS_ACRONYM[cur_alg])-0.1], [mn], yerr=np.float32(std), lolims=True, capsize=5, ls='None', color='black', elinewidth=0.5, label='_nolegend_')
        err1[1][0].set_marker('_')
        err1[1][0].set_markersize(10)
        ax.text(list(my_order).index(constants.ALGOS_ACRONYM[cur_alg]), mx, "inf={}".format(inf_runtimes[cur_alg]))
        ax.scatter(list(my_order).index(constants.ALGOS_ACRONYM[cur_alg]), [mdn], marker='P', color='red', edgecolors='b', s=[200], alpha=0.7, zorder=5)
        i+=1

    ax.set_xlabel(ax.get_xlabel(), fontdict={"size" : 22})
    ax.set_ylabel(ax.get_ylabel(), fontdict={"size": 22})
    ax.set_title(title, fontdict={"size":22})
    ax.legend(handles=patches, loc=(0.0,1.0), fontsize=20, facecolor='#ffffff')
    ax.set_xticklabels([a for a in ax.get_xticklabels()], rotation=45)


def main(algos, network):

    fig_4, axs_4 = plt.subplots(1,2, figsize=(30, 15))
    main_path = os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation")

    prefix="GE"
    df_avg_ge=pd.DataFrame()
    df_summary=pd.read_csv(os.path.join(main_path, "timing_{}_{}.tsv".format(network, prefix)), sep='\t', index_col=0)
    dot_plot(df_summary, "Seconds", ax=axs_4[0], algos=algos, title="Alg. running times, {}".format(prefix), df_avg=df_avg_ge)
    df_avg_ge["omics"]=prefix
    df_avg_ge["network"]=network

    prefix = "PASCAL_SUM"
    df_avg_gwas=pd.DataFrame()
    df_summary = pd.read_csv(os.path.join(main_path, "timing_{}_{}.tsv".format(network, prefix)), sep='\t', index_col=0)
    dot_plot(df_summary, "Seconds", ax=axs_4[1], algos=algos, title="Alg. running times, {}".format(prefix), df_avg=df_avg_gwas)
    df_avg_gwas["omics"]=prefix
    df_avg_gwas["network"]=network

    fig_4.tight_layout()
    fig_4.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "figure_21_{}.png".format(network)))
    
    return pd.concat([df_avg_ge,df_avg_gwas])


if __name__=='__main__':
    df_dip=pd.DataFrame()
    df_huri=pd.DataFrame()
    df_string_th_901=pd.DataFrame()

    df_dip=main(["netbox", "DOMINO2", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "hotnet2", "keypathwayminer_INES_GREEDY"], "dip")
    df_huri=main(["netbox", "DOMINO2"], "huri")
    df_string_th_901=main(["netbox", "DOMINO2"], "string_th_901")
 
    pd.concat([df_dip,df_huri,df_string_th_901]).to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation", "df_runtime_avg.tsv"), sep='\t')
