import sys
sys.path.insert(0, '../')

import pandas as pd

from fastsemsim.SemSim import *

import matplotlib
matplotlib.use("Agg")

import constants

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def plot_mEHR(df, algos, ax, max_n_modules_th=20, title="mEHR"):
    for i, cur_row in df.loc[constants.ALGOS].iterrows():

        ax.scatter(np.arange(1, max_n_modules_th + 1), cur_row,c=constants.COLORDICT[i])
        
    ax.set_xlabel("# top ranked modules", fontsize=25)
    ax.set_ylabel("average mEHR", fontsize=25)
    # ax.set_title("average EHR as function of modules' head threshold")
    ax.set_facecolor('#ffffff')
    ax.grid(color='gray')
    patches=[Line2D([0], [0], marker='o', color='gray', label=constants.ALGOS_ACRONYM[a], markersize=12, markerfacecolor=constants.COLORDICT[a], alpha=0.7) for a in algos]
    ax.legend(handles=patches, fontsize=22, loc=(0,1.1), ncol=3, facecolor='#ffffff')
    ax.set_title("mEHR", size=25)


def main(prefix, algos, datasets):
    fig,axs=plt.subplots(3,2,figsize=(20,30))
    base_folder = os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation")
    cutoffs = [1.0, 2.0, 3.0, 4.0]
    summary_m_ehr_mean = pd.read_csv(os.path.join(base_folder, "mEHR_mean_{}.tsv".format(prefix)),
                                     sep='\t', index_col=0)
    plot_mEHR(summary_m_ehr_mean, algos, axs[0][0], title="mEHR")
    import richness_curve_plot
    richness_curve_plot.main(algos=algos, datasets=datasets, prefix=prefix, cutoffs=cutoffs, ax=axs[0][1],
                             title="Biological richness")
    homogeneity_file_format = "homogeneity_avg_matrix_{}_{}.tsv"
    homogeneity_std_file_format = "homogeneity_std_matrix_{}_{}.tsv"
    heterogeneity_file_format = "heterogeneity_avg_matrix_{}_{}.tsv"
    import homogeneity_plot
    homogeneity_plot.variability_plot(algos, datasets, base_folder, homogeneity_file_format,
                                       prefix, cutoffs,
                                      axs[1][0], None, title="Intra-module homogeneity")
    auc_file_format = "robustness_auc_{}_{}.tsv"
    p_file_format = "robustness_{}_{}_matrix_p.tsv"
    r_file_format = "robustness_{}_{}_matrix_r.tsv"
    f1_file_format = "robustness_{}_{}_matrix_f1.tsv"
    empty_file_format = "robustness_{}_{}_matrix_empty.tsv"
    average_file_format = "robustness_{}_{}.tsv"
    zeros_file_format = os.path.join(base_folder, "count_matrix_{}.tsv")
    ss_ratios = [0.1, 0.2, 0.3, 0.4]
    suffix = "{}_100".format(prefix)
    zeros_file_name = zeros_file_format.format(prefix)
    import robustness_integrated_plot
    robustness_integrated_plot.plot_itegrated_recovery(ss_ratios, base_folder, auc_file_format,
                                                       p_file_format, r_file_format, f1_file_format,
                                                       zeros_file_name, suffix, axs=[axs[1, 1], axs[2, 0]],
                                                       title="Robustness", algos=algos, datasets=datasets)
    axs[2][1].axis('off')
    axs[0][0].legend([])
    axs[0][1].legend([])
    axs[1][0].legend([])
    axs[1][1].legend([])
    axs[2][0].legend([])
    patches = [Line2D([0], [0], marker='o', color='gray', label=constants.ALGOS_ACRONYM[a], markersize=12,
                      markerfacecolor=constants.COLORDICT[a], alpha=0.7) for a in algos]
    fig.legend(handles=patches, fontsize=22, loc='lower right', ncol=3, facecolor='#ffffff')
    fig.tight_layout()
    fig.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "aggregated_plot_{}.png".format(prefix)))


if __name__ == "__main__":

    prefix = "GE" # "GE" # PASCAL_SUM
    algos = ['netbox2_string', 'DOMINO4'] # ["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox",  "DOMINO", "keypathwayminer_INES_GREEDY", "hotnet2"] # "keypathwayminer_INES_GREEDY",
    datasets = ["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"]
    main(prefix, algos, datasets)

    prefix = "PASCAL_SUM" # "GE" # PASCAL_SUM
    algos = ['netbox2_string', 'DOMINO4'] # ["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox",  "DOMINO", "keypathwayminer_INES_GREEDY", "hotnet2"]
    datasets = ["brca", "crh", "scz", "tri", "t2d", "cad", "bmd", "hgt", "amd", "af"]
    main(prefix, algos, datasets)
