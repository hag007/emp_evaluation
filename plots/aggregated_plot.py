import math
import random
import sys
sys.path.insert(0, '../')

import pandas as pd

from fastsemsim.SemSim import *
from fastsemsim.Ontology import ontologies
from fastsemsim.Ontology import AnnotationCorpus
from fastsemsim.SemSim.SetSemSim import SetSemSim

import matplotlib
matplotlib.use("Agg")
import seaborn as sns
from rpy2.robjects import pandas2ri
pandas2ri.activate()

import constants

import multiprocessing

from utils.daemon_multiprocessing import func_star

import argparse

import utils.go_hierarcies as go_hierarcies

import simplejson as json

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.lines import Line2D

def plot_mEHR(df, ax, max_n_modules_th=20, title="mEHR"):
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
    base_folder = os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr", "MAX")
    cutoffs = [1.0, 2.0, 3.0, 4.0]
    summary_m_ehr_mean = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "mEHR_mean_{}.tsv".format(prefix)),
                                     sep='\t', index_col=0)
    plot_mEHR(summary_m_ehr_mean, axs[0][0], title="mEHR")
    import richness_curve_plot
    richness_curve_plot.main(algos=algos, datasets=datasets, prefix=prefix, cutoffs=cutoffs, ax=axs[0][1],
                             title="Biological richness")
    homogeneity_file_format = "homogeneity_avg_matrix_{}_{}.tsv"
    homogeneity_std_file_format = "homogeneity_std_matrix_{}_{}.tsv"
    heterogeneity_file_format = "heterogeneity_avg_matrix_{}_{}.tsv"
    import homogeneity_plot
    homogeneity_plot.variability_plot(algos, datasets, constants.OUTPUT_GLOBAL_DIR, homogeneity_file_format,
                                      homogeneity_std_file_format, heterogeneity_file_format, prefix, cutoffs,
                                      axs[1][0], None, title="Intra-module homogeneity")
    base_folder = "/media/hag007/Data/bnet/output/emp_fdr/MAX"
    auc_file_format = "pr_auc_recovery_summary_{}_{}.tsv"
    p_file_format = "recovery_results_{}_{}_matrix_p.tsv"
    r_file_format = "recovery_results_{}_{}_matrix_r.tsv"
    f1_file_format = "recovery_results_{}_{}_matrix_f1.tsv"
    empty_file_format = "recovery_results_{}_{}_matrix_empty.tsv"
    average_file_format = "recovery_results_{}_{}.tsv"
    zeros_file_format = os.path.join('/home/hag007/Desktop/aggregate{}_report/venn', "count_matrix.tsv")
    ss_ratios = [0.1, 0.2, 0.3, 0.4]
    suffix = "{}_100".format(prefix)
    omic_type = ""
    zeros_file_name = zeros_file_format.format(omic_type)
    import robustness_integrated_plot
    robustness_integrated_plot.plot_itegrated_recovery(ss_ratios, base_folder, average_file_format, auc_file_format,
                                                       p_file_format, r_file_format, f1_file_format, empty_file_format,
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
    fig.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "aggregated_plot_{}.png".format(prefix)))


if __name__ == "__main__":

    prefix = "GE" # "GE" # PASCAL_SUM
    algos = ["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox",  "domino_original", "keypathwayminer_INES_GREEDY"] # "keypathwayminer_INES_GREEDY",
    datasets = ["TNFa_2", "HC12", "SHERA", "SHEZH_1", "ROR_1", "ERS_1", "IEM" , "APO", "CBX", "IFT"]# ["Breast_Cancer.G50", "Crohns_Disease.G50", "Schizophrenia.G50", "Triglycerides.G50", "Type_2_Diabetes.G50" ,"Coronary_Artery_Disease.G50" , "Bone_Mineral_Density.G50", "Height1.G50", "Age_Related_Macular_Degeneration.G50", "Atrial_Fibrillation.G50"] # ["TNFa_2", "HC12", "SHERA", "SHEZH_1", "ROR_1", "ERS_1", "IEM" , "APO", "CBX", "IFT"]
    main(prefix, algos, datasets)

    prefix = "PASCAL_SUM" # "GE" # PASCAL_SUM
    algos = ["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox",  "domino_original"] # "keypathwayminer_INES_GREEDY",
    datasets = ["Breast_Cancer.G50", "Crohns_Disease.G50", "Schizophrenia.G50", "Triglycerides.G50", "Type_2_Diabetes.G50" ,"Coronary_Artery_Disease.G50" , "Bone_Mineral_Density.G50", "Height1.G50", "Age_Related_Macular_Degeneration.G50", "Atrial_Fibrillation.G50"]
    main(prefix, algos, datasets)
