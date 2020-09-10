import sys
sys.path.insert(0, '../../')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fastsemsim.SemSim import *
import constants
from plots.homogeneity_plot import variability_plot

if __name__=="__main__":

    base_folder=os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation")
    homogeneity_file_format = "homogeneity_avg_matrix_{}_{}.tsv"
    cutoffs = [1.0, 2.0, 3.0, 4.0]

    fig,axs=plt.subplots(1,2, figsize=(22,10))
    fig_violin, axs_violin = plt.subplots(2, len(cutoffs), figsize=(4 * len(cutoffs) * 2, 10))

    datasets = ["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"]
    algos = ["DOMINO2", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    prefix = "GE"
    variability_plot(algos, datasets, base_folder, homogeneity_file_format, prefix, cutoffs, axs[0], axs_violin[0], title="GE")

    datasets = ["brca", "crh", "scz", "tri", "t2d", "cad", "bmd", "hgt", "amd", "af"]
    algos = ["DOMINO2", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    prefix = "PASCAL_SUM"
    variability_plot(algos, datasets, base_folder, homogeneity_file_format, prefix, cutoffs, axs[1], axs_violin[1], title="GWAS")

    fig.text(0.01,0.97, "A:", weight='bold',fontsize=22)
    fig.text(0.5, 0.97, "B:", weight='bold',fontsize=22)
    fig.tight_layout()
    fig.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "figure_15.png"))

