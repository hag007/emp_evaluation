import sys
sys.path.insert(0,'../../')


import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import os
import constants

from plots.robustness_integrated_plot import plot_itegrated_recovery

if __name__=="__main__":

    base_folder=os.path.join(constants.OUTPUT_GLOBAL_DIR,"evaluation")
    auc_file_format = "robustness_auc_{}_{}.tsv"
    p_file_format = "robustness_{}_{}_matrix_p.tsv"
    r_file_format = "robustness_{}_{}_matrix_r.tsv"
    f1_file_format = "robustness_{}_{}_matrix_f1.tsv"
    zeros_file_format = os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation","count_matrix_{}.tsv")
    ss_ratios = [0.1, 0.2, 0.3, 0.4]

    fig,axs=plt.subplots(2,2,figsize=(18,16))
    prefix="GE"
    suffix = "{}_100".format(prefix)
    zeros_file_name =  zeros_file_format.format(prefix)
    algos = ["DOMINO2", "netbox"] # , "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    datasets = ["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"]
    plot_itegrated_recovery(ss_ratios, base_folder, auc_file_format, p_file_format, r_file_format, f1_file_format, zeros_file_name, suffix, axs=axs[:,0], title="GE", algos=algos, datasets=datasets)

    prefix = "PASCAL_SUM"
    suffix = "{}_100".format(prefix)
    zeros_file_name =  zeros_file_format.format(prefix)
    algos = ["DOMINO2", "netbox"] # , "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    datasets = ["brca", "crh", "scz", "tri", "t2d", "cad", "bmd", "hgt", "amd", "af"]
    plot_itegrated_recovery(ss_ratios, base_folder, auc_file_format, p_file_format, r_file_format,
                            f1_file_format, zeros_file_name, suffix, axs=axs[:,1], title="GWAS", algos=algos, datasets=datasets)

    plt.tight_layout()
    plt.figtext(0.01, 0.97, "A:", weight='bold', fontsize=22)
    plt.figtext(0.01, 0.5, "B:", weight='bold', fontsize=22)
    plt.figtext(0.5, 0.97, "C:", weight='bold', fontsize=22)
    plt.figtext(0.5, 0.5, "D:", weight='bold', fontsize=22)
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "figure_14.png"))




