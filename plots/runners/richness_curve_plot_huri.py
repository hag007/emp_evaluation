import sys
sys.path.insert(0, '../../')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import constants
import os
from plots.richness_curve_plot import main


if __name__ == "__main__":

    cutoffs=[1.0,2.0,3.0,4.0]

    fig, axs = plt.subplots(1,2,figsize=(20,10))
    fig_violin, axs_violin = plt.subplots(2,len(cutoffs),figsize=(4*len(cutoffs)*2,10))
    datasets = ["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"]
    algos = ["DOMINO3", "netbox3"]
    prefix = "GE"
    main(algos=algos, datasets=datasets, prefix=prefix, cutoffs=cutoffs, ax=axs[0], axs_violin=axs_violin[0], title="GE")

    datasets = ["brca", "crh", "scz", "tri", "t2d", "cad", "bmd", "hgt", "amd", "af"]
    algos = ["DOMINO3", "netbox3"]
    prefix = "PASCAL_SUM"
    main(algos=algos, datasets=datasets, prefix=prefix, cutoffs=cutoffs, ax=axs[1], axs_violin=axs_violin[1], title="GWAS")

    fig.text(0.01, 0.97, "A:", weight='bold', fontsize=22)
    fig.text(0.5, 0.97, "B:", weight='bold', fontsize=22)
    fig.tight_layout()
    fig.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots","figure_13.png"))