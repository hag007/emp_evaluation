import sys
sys.path.insert(0,'../../')
import matplotlib
matplotlib.use("Agg")
from plots.robustness_empties_plot import main
import matplotlib.pyplot as plt
import os
import constants

if __name__=="__main__":
    algos=["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox", "keypathwayminer_INES_GREEDY", "DOMINO2", "hotnet2"]
    main(algos)
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "figure_19.png"))
