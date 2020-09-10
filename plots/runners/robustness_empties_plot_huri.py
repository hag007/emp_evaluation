import sys
sys.path.insert(0,'../../')
from plots.robustness_empties_plot import main
import matplotlib.pyplot as plt
import os
import constants

if __name__=="__main__":
    algos=["DOMINO3", "netbox3"]
    main(algos)
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "figure_19.png"))
