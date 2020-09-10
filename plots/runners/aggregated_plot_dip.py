import sys
sys.path.insert(0, '../../')

import matplotlib
matplotlib.use("Agg")

from plots.aggregated_plot import main


if __name__ == "__main__":

    prefix = "GE"
    algos = ["DOMINO2", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox", "keypathwayminer_INES_GREEDY", "hotnet2"]
    datasets = ["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"]
    main(prefix, algos, datasets)

    prefix = "PASCAL_SUM"
    algos = ["DOMINO2", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox", "keypathwayminer_INES_GREEDY", "hotnet2"]
    datasets = ["brca", "crh", "scz", "tri", "t2d", "cad", "bmd", "hgt", "amd", "af"]
    main(prefix, algos, datasets)
