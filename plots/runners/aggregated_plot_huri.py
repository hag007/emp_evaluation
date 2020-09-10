import sys
sys.path.insert(0, '../../')

import matplotlib
matplotlib.use("Agg")

from plots.aggregated_plot import main


if __name__ == "__main__":

    prefix = "GE"
    algos = ['DOMINO3', 'netbox3']
    datasets = ["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"]
    main(prefix, algos, datasets)

    prefix = "PASCAL_SUM"
    algos = ['DOMINO3', 'netbox3']
    datasets = ["brca", "crh", "scz", "tri", "t2d", "cad", "bmd", "hgt", "amd", "af"]
    main(prefix, algos, datasets)
