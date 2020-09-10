import sys
sys.path.insert(0,'../../')
import matplotlib
matplotlib.use('Agg')
from plots.ehr_count_richness_dot_plots import main


if __name__=='__main__':
    main(["DOMINO2", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"])
