import sys
sys.path.insert(0,'../../')
from plots.summary_table import main
import matplotlib
matplotlib.use('Agg')

if __name__=="__main__":
    main(["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox", "keypathwayminer_INES_GREEDY", "DOMINO2", "hotnet2"])