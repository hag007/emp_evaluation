import sys
sys.path.insert(0,'../../')
import matplotlib
matplotlib.use('Agg')
from plots.ehr_count_richness_dot_plots import main


if __name__=='__main__':
    main(["DOMINO4", "netbox2_string"])
