import sys
sys.path.insert(0,'../../')
from plots.summary_table import main
import matplotlib
matplotlib.use('Agg')

if __name__=="__main__":
    main(["netbox3", "DOMINO3"])