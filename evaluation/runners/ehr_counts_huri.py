import sys
sys.path.insert(0, '../../')
import argparse
from evaluation.ehr_counts import main

if __name__ == "__main__":
    datasets=["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"]
    algos=['DOMINO3', 'netbox3']
    prefix="GE"
    main(datasets,algos,prefix)

    datasets=["brca", "crh", "scz", "tri", "t2d", "cad", "bmd", "hgt", "amd", "af"]
    algos=['DOMINO3', 'netbox3']
    prefix="PASCAL_SUM"
    main(datasets,algos,prefix)