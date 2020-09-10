import sys
sys.path.insert(0, '../../')
import argparse
from evaluation.ehr_counts import main

if __name__ == "__main__":
    datasets=["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"]
    algos=['DOMINO4', 'netbox2_string']
    prefix="GE"
    main(datasets,algos,prefix)

    datasets=["brca", "crh", "scz", "tri", "t2d", "cad", "bmd", "hgt", "amd", "af"]
    algos=['DOMINO4', 'netbox2_string']
    prefix="PASCAL_SUM"
    main(datasets,algos,prefix)