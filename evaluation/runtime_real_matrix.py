import sys
sys.path.insert(0, '../')

import constants
import os
import argparse
import pandas as pd
import numpy as np
from statsmodels.sandbox.stats.multicomp import fdrcorrection0

def get_timing(timing_file_name):
   
    try:
        return float(open(timing_file_name,'r').read())
    except Exception, e:
        print e
        return 9999 
    
def main(datasets, algos, indices, network, prefix):

    df_summary = pd.DataFrame()
    for cur_ds in datasets:
        for cur_alg in algos:
            timing=get_timing("/specific/netapp5/gaga/hagailevi/emp_test/timings/{}/{}/{}.txt".format(network,cur_ds,cur_alg))
            df_summary.loc[cur_alg,cur_ds]=timing

    df_summary.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation", "timing_real_{}_{}.tsv".format(network, prefix)), sep='\t')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="TNFa_2,HC12,ROR_1,SHERA,SHEZH_1,ERS_1,IEM,APO,CBX,IFT") # "TNFa_2,HC12,ROR_1,SHERA,SHEZH_1,ERS_1,IEM,APO,CBX,IFT" Breast_Cancer.G50,Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50
    parser.add_argument('--prefix', dest='prefix', default="GE") # PASCAL_SUM   GE
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy,jactivemodules_sa,bionet,netbox,keypathwayminer_INES_GREEDY,dcem2") # ,dcem2,dcem3,dcem4,my_netbox_td,hotnet2

    args = parser.parse_args()
    datasets = args.datasets.split(",")
    algos = args.algos.split(",")
    prefix = args.prefix
    
    datasets=["brca", "crh", "scz", "tri", "t2d", "bmd", "amd", "af", "hgt", "cad"]
    algos=["netbox", "DOMINO2" , "jactivemodules_greedy", "jactivemodules_sa", "bionet", "hotnet2", "keypathwayminer_INES_GREEDY"]
    network="dip"
    prefix="PASCAL_SUM"
    indices=np.arange(5100,5110)
    main(datasets,algos,indices,network,prefix)

    datasets=["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"]
    algos=["netbox", "DOMINO2", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "hotnet2", "keypathwayminer_INES_GREEDY"]
    network="dip"    
    prefix="GE"
    indices=np.arange(5100,5110)
    main(datasets,algos,indices,network,prefix)

    datasets=["brca", "crh", "scz", "tri", "t2d", "bmd", "amd", "af", "hgt", "cad"]
    algos=["netbox2_string", "DOMINO2"]
    network="huri"
    prefix="PASCAL_SUM"
    indices=np.arange(5100,5110)
    main(datasets,algos,indices,network,prefix)

    datasets=["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"]
    algos=["netbox2_string", "DOMINO2"]
    network="huri"    
    prefix="GE"
    indices=np.arange(5100,5110)
    main(datasets,algos,indices,network,prefix)

 
    datasets=["brca", "crh", "scz", "tri", "t2d", "bmd", "amd", "af", "hgt", "cad"]
    algos=["netbox2_string", "DOMINO2"]
    network="string_th_901"
    prefix="PASCAL_SUM"
    indices=np.arange(5100,5110)
    main(datasets,algos,indices,network,prefix)

    datasets=["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"]
    algos=["netbox2_string", "DOMINO2"]
    network="string_th_901"    
    prefix="GE"
    indices=np.arange(5100,5110)
    main(datasets,algos,indices,network,prefix)
# 
