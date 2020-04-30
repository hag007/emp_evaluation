from scipy.stats import pearsonr
import pandas as pd
import constants
import os
import numpy as np

def extract_pvals(dataset, algo,index):
    return pd.read_csv(os.path.join(constants.DATASETS_DIR, "GE_random_{}_{}_{}".format(dataset,algo,index),"cache",'deg_edger.tsv'),index_col=0,sep='\t').loc[:,'pval']

def extract_first_sample(dataset, algo,index):
    return pd.read_csv(os.path.join(constants.DATASETS_DIR, "GE_random_{}_{}_{}".format(dataset,algo,index),"data",'ge.tsv'),index_col=0,sep='\t').iloc[:,0]


dataset="TNFa_2"
algo="jactivemodules_greedy"
agg= pd.read_csv(os.path.join(constants.DATASETS_DIR, "GE_random_{}_{}_1".format(dataset,algo),"data","ge.tsv"),index_col=0,sep='\t').iloc[:,0]
for a in np.arange(3):
    temp=pd.concat([agg,extract_first_sample(dataset,algo,a)], axis=1)

    print pearsonr(temp.iloc[:,0], temp.iloc[:,1])

# print agg



