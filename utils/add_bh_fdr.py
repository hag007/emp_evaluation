import sys
sys.path.append('../')
import pandas as pd
import numpy as np
import constants
import os
from statsmodels.sandbox.stats.multicomp import fdrcorrection0


def add_bh_correction(dataset_name):
    constants.update_dirs(DATASET_NAME_u=dataset_name)
    df_raw = pd.read_csv(os.path.join(constants.DATA_DIR, "score.tsv"), sep="\t", index_col=0)
    df_raw=df_raw[df_raw.index!="nan"]
    df_raw=df_raw[(df_raw.T != 0).any()]
    # df_raw=df_raw.dropna()
    df_raw= df_raw[~df_raw.index.duplicated(keep="first")]
    print df_raw.columns
    df_raw.loc[:,'qval'] =fdrcorrection0(df_raw.loc[:,'pval'].values)[1]
    df_raw.to_csv(os.path.join(constants.DATA_DIR, "score.tsv"), sep="\t", index_label="id")



if __name__ =='__main__':
    for dataset_name in ['PASCAL_SUM_Age_Related_Macular_Degeneration.G50','PASCAL_SUM_Atrial_Fibrillation.G50']:
        add_bh_correction(dataset_name)