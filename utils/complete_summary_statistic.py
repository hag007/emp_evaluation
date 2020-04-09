import pandas as pd
import numpy as np
import os

import constants



prefix="PASCAL_SUM"
df_summary=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "summary_statistic_{}.tsv".format(prefix)), index_col=0, sep='\t')

for i, cur_line in df_summary.iterrows():
    n_modules=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"{}_{}".format(prefix, cur_line["dataset"][len(prefix)+1:]),cur_line['algo'],"modules_summary.tsv"), sep='\t').shape[0]
    df_summary.loc[(df_summary.index==i) & (df_summary.dataset==cur_line["dataset"]), 'n_modules']=n_modules
    print i, cur_line["dataset"], n_modules

df_summary.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "summary_statistic_{}.tsv".format(prefix)), sep='\t')


