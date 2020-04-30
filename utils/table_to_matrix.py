import pandas as pd
import numpy as np
import os
import constants

def table2matrix(df_table, row_name, col_name, value_name):
    df_matrix=pd.DataFrame()
    for index, row in df_table.iterrows():
        df_matrix.loc[row[row_name], row[col_name]]= row[value_name]

    return df_matrix

if __name__=='__main__':

    for cur_prefix in ["GE", "PASCAL_SUM"]:
        for cur_ss_ratio in [0.4, 0.3, 0.2, 0.1]:
            suffix="{}_100_{}".format(cur_prefix,cur_ss_ratio)
            # table2matrix(pd.read_csv('/home/hag007/Desktop/recovery_results_PASCAL_SUM_100.tsv',sep='\t'), 'algo', 'dataset', 'p_mean').to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "recovery_results_PASCAL_SUM_100_matrix_p.tsv"),sep='\t')
            # table2matrix(pd.read_csv('/home/hag007/Desktop/recovery_results_PASCAL_SUM_100.tsv',sep='\t'), 'algo', 'dataset', 'r_mean').to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "recovery_results_PASCAL_SUM_100_matrix_r.tsv"),sep='\t')
            # table2matrix(pd.read_csv('/home/hag007/Desktop/recovery_results_PASCAL_SUM_100.tsv',sep='\t'), 'algo', 'dataset', 'f1_mean').to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "recovery_results_PASCAL_SUM_100_matrix_f1.tsv"),sep='\t')

            base_folder="/media/hag007/Data/bnet/output/emp_fdr/MAX"
            df=pd.read_csv(os.path.join(base_folder, 'recovery_results_{}.tsv'.format(suffix)), sep='\t')
            table2matrix(df, 'algo', 'dataset', 'p_mean').to_csv(os.path.join(base_folder, "recovery_results_{}_matrix_p.tsv".format(suffix)),sep='\t')
            table2matrix(df, 'algo', 'dataset', 'r_mean').to_csv(os.path.join(base_folder, "recovery_results_{}_matrix_r.tsv".format(suffix)),sep='\t')
            table2matrix(df, 'algo', 'dataset', 'f1_mean').to_csv(os.path.join(base_folder, "recovery_results_{}_matrix_f1.tsv".format(suffix)),sep='\t')
            # df['non_empty_solutions']=df.apply(lambda a: np.sum(np.logical_and(np.array(a["precisions"][1:-1].split(','),dtype=np.float)!= 0, np.array(a["recalls"][1:-1].split(','),dtype=np.float) != 0)), axis=1)
            # table2matrix(df, 'algo', 'dataset', 'non_empty_solutions').to_csv(os.path.join(base_folder, "recovery_results_{}_matrix_empty.tsv".format(suffix)), sep='\t')

