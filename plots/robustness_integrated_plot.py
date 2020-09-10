import sys
sys.path.insert(0,'../')

import matplotlib
matplotlib.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import constants

from matplotlib.lines import Line2D


def plot_itegrated_recovery(ss_ratios, base_folder, auc_file_format, p_file_format, r_file_format, f1_file_format, zero_file_name, suffix, axs, title="", algos=constants.ALGOS_ACRONYM.keys(), datasets=[]):

    axs[0].set_facecolor('#ffffff')
    axs[0].grid(color='gray')
    axs[1].set_facecolor('#ffffff')
    axs[1].grid(color='gray')

    ps = pd.DataFrame()
    rs = pd.DataFrame()
    f1s = pd.DataFrame()
    prs = pd.DataFrame()
    empties = pd.DataFrame()

    for cur_ss_ratio in ss_ratios:
        print pd.read_csv(os.path.join(base_folder, p_file_format.format(suffix, cur_ss_ratio)), sep='\t',
                        index_col=0).index

        df_zeros = pd.read_csv(zero_file_name, sep='\t', index_col=0).loc[pd.read_csv(os.path.join(base_folder, p_file_format.format(suffix, cur_ss_ratio)), sep='\t',
                        index_col=0).index]
        df_zeros.loc[:,:] = 1
        n_iterations = 3 
        print "n_iteration for ss_ratio={}: {}".format(cur_ss_ratio, n_iterations)
        df_cur_pr = pd.read_csv(os.path.join(base_folder, auc_file_format.format(suffix, cur_ss_ratio)), sep='\t',
                                index_col=0).loc[:,datasets]

        df_cur_pr[pd.isnull(df_cur_pr)]=0
        df_cur_pr=df_cur_pr.mean(axis=1)
        df_cur_p = pd.read_csv(os.path.join(base_folder, p_file_format.format(suffix, cur_ss_ratio)), sep='\t',
                               index_col=0).loc[:,datasets]

        df_cur_p[pd.isnull(df_cur_p)]=0
        df_cur_p=df_cur_p.mean(axis=1)
        df_cur_r = pd.read_csv(os.path.join(base_folder, r_file_format.format(suffix, cur_ss_ratio)), sep='\t',
                               index_col=0).loc[:,datasets]

        df_cur_r[pd.isnull(df_cur_r)]=0
        df_cur_r=df_cur_r.mean(axis=1)
        df_cur_f1 = pd.read_csv(os.path.join(base_folder, f1_file_format.format(suffix, cur_ss_ratio)), sep='\t',
                                index_col=0).loc[:,datasets]

        df_cur_f1[pd.isnull(df_cur_f1)]=0
        df_cur_f1=df_cur_f1.mean(axis=1)


        prs = pd.concat([prs, df_cur_pr], axis=1)
        ps = pd.concat([ps, df_cur_p], axis=1)
        rs = pd.concat([rs, df_cur_r], axis=1)
        f1s = pd.concat([f1s, df_cur_f1], axis=1)

    prs = prs.loc[np.sort(ps.index.values)]
    ps = ps.loc[np.sort(ps.index.values)]
    rs = rs.loc[np.sort(ps.index.values)]
    f1s = f1s.loc[np.sort(ps.index.values)]


    patches_0 = [Line2D([0], [0], marker='o', color='gray', label=constants.ALGOS_ACRONYM[a], markersize=12, markerfacecolor=constants.COLORDICT[a], alpha=0.7) for a in algos] # + \

    patches_1 = [Line2D([0], [0], marker='o', color='gray', label=constants.ALGOS_ACRONYM[a], markersize=12, markerfacecolor=constants.COLORDICT[a], alpha=0.7) for a in algos]
    ss_ratios = 1 - np.array(ss_ratios)

#     prs.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"prs_{}.tsv".format(title)), sep='\t')
#     empties.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"empties_{}.tsv".format(title)), sep='\t')
#     f1s.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"f1s_{}.tsv".format(title)), sep='\t')

    for cur in set(prs.index).intersection(constants.ALGOS):
        axs[0].plot([a for a,b in zip(ss_ratios,prs.loc[cur]) if not np.isnan(b)], [a for a in prs.loc[cur] if not np.isnan(a) ], c=constants.COLORDICT[cur])
        axs[1].plot(ss_ratios, f1s.loc[cur], c=constants.COLORDICT[cur], linewidth=2.0)

    axs[0].set_xlabel("subsample fraction", fontsize=22)
    axs[0].set_ylabel("AUPR", fontsize=22)
    axs[0].set_title(title, fontdict={"size":22})

    axs[1].set_xlabel("subsample fraction", fontsize=22)
    axs[1].set_ylabel("F1", fontsize=22)
    axs[1].set_title(title, fontdict={"size":22})

    axs[0].legend(handles=patches_0, loc=(0,1.1), ncol=3, fontsize=20, facecolor='#ffffff')
    axs[1].legend(handles=patches_1, loc=(0,1.1), ncol=3, fontsize=20, facecolor='#ffffff')



if __name__=="__main__":

    base_folder=os.path.join(constants.OUTPUT_GLOBAL_DIR,"evaluation")
    auc_file_format = "robustness_auc_{}_{}.tsv"
    p_file_format = "robustness_{}_{}_matrix_p.tsv"
    r_file_format = "robustness_{}_{}_matrix_r.tsv"
    f1_file_format = "robustness_{}_{}_matrix_f1.tsv"
    empty_file_format = "robustness_{}_{}_matrix_empty.tsv"
    average_file_format="robustness_{}_{}.tsv"
    zeros_file_format = os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation","count_matrix_{}.tsv")
    ss_ratios = [0.1, 0.2, 0.3, 0.4]

    fig,axs=plt.subplots(2,2,figsize=(18,16))
    prefix="GE"
    suffix = "{}_3".format(prefix)
    zeros_file_name =  zeros_file_format.format(prefix)
    algos = ['DOMINO4', 'netbox2_string']
    datasets = ["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"]
    plot_itegrated_recovery(ss_ratios, base_folder, auc_file_format, p_file_format, r_file_format, f1_file_format, zeros_file_name, suffix, axs=axs[:,0], title="GE", algos=algos, datasets=datasets)

    prefix = "PASCAL_SUM"
    suffix = "{}_3".format(prefix)
    zeros_file_name =  zeros_file_format.format(prefix)
    algos = ['DOMINO4', 'netbox2_string']
    datasets = ["brca", "crh", "scz", "tri", "t2d", "cad", "bmd", "hgt", "amd", "af"]
    plot_itegrated_recovery(ss_ratios, base_folder, auc_file_format, p_file_format, r_file_format,
                            f1_file_format, zeros_file_name, suffix, axs=axs[:,1], title="GWAS", algos=algos, datasets=datasets)


    plt.tight_layout()
    plt.figtext(0.01, 0.97, "A:", weight='bold', fontsize=22)
    plt.figtext(0.01, 0.5, "B:", weight='bold', fontsize=22)
    plt.figtext(0.5, 0.97, "C:", weight='bold', fontsize=22)
    plt.figtext(0.5, 0.5, "D:", weight='bold', fontsize=22)
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "figure_14.png"))




