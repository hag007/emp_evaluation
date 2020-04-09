import pandas as pd

from infra import *

from rpy2.robjects import pandas2ri
pandas2ri.activate()
import constants

from matplotlib.lines import Line2D


def plot_itegrated_recovery(ss_ratios, base_folder, average_file_format, auc_file_format, p_file_format, r_file_format, f1_file_format, empty_file_format, zero_file_name, suffix, axs, title="", algos=constants.ALGOS_ACRONYM.keys(), datasets=[]):

    axs[0].set_facecolor('#ffffff')
    axs[0].grid(color='gray')
    axs[1].set_facecolor('#ffffff')
    axs[1].grid(color='gray')
      #

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
        n_iterations = 100 # len(
            # pd.read_csv(os.path.join(base_folder, average_file_format.format(suffix, cur_ss_ratio)), sep='\t',
            #             index_col=0)["precisions"][0][1:-1].split(", "))
        print "n_teration for ss_ratio={}: {}".format(cur_ss_ratio, n_iterations)
        df_cur_pr = pd.read_csv(os.path.join(base_folder, auc_file_format.format(suffix, cur_ss_ratio)), sep='\t',
                                index_col=0).loc[:,datasets]#[~np.logical_or(df_zeros == 0, np.isnan(df_zeros))]
        df_cur_pr=df_cur_pr.mean(axis=1)
        df_cur_p = pd.read_csv(os.path.join(base_folder, p_file_format.format(suffix, cur_ss_ratio)), sep='\t',
                               index_col=0).loc[:,datasets]#[~np.logical_or(df_zeros == 0, np.isnan(df_zeros))]
        df_cur_p=df_cur_p.mean(axis=1)
        df_cur_r = pd.read_csv(os.path.join(base_folder, r_file_format.format(suffix, cur_ss_ratio)), sep='\t',
                               index_col=0).loc[:,datasets]#[~np.logical_or(df_zeros == 0, np.isnan(df_zeros))]
        df_cur_r=df_cur_r.mean(axis=1)
        df_cur_f1 = pd.read_csv(os.path.join(base_folder, f1_file_format.format(suffix, cur_ss_ratio)), sep='\t',
                                index_col=0).loc[:,datasets]#[~np.logical_or(df_zeros == 0, np.isnan(df_zeros))]
        df_cur_f1=df_cur_f1.mean(axis=1)
        df_cur_empty = pd.read_csv(os.path.join(base_folder, empty_file_format.format(suffix, cur_ss_ratio)), sep='\t',
                                   index_col=0).loc[:,datasets]#[~np.logical_or(df_zeros == 0, np.isnan(df_zeros))].mean(axis=1) / float(n_iterations)
        df_cur_empty=df_cur_empty.mean(axis=1)


        prs = pd.concat([prs, df_cur_pr], axis=1)
        ps = pd.concat([ps, df_cur_p], axis=1)
        rs = pd.concat([rs, df_cur_r], axis=1)
        f1s = pd.concat([f1s, df_cur_f1], axis=1)
        empties = pd.concat([empties, df_cur_empty], axis=1)

    prs = prs.loc[np.sort(ps.index.values)]
    ps = ps.loc[np.sort(ps.index.values)]
    rs = rs.loc[np.sort(ps.index.values)]
    f1s = f1s.loc[np.sort(ps.index.values)]
    empties = empties.loc[np.sort(ps.index.values)]


    patches_0 = [Line2D([0], [0], marker='o', color='gray', label=constants.ALGOS_ACRONYM[a], markersize=12, markerfacecolor=constants.COLORDICT[a], alpha=0.7) for a in algos] # + \
                # [Line2D([0], [0], linestyle='-', color='black', label='AUPR', markersize=12, markerfacecolor='gray',
               # alpha=0.7)] + \
               #  [Line2D([0], [0], linestyle=(0, (5, 10)), color='black', label='fraction of\nnon-empty solutions', markersize=12,
               # markerfacecolor='gray', alpha=0.7)]

    patches_1 = [Line2D([0], [0], marker='o', color='gray', label=constants.ALGOS_ACRONYM[a], markersize=12, markerfacecolor=constants.COLORDICT[a], alpha=0.7) for a in algos]
        # + [
        # Line2D([0], [0], linestyle='-', linewidth=3.0, color='black', label='f1', markersize=12, markerfacecolor='gray',
        #        alpha=0.7)]+ \
        #         [Line2D([0], [0], linestyle='dashed', color='black', label='precision', markersize=12, markerfacecolor='gray',
        #        alpha=0.7)] + \
        #         [Line2D([0], [0], linestyle='dotted', color='black', label='recall', markersize=12,
        #                           markerfacecolor='gray', alpha=0.7)]

    ss_ratios = 1 - np.array(ss_ratios)

    prs.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"prs_{}.tsv".format(title)), sep='\t')
    empties.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"empties_{}.tsv".format(title)), sep='\t')
    f1s.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"f1s_{}.tsv".format(title)), sep='\t')


    for cur in set(prs.index).intersection(constants.ALGOS):
        # if prs.index[cur_i] == "my_netbox_td": continue
        axs[0].plot([a for a,b in zip(ss_ratios,prs.loc[cur]) if not np.isnan(b)], [a for a in prs.loc[cur] if not np.isnan(a) ], c=constants.COLORDICT[cur])
        # axs[0].plot(ss_ratios, empties.loc[cur], c=constants.COLORDICT[cur], linestyle=(0, (5, 10)), linewidth=3.0)
        # plt.plot(ss_ratios, np.multiply(empties.iloc[cur_i] / 100.0, prs.iloc[cur_i]) , c=colorlist[cur_i], linestyle='dashed')
        # axs[1].plot(ss_ratios, ps.iloc[cur_i], c=colorlist[cur_i], linestyle='dashed')
        # axs[1].plot(ss_ratios, rs.iloc[cur_i], c=colorlist[cur_i], linestyle='dotted')
        axs[1].plot(ss_ratios, f1s.loc[cur], c=constants.COLORDICT[cur], linewidth=2.0)

    axs[0].set_xlabel("subsample fraction", fontsize=22)
    axs[0].set_ylabel("AUPR", fontsize=22)
    # ax2 = axs[0].twinx()
    # ax2.set_ylabel("Fration of non-empty iterations", fontsize=22)  #
    axs[0].set_title(title, fontdict={"size":22})

    axs[1].set_xlabel("subsample fraction", fontsize=22)
    axs[1].set_ylabel("F1", fontsize=22)
    axs[1].set_title(title, fontdict={"size":22})

    axs[0].legend(handles=patches_0, loc=(0,1.1), ncol=3, fontsize=20, facecolor='#ffffff')
    axs[1].legend(handles=patches_1, loc=(0,1.1), ncol=3, fontsize=20, facecolor='#ffffff')



if __name__=="__main__":

    base_folder="/media/hag007/Data/bnet/output/emp_fdr/MAX"
    auc_file_format = "pr_auc_recovery_summary_{}_{}.tsv"
    p_file_format = "recovery_results_{}_{}_matrix_p.tsv"
    r_file_format = "recovery_results_{}_{}_matrix_r.tsv"
    f1_file_format = "recovery_results_{}_{}_matrix_f1.tsv"
    empty_file_format = "recovery_results_{}_{}_matrix_empty.tsv"
    average_file_format="recovery_results_{}_{}.tsv"
    zeros_file_format = os.path.join('/home/hag007/Desktop/aggregate{}_report/venn',"count_matrix.tsv")
    ss_ratios = [0.1, 0.2, 0.3, 0.4]

    fig,axs=plt.subplots(2,2,figsize=(18,16))
    suffix = "GE_100"
    omic_type=""
    zeros_file_name = zeros_file_format.format(omic_type)
    algos = ["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox", "keypathwayminer_INES_GREEDY","domino_original"]
    datasets = ["TNFa_2", "HC12", "SHERA", "SHEZH_1", "ROR_1", "ERS_1", "IEM" , "APO", "CBX", "IFT"]
    plot_itegrated_recovery(ss_ratios, base_folder, average_file_format, auc_file_format, p_file_format, r_file_format, f1_file_format, empty_file_format, zeros_file_name, suffix, axs=axs[:,0], title="GE", algos=algos, datasets=datasets)
    suffix = "PASCAL_SUM_100"
    omic_type="_gwas"
    zeros_file_name = zeros_file_format.format(omic_type)
    algos = ["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox","domino_original"]
    datasets=["Breast_Cancer.G50", "Crohns_Disease.G50", "Schizophrenia.G50", "Triglycerides.G50", "Type_2_Diabetes.G50" ,"Coronary_Artery_Disease.G50"  , "Bone_Mineral_Density.G50", "Height1.G50", "Age_Related_Macular_Degeneration.G50", "Atrial_Fibrillation.G50"]
    plot_itegrated_recovery(ss_ratios, base_folder, average_file_format, auc_file_format, p_file_format, r_file_format,
                            f1_file_format, empty_file_format, zeros_file_name, suffix, axs=axs[:,1], title="GWAS", algos=algos, datasets=datasets)


    plt.tight_layout()
    plt.figtext(0.01, 0.97, "A:", weight='bold', fontsize=22)
    plt.figtext(0.01, 0.5, "B:", weight='bold', fontsize=22)
    plt.figtext(0.5, 0.97, "C:", weight='bold', fontsize=22)
    plt.figtext(0.5, 0.5, "D:", weight='bold', fontsize=22)
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_14.png"))




