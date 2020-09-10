import matplotlib
matplotlib.use('Agg')

import sys
sys.path.insert(0, '../')

import argparse
import seaborn as sns
sns.set(color_codes=True)
from infra import *
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np
import scipy
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
from evaluation.modules_report import modules_ehr_for_solution
import matplotlib.gridspec as gridspec
MAX_N_MODULES_TH=20
SCATTER_FACTOR=40

def calc_emp_pval(cur_rv, cur_dist):
    pos = np.size(cur_dist) - np.searchsorted(np.sort(cur_dist), cur_rv, side='left')


    return pos / float(np.size(cur_dist))


def plot_modules_ehr_summary(prefix,datasets,algos):


    df_statistics = pd.DataFrame()
    for cur_ds in datasets:
        for cur_alg in algos:
            try:
                print("{} {}".format(cur_alg, cur_ds))
                res=modules_ehr_for_solution(cur_alg, cur_ds, prefix=prefix, terms_file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR, "oob", "emp_diff_modules_{}_{}_passed_oob.tsv"),
                                  modules_file_name=os.path.join(constants.TRUE_SOLUTIONS_DIR, "{}_{}/report/modules_summary.tsv"), emp_ratio_th=0.5)
            except Exception, e:
                print "error: "+str(e)
                continue

            if res is None:
                continue
            tps, fps, sig_hg_genes, sig_emp_genes, statistics, full_data = res

            statistics["algo"]=cur_alg
            statistics["dataset"]=cur_ds
            statistics["id"]="{}_{}".format(cur_alg, cur_ds)

            df_statistics=df_statistics.append(statistics, ignore_index=True)

    df_statistics=df_statistics.set_index("id")
    df_statistics.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "mehr_cache_files", "modules_statistics_{}.tsv".format(prefix)), sep='\t')

    df_statistics = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "mehr_cache_files", "modules_statistics_{}.tsv".format(prefix)), sep='\t',
                                index_col=0)

    figure = plt.figure(constrained_layout=True, figsize=(15,24))
    spec2 = gridspec.GridSpec(ncols=15, nrows=24, figure=figure)
    h_ax_1 = figure.add_subplot(spec2[0:2, 0:9])
    h_ax_1.axis('off')
    ax_leg = figure.add_subplot(spec2[0:2, 9:14])
    ax_leg.axis('off')
    ax_cb = figure.add_subplot(spec2[0:2, 14:])
    ax_main = figure.add_subplot(spec2[2:, :])
    ax_main.set_facecolor('#ffffff')
    ax_main.grid(color='#cccccc')

    default_font_size=16

    y_axis = -1
    ids = []
    xs = []
    ys = []
    cs = []
    ss = []

    y_tick_labels=[]
    for i_a, cur_alg in enumerate(set(constants.ALGOS).intersection(algos)):
        for i_d, cur_ds in enumerate(datasets):

            ax_main.plot([-0.5, 13.5], [(len(datasets))*i_a -0.3, (len(datasets))*i_a-0.3], linestyle=(0, (5, 10)), linewidth  =1, c='black')
            id = "{}_{}".format(cur_alg, cur_ds)
            if id not in list(df_statistics.index): continue
            y_axis += 1
            y_tick_labels.append("{} / {}".format(constants.ALGOS_ACRONYM[cur_alg],cur_ds))
            ids.append(id)
            cur_series = df_statistics.loc[id, :]
            l_modules = []
            for cur_modules in range(int(cur_series["n_modules"])):
                l_modules.append((cur_series['module_{}_emp_ratio'.format(cur_modules)],
                                  cur_series['module_{}_tp'.format(cur_modules)],
                                  cur_series['module_{}_fp'.format(cur_modules)],
                                  cur_series['module_{}_total'.format(cur_modules)],
                                  cur_series['module_{}_size'.format(cur_modules)]))
            l_modules = sorted(l_modules, key=lambda x: (x[0], x[3]), reverse=True)
            for cur in range(min(10, len(l_modules))):

                if l_modules[cur][4] < 1: continue
                xs.append(cur + 1)
                ys.append(y_axis)
                cs.append(l_modules[cur][0])
                ss.append(np.max([np.log2(l_modules[cur][4]), 1]) * SCATTER_FACTOR)

            txt = "{}/{}".format(int(cur_series['emp_corrected']), int(cur_series['hg_corrected']))
            ax_main.scatter(13, y_axis, s=0)
            ax_main.annotate(txt, (10.3, y_axis), color='green', size=default_font_size)

    cs = cs
    sc = ax_main.scatter(xs, ys, s=ss, c=cs, cmap='coolwarm')
    ax_main.set_yticks(np.arange(len(y_tick_labels)))
    ax_main.set_yticklabels(tuple(y_tick_labels), size=default_font_size)
    red_list=["jAM_greedy_TNFa_2", "netbox_TNFa_2", "jAM_greedy_Type_2_Diabetes.G50", "netbox_Crohns_Disease.G50"]

    for a in ax_main.get_yticklabels():
        a.set_color("red" if a._text in red_list else "black")

    ax_main.set_xticks(np.arange(MAX_N_MODULES_TH/2+1))
    ax_main.set_xticklabels(np.append("",np.arange(MAX_N_MODULES_TH/2+1)[1:]), size=default_font_size)

    ax_main.set_xlabel("module index", size=default_font_size)
    ax_main.set_ylabel("solution\n(alg-dataset combination)", size=default_font_size)
    ax_main.set_title("mEHR", size=25)

    # produce a legend with a cross section of sizes from the scatter
    l1 = plt.scatter([], [], s=np.max([np.log2(10), 1]) * SCATTER_FACTOR, edgecolors='none', c='gray')
    l2 = plt.scatter([], [], s=np.max([np.log2(50), 1]) * SCATTER_FACTOR, edgecolors='none', c='gray')
    l3 = plt.scatter([], [], s=np.max([np.log2(100), 1]) * SCATTER_FACTOR, edgecolors='none', c='gray')
    l4 = plt.scatter([], [], s=np.max([np.log2(200), 1]) * SCATTER_FACTOR, edgecolors='none', c='gray')

    labels = ["10", "50", "100", "200"]



    plt.colorbar(mappable=sc , cax=ax_cb)
    ax_cb.tick_params(labelsize=default_font_size)
    ax_cb.set_ylabel("mEHR", size=25)

    summary_m_ehr_mean = pd.DataFrame()
    summary_m_ehr_std = pd.DataFrame()
    summary_sig = pd.DataFrame()
    n_modules_fraction = []
    for n_modules_th in range(1, MAX_N_MODULES_TH + 1):

        l_n_modules = []
        EHRs = pd.DataFrame()
        for cur_ds in datasets:
            for cur_alg in algos:
                id = "{}_{}".format(cur_alg, cur_ds)
                y_axis += 1
                if id not in list(df_statistics.index): continue
                ids.append(id)
                cur_series = df_statistics.loc[id, :]
                l_modules = []
                for cur_modules in range(int(cur_series["n_modules"])):
                    l_modules.append((cur_series['module_{}_emp_ratio'.format(cur_modules)],
                                      cur_series['module_{}_tp'.format(cur_modules)],
                                      cur_series['module_{}_fp'.format(cur_modules)],
                                      cur_series['module_{}_total'.format(cur_modules)],
                                      cur_series['module_{}_size'.format(cur_modules)]))
                l_modules = sorted(l_modules, key=lambda x: (x[0], x[3]), reverse=True)
                l_n_modules.append(len(l_modules))
                for cur in range(min(n_modules_th, len(l_modules))):

                    if l_modules[cur][4] < 1: continue
                    cs.append(l_modules[cur][0])
                    EHRs.loc["{}_{}_{}".format(cur_ds, cur_alg, cur), "ehr"] = l_modules[cur][0]
                    EHRs.loc["{}_{}_{}".format(cur_ds, cur_alg, cur), "algo"] = cur_alg
                    EHRs.loc["{}_{}_{}".format(cur_ds, cur_alg, cur), "dataset"] = cur_ds
                    EHRs.loc["{}_{}_{}".format(cur_ds, cur_alg, cur), "rank"] = cur+1

        pvals = []
        for cur_ds in datasets:
            for cur_alg in algos:
                try:
                    pvals.append(scipy.stats.mannwhitneyu(
                    EHRs[(EHRs["dataset"] == cur_ds) & (EHRs["algo"] == cur_alg)]["ehr"].values,
                    EHRs[(EHRs["dataset"] == cur_ds) & (EHRs["algo"] != cur_alg)]["ehr"].values)[1])
                except Exception,e:
                    print e
                    pvals.append(1)

        # fdrs = fdrcorrection0(pvals, is_sorted=False)[1]
        # i = 0
        # sig = pd.DataFrame()
        #

        summary_m_ehr_mean[n_modules_th] = EHRs.groupby(by=['dataset','algo'])['ehr'].mean().groupby(by=['algo']).mean()
        # summary_m_ehr_std[n_modules_th] = EHRs.groupby(by=['dataset','algo'])['ehr'].mean().groupby(by=['algo']).std()

    summary_m_ehr_mean = summary_m_ehr_mean[np.sort(summary_m_ehr_mean.columns.values)]
    # summary_m_ehr_std = summary_m_ehr_std[np.sort(summary_m_ehr_std.columns.values)]

    figure_average, ax_average = plt.subplots(figsize=(10, 12))
    summary_m_ehr_mean.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation", "mEHR_mean_{}.tsv".format(prefix)), sep='\t')
    for i, cur_row in summary_m_ehr_mean.loc[constants.ALGOS].iterrows():

        ax_average.scatter(np.arange(1, MAX_N_MODULES_TH + 1), cur_row,c=constants.COLORDICT[i])



    for i in np.arange(MAX_N_MODULES_TH/2):
        df_mEHR = pd.DataFrame(index=constants.ALGOS, columns=datasets)
        for j, cur_row in summary_m_ehr_mean.loc[constants.ALGOS].iterrows():
            df_mEHR.loc[j,:]=cur_row.iloc[i]
        df_mEHR.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "mehr_cache_files","summary_mEHR_mean_{}_{}.tsv".format(i+1, prefix)), sep='\t')

        # df_mEHR = pd.DataFrame(index=constants.ALGOS, columns=datasets)
        # for j, cur_row in summary_m_ehr_std.loc[constants.ALGOS].iterrows():
        #     df_mEHR.loc[j,:]=cur_row.iloc[i]
        # df_mEHR.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"summary_mEHR_std_{}_{}.tsv".format(i+1, prefix)), sep='\t')


    ax_average.set_xlabel("# top ranked modules", fontsize=25)
    ax_average.set_ylabel("average mEHR", fontsize=25)
    ax_average.set_facecolor('#ffffff')
    ax_average.grid(color='gray')
    patches=[Line2D([0], [0], marker='o', color='gray', label=constants.ALGOS_ACRONYM[a], markersize=12, markerfacecolor=constants.COLORDICT[a], alpha=0.7) for a in algos]
    ax_average.legend(handles=patches, fontsize=22, loc=(0,1.1), ncol=3, facecolor='#ffffff')
    ax_average.set_title("mEHR results per module", size=25)
    figure.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "figure_5_A_{}.png".format(prefix)))
    figure_average.tight_layout()
    figure_average.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "figure_5_B_{}.png".format(prefix)))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="Breast_Cancer.G50,Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50")  # TNFa_2,HC12,SHERA,SHEZH_1,ROR_1,ERS_1,IEM Breast_Cancer.G50,Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50
    parser.add_argument('--prefix', dest='prefix', default="PASCAL_SUM") # PASCAL_SUM   GE
    parser.add_argument('--base_folder_format', dest='base_folder_format', default=os.path.join(constants.OUTPUT_GLOBAL_DIR, "oob"))
    parser.add_argument('--terms_file_name_format', dest='terms_file_name_format', default="emp_diff_modules_{}_{}_passed_oob.tsv")
    parser.add_argument('--algos', dest='algos',
                        default="jactivemodules_greedy,jactivemodules_sa,netbox,bionet,dcem")  # ,keypathwayminer_INES_GREEDY,hotnet2,my_netbox_td

    args = parser.parse_args()

    datasets = args.datasets.split(",")
    algos = args.algos.split(",")
    prefix = args.prefix
    base_folder_format = args.base_folder_format
    terms_file_name_format = args.terms_file_name_format

    datasets = ["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"]
    algos = ["netbox3", "DOMINO3"] # ["DOMINO", "DOMINO2", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    prefix = "GE"
    plot_modules_ehr_summary(prefix, datasets, algos)

    datasets = ["brca", "crh", "scz", "tri", "t2d", "cad", "bmd", "hgt", "amd", "af"]
    algos = ["netbox3", "DOMINO3"] # ["DOMINO", "DOMINO2", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"]
    prefix = "PASCAL_SUM"
    plot_modules_ehr_summary(prefix, datasets, algos)


