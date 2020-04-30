import json
import matplotlib
# from matplotlib_venn import venn2
matplotlib.use('Agg')

from pandas._libs.parsers import k
import sys
sys.path.insert(0, '../')

import argparse
import seaborn as sns
sns.set(color_codes=True)
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from constants import *
from infra import *
from utils.param_builder import build_gdc_params
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as ml_colors
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np
import networkx as nx
import shutil
from pandas.errors import EmptyDataError
from utils.permute_network import EdgeSwapGraph
import scipy
from scipy.optimize import least_squares
from runners.FDR_runner import run_FDR
import random
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
import multiprocessing
from utils.daemon_multiprocessing import func_star

from matplotlib.patches import Circle

MAX_N_MODULES_TH=20
SCATTER_FACTOR=40

def calc_emp_pval(cur_rv, cur_dist):
    pos = np.size(cur_dist) - np.searchsorted(np.sort(cur_dist), cur_rv, side='left')


    return pos / float(np.size(cur_dist))



def calc_modules_ehr(terms_file_name, modules_file_name, algo_sample=None, dataset_sample=None, emp_ratio_th=0.5):
    statistics = {}
    full_data = pd.DataFrame()
    try:
        output_terms = pd.read_csv(terms_file_name, sep='\t', index_col=0).dropna()
        output_modules = pd.read_csv(modules_file_name, sep='\t', index_col=0).dropna()
    except IOError,e:
        print e
        return None



    output_terms = output_terms.rename(columns={"filtered_pval": "hg_pval_max"})
    filtered_genes=output_terms.loc[np.logical_and.reduce([output_terms["n_genes"].values > 5, output_terms["n_genes"].values < 500]), ["GO name","hg_pval_max", "emp_pval_max", "passed_oob_permutation_test"]]


    print "total n_genes with pval:{}/{}".format(np.size(filtered_genes["hg_pval_max"].values), 7435)

    sorted_genes_hg = filtered_genes.sort_values(by=['hg_pval_max'], ascending=False)
    sig_genes_hg_pval=np.append(sorted_genes_hg["hg_pval_max"].values,np.zeros(7435-np.size(sorted_genes_hg["hg_pval_max"].values)))
    sig_genes_hg_pval = [10**(-x) for x in sig_genes_hg_pval]
    fdr_results = fdrcorrection0(sig_genes_hg_pval, alpha=0.05, method='indep', is_sorted=False)
    n_hg_true = len([cur for cur in fdr_results[0] if cur])
    sig_hg_genes= sorted_genes_hg.iloc[:n_hg_true, :] if n_hg_true > 0 else 0
    if type(sig_hg_genes)==int:
        HG_CUTOFF = 0
        statistics["hg_cutoff"] = 0
        statistics["hg_corrected"] = 0
        print "HG cutoff: {}, n={}".format(HG_CUTOFF, 0)
    else:
        HG_CUTOFF = sig_hg_genes.iloc[- 1]["hg_pval_max"]
        statistics["hg_cutoff"]=round(HG_CUTOFF,2)
        statistics["hg_corrected"]=len(sig_hg_genes.index)
        print "HG cutoff: {}, n={}".format(HG_CUTOFF, len(sig_hg_genes.index))

    sorted_genes_emp = filtered_genes.sort_values(by=['emp_pval_max'])
    sorted_genes_emp.loc[sorted_genes_emp['emp_pval_max']==0,'emp_pval_max']=1.0/5000
    sig_genes_emp_pval = sorted_genes_emp["emp_pval_max"].values
    fdr_results = fdrcorrection0(sig_genes_emp_pval, alpha=0.05, method='indep', is_sorted=False)
    n_emp_true = len([cur for cur in fdr_results[0] if cur])
    sig_emp_genes = sorted_genes_emp.iloc[:n_emp_true, :]
    EMP_CUTOFF = sig_emp_genes.iloc[- 1]["emp_pval_max"] if n_emp_true > 0 else 0
    statistics["emp_cutoff"]=EMP_CUTOFF
    statistics["emp_corrected"]=len(sig_emp_genes.index)
    print "EMP cutoff: {}, n={}".format(EMP_CUTOFF, len(sig_emp_genes.index))

    tps={}
    fps = {}
    hg_separated={}
    if output_terms.shape[0]==0:
        n_modules=0
    else:
        n_modules=output_terms.iloc[0]["hg_pval"].count(',')+1
    statistics["n_modules"]=n_modules
    for a in range(n_modules):
        tps[a]=[]
        fps[a]=[]
        # hg_file_name = os.path.join(constants.OUTPUT_GLOBAL_DIR, "{}_{}".format(prefix,dataset_sample), algo_sample,
        #                                     "module_{}_separated_modules_hg_samples.tsv".format(a))
        # hg_separated[a]=pd.read_csv(hg_file_name,sep='\t', index_col=0)
    for go_id , cur in output_terms.iterrows():
        hgs=list(np.array(cur["hg_pval"][1:-1].split(", "),dtype=np.float32))
        emps=list(np.array(cur["emp_pval"][1:-1].split(", "),dtype=np.float32))
        for i, v in enumerate(zip(hgs,emps)):
            hg,emp=v
            if emp <= EMP_CUTOFF and hg >= HG_CUTOFF : # go_id in hg_separated[i].index
                tps[i].append("{}: {}".format(go_id, cur["GO name"]))
            elif emp > EMP_CUTOFF and hg >= HG_CUTOFF : # go_id in hg_separated[i].index
                fps[i].append("{}: {}".format(go_id, cur["GO name"]))

    real_modules_counter=0
    for a in range(n_modules):
        n_hg_terms=len(tps[a]) + len(fps[a])
        # print "module_{}: {} ({}/{})".format(a, float(len(tps[a]))/n_hg_terms, len(tps[a]), n_hg_terms)
        statistics["module_{}_emp_ratio".format(a)]=round(float(len(tps[a]))/max(len(tps[a])+len(fps[a]),1),2)
        statistics["module_{}_tp".format(a)] = len(tps[a])
        statistics["module_{}_fp".format(a)] = len(fps[a])
        statistics["module_{}_total".format(a)] = len(tps[a])+len(fps[a])
        statistics["module_{}_size".format(a)] = output_modules.loc[a, '#_genes']

        if statistics["module_{}_emp_ratio".format(a)] >emp_ratio_th:
            real_modules_counter+=1
            # for cur in tps[a]:
            #     print cur
        full_data.loc["module_{}".format(a), "EHR"] = round(float(len(tps[a])) / max(len(tps[a]) + len(fps[a]), 1),2)
        full_data.loc["module_{}".format(a),"tp"]="\n".join(tps[a])
        full_data.loc["module_{}".format(a),"fp"]="\n".join(fps[a])
        full_data.sort_index(inplace=True)
    statistics["real_modules_ratio"]=round(float(real_modules_counter)/max(n_modules,1),2)
    # print "real_modules_counter: {}. ratio: {}".format(real_modules_counter, statistics["real_modules_ratio"])


    for cur_key in np.sort(statistics.keys()):
        print "{}: {}".format(cur_key, statistics[cur_key])

    file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR, "full_modules_report_{}_{}.tsv".format(algo_sample, dataset_sample))
    full_data.to_csv(file_name,sep='\t')

    with file(file_name,'a+') as f:
        f.write("\nsolution statistics: \n")
        for k in np.sort(statistics.keys()):
            if not k.startswith("modules_"):
                f.write("{}\t{}\n".format(k,statistics[k]))

        for k in np.sort(statistics.keys()):
            if k.startswith("modules_"):
                f.write("{}\t{}\n".format(k,statistics[k]))


    return tps, fps, sig_hg_genes, sig_emp_genes, statistics



def plot_modules_ehr_summary(prefix,datasets,algos, base_folder, terms_file_name_format):



    terms_limit = 0
    df_rank_matrix = pd.DataFrame()
    df_matrix = pd.DataFrame()
    df_count_matrix = pd.DataFrame()
    df_all = pd.DataFrame()
    df_statistics = pd.DataFrame()
    for cur_ds in datasets:
        df_ds=pd.DataFrame()
        for cur_alg in algos:
            terms_file_name=os.path.join(base_folder, terms_file_name_format.format(cur_ds,cur_alg))
            modules_file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR, "{}_{}".format(prefix,cur_ds), cur_alg, "modules_summary.tsv")
            res=calc_modules_ehr(terms_file_name, modules_file_name, cur_alg, cur_ds)
            if res is None:
                continue
            tps, fps, sig_hg_genes, sig_emp_genes, statistics = res
            statistics["algo"]=cur_alg
            statistics["dataset"]=cur_ds
            statistics["id"]="{}_{}".format(cur_alg, cur_ds)

            df_statistics=df_statistics.append(statistics, ignore_index=True)

    df_statistics=df_statistics.set_index("id")
    df_statistics.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"modules_statistics_{}.tsv".format(prefix)), sep='\t')

    df_statistics = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "modules_statistics_{}.tsv".format(prefix)), sep='\t',
                                index_col=0)

    figure = plt.figure(constrained_layout=True, figsize=(15,24))
    import matplotlib.gridspec as gridspec
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

            y_tick_labels.append("{}_{}".format(constants.ALGOS_ACRONYM[cur_alg],constants.DATASETS_ACRONYM[cur_ds]))
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

    # ax.tick_params(axis='y', colors=['red' if a in red_list else 'black' ])

    ax_main.set_xticks(np.arange(MAX_N_MODULES_TH/2+1))
    ax_main.set_xticklabels(np.append("",np.arange(MAX_N_MODULES_TH/2+1)[1:]), size=default_font_size)

    ax_main.set_xlabel("module index", size=default_font_size)
    ax_main.set_ylabel("solution\n(alg-dataset combination)", size=default_font_size)
    ax_main.set_title("mEHR results per module", size=25)

    # produce a legend with a cross section of sizes from the scatter
    l1 = plt.scatter([], [], s=np.max([np.log2(10), 1]) * SCATTER_FACTOR, edgecolors='none', c='gray')
    l2 = plt.scatter([], [], s=np.max([np.log2(50), 1]) * SCATTER_FACTOR, edgecolors='none', c='gray')
    l3 = plt.scatter([], [], s=np.max([np.log2(100), 1]) * SCATTER_FACTOR, edgecolors='none', c='gray')
    l4 = plt.scatter([], [], s=np.max([np.log2(200), 1]) * SCATTER_FACTOR, edgecolors='none', c='gray')

    labels = ["10", "50", "100", "200"]


    leg = ax_leg.legend([l1, l2, l3, l4], labels, ncol=2, frameon=True, fontsize=default_font_size,
                     handlelength=2, loc='lower left',
                     handletextpad=1, title='# of genes in module', scatterpoints=1, facecolor='#ffffff')
    # leg = h_ax_1.legend([l1, l2, l3, l4], labels, ncol=2, frameon=True, fontsize=default_font_size,
    #                  handlelength=2, loc=(10,10),
    #                  handletextpad=1, title='# of genes in module', scatterpoints=1, facecolor='#ffffff')


    # # cax = figure.add_axes([1.0, 0.01,0.01, 1.1])
    # # aspect = 30
    # # pad_fraction = 0.5
    # from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
    # divider = make_axes_locatable(ax)
    # # width = axes_size.AxesY(ax_, aspect=1. / aspect)
    # # pad = axes_size.Fraction(pad_fraction, width)
    # cax = divider.append_axes("right", size="3%", pad=0.3)

    plt.colorbar(mappable=sc , cax=ax_cb)
    ax_cb.tick_params(labelsize=default_font_size)
    ax_cb.set_ylabel("mEHR", size=25)

    # plt.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "modules_statistics_plot.png"))

    summary_m_ehr_mean = pd.DataFrame()
    summary_m_ehr_std = pd.DataFrame()
    summary_sig = pd.DataFrame()
    n_modules_fraction = []
    for n_modules_th in range(1, MAX_N_MODULES_TH + 1):

        l_n_modules = []
        EHRs = pd.DataFrame()
        for cur_ds in datasets:
            df_ds = pd.DataFrame()
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

        fdrs = fdrcorrection0(pvals, is_sorted=False)[1]
        i = 0
        sig = pd.DataFrame()

        for cur_ds in datasets:
            for cur_alg in algos:
                mn_alg = EHRs[(EHRs["dataset"] == cur_ds) & (EHRs["algo"] == cur_alg)]["ehr"].values.mean()
                mn_all = EHRs[(EHRs["dataset"] == cur_ds)]["ehr"].values.mean()
                print "{} > {} ({})".format(mn_alg, mn_all, mn_alg > mn_all and fdrs[i] < 0.05)

                sig.loc["{}_{}".format(cur_ds, cur_alg), "dataset"] = cur_ds
                sig.loc["{}_{}".format(cur_ds, cur_alg), "algo"] = cur_alg
                sig.loc["{}_{}".format(cur_ds, cur_alg), "pval"] = pvals[i]
                sig.loc["{}_{}".format(cur_ds, cur_alg), "qval"] = fdrs[i]
                sig.loc["{}_{}".format(cur_ds, cur_alg), "qval_bool"] = int((fdrs[i] < 0.05) & (mn_alg > mn_all))

                print "{} {} mwu pval: {} (qval={})".format(cur_ds, cur_alg, pvals[i], fdrs[i])
                i += 1

        print sig

        n_modules_fraction.append(np.sum(np.array(l_n_modules) >= n_modules_th) / float(len(np.unique(EHRs["dataset"])) * len(np.unique(EHRs["algo"]))))
        print "fraction over n_modules th:", n_modules_fraction[-1]

        summary_m_ehr_mean[n_modules_th] = EHRs.groupby(by=['algo'])['ehr'].mean()
        summary_m_ehr_std[n_modules_th] = EHRs.groupby(by=['algo'])['ehr'].std()
        print summary_m_ehr_mean[n_modules_th]

        summary_sig[n_modules_th] = sig.groupby(by=['algo'])['qval_bool'].sum()
        print summary_sig[n_modules_th]

    summary_m_ehr_mean = summary_m_ehr_mean[np.sort(summary_m_ehr_mean.columns.values)]
    summary_m_ehr_std = summary_m_ehr_std[np.sort(summary_m_ehr_std.columns.values)]
    summary_sig = summary_sig[np.sort(summary_m_ehr_mean.columns.values)]

    figure_average, ax_average = plt.subplots(figsize=(10, 12))
    summary_m_ehr_mean.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "mEHR_mean_{}.tsv".format(prefix)), sep='\t')
    summary_m_ehr_std.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "mEHR_std_{}.tsv".format(prefix)), sep='\t')
    for i, cur_row in summary_m_ehr_mean.loc[constants.ALGOS].iterrows():

        ax_average.scatter(np.arange(1, MAX_N_MODULES_TH + 1), cur_row,c=constants.COLORDICT[i])



    for i in np.arange(MAX_N_MODULES_TH/2):
        df_mEHR = pd.DataFrame(index=constants.ALGOS, columns=datasets)
        for j, cur_row in summary_m_ehr_mean.loc[constants.ALGOS].iterrows():
            df_mEHR.loc[j,:]=cur_row.iloc[i]
        df_mEHR.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"summary_mEHR_mean_{}_{}.tsv".format(i+1, prefix)), sep='\t')

        df_mEHR = pd.DataFrame(index=constants.ALGOS, columns=datasets)
        for j, cur_row in summary_m_ehr_std.loc[constants.ALGOS].iterrows():
            df_mEHR.loc[j,:]=cur_row.iloc[i]
        df_mEHR.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"summary_mEHR_std_{}_{}.tsv".format(i+1, prefix)), sep='\t')


    ax_average.set_xlabel("# top ranked EHR modules", fontsize=25)
    ax_average.set_ylabel("average mEHR", fontsize=25)
    # ax_average.set_title("average EHR as function of modules' head threshold")
    ax_average.set_facecolor('#ffffff')
    ax_average.grid(color='gray')
    patches=[Line2D([0], [0], marker='o', color='gray', label=constants.ALGOS_ACRONYM[a], markersize=12, markerfacecolor=constants.COLORDICT[a], alpha=0.7) for a in algos]
    ax_average.legend(handles=patches, fontsize=22, loc=(0,1.1), ncol=3, facecolor='#ffffff')
    ax_average.set_title("mEHR results per module", size=25)
    # subplots[0].savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,"summary_m_ehr.png"))


    # figure.text(0.02, 0.97, "A:", weight='bold', size=25)
    # figure_average.text(0.02, 0.97, "B:", weight='bold', size=18)

    # figure.tight_layout()
    figure.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_5_A_{}.png".format(prefix)))
    figure_average.tight_layout()
    figure_average.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_5_B_{}.png".format(prefix)))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="Breast_Cancer.G50,Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50")  # TNFa_2,HC12,SHERA,SHEZH_1,ROR_1,ERS_1,IEM Breast_Cancer.G50,Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50
    parser.add_argument('--prefix', dest='prefix', default="PASCAL_SUM") # PASCAL_SUM   GE
    parser.add_argument('--base_folder_format', dest='base_folder_format', default=os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr","MAX"))
    parser.add_argument('--terms_file_name_format', dest='terms_file_name_format', default="emp_diff_modules_{}_{}_passed_oob.tsv")
    parser.add_argument('--algos', dest='algos',
                        default="jactivemodules_greedy,jactivemodules_sa,netbox,bionet,dcem")  # ,keypathwayminer_INES_GREEDY,hotnet2,my_netbox_td

    args = parser.parse_args()

    datasets = args.datasets.split(",")
    algos = args.algos.split(",")
    prefix = args.prefix
    base_folder_format = args.base_folder_format
    terms_file_name_format = args.terms_file_name_format




    prefix = "GE"
    algos = ["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox", "keypathwayminer_INES_GREEDY", "domino_original"] # "dcem", "dcem2", "dcem3", "dcem4", "my_netbox_td"
    datasets = ["TNFa_2", "HC12", "SHERA", "SHEZH_1", "ROR_1", "ERS_1", "IEM", "APO", "CBX", "IFT"]
    omic_type = ""
    plot_modules_ehr_summary(prefix, datasets, algos, base_folder_format.format(omic_type), terms_file_name_format)

    prefix = "PASCAL_SUM"
    algos = ["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox", "domino_original"] # "dcem", "dcem2", "dcem3", "dcem4", "my_netbox_td"
    datasets=["Breast_Cancer.G50", "Crohns_Disease.G50", "Schizophrenia.G50", "Triglycerides.G50", "Type_2_Diabetes.G50" ,"Coronary_Artery_Disease.G50" , "Bone_Mineral_Density.G50", "Height1.G50", "Age_Related_Macular_Degeneration.G50", "Atrial_Fibrillation.G50"] #  , "Alzheimer.G50"
    omic_type = "_gwas"
    plot_modules_ehr_summary(prefix, datasets, algos, base_folder_format.format(omic_type), terms_file_name_format)


