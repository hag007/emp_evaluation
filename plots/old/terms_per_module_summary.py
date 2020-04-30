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
    if output_terms.shape[0]==0:
        n_modules=0
    else:
        n_modules=output_terms.iloc[0]["hg_pval"].count(',')+1
    statistics["n_modules"]=n_modules
    for a in range(n_modules):
        tps[a]=[]
        fps[a]=[]

    for go_id , cur in output_terms.iterrows():
        hgs=list(np.array(cur["hg_pval"][1:-1].split(", "),dtype=np.float32))
        emps=list(np.array(cur["emp_pval"][1:-1].split(", "),dtype=np.float32))
        for i, v in enumerate(zip(hgs,emps)):
            hg,emp=v
            if hg >= HG_CUTOFF and emp <= EMP_CUTOFF:
                tps[i].append("{}: {}".format(go_id, cur["GO name"]))
            elif hg >= HG_CUTOFF and emp > EMP_CUTOFF:
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



    # terms_limit = 0
    # df_rank_matrix = pd.DataFrame()
    # df_matrix = pd.DataFrame()
    # df_count_matrix = pd.DataFrame()
    # df_all = pd.DataFrame()
    # df_statistics = pd.DataFrame()
    # for cur_ds in datasets:
    #     df_ds=pd.DataFrame()
    #     for cur_alg in algos:
    #         terms_file_name=os.path.join(base_folder, terms_file_name_format.format(cur_ds,cur_alg))
    #         modules_file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR, "{}_{}".format(prefix,cur_ds), cur_alg, "modules_summary.tsv")
    #         res=calc_modules_ehr(terms_file_name, modules_file_name, cur_alg, cur_ds)
    #         if res is None:
    #             continue
    #         tps, fps, sig_hg_genes, sig_emp_genes, statistics = res
    #         statistics["algo"]=cur_alg
    #         statistics["dataset"]=cur_ds
    #         statistics["id"]="{}_{}".format(cur_alg, cur_ds)
    #
    #         df_statistics=df_statistics.append(statistics, ignore_index=True)
    #
    # df_statistics=df_statistics.set_index("id")
    # df_statistics.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"modules_statistics.tsv"), sep='\t')

    df_statistics = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "modules_statistics.tsv"), sep='\t',
                                index_col=0)

    figure = plt.figure(figsize=(35, 35))
    # gs = figure.add_gridspec(3, 3, figsize=(30, 30))
    subplots=[]
    subplots.append(plt.subplot2grid((30, 30), (0, 0), rowspan=30, colspan=17))
    subplots.append(plt.subplot2grid((30, 30), (1, 20), rowspan=9, colspan=9))
    subplots.append(plt.subplot2grid((30, 30), (11, 20), rowspan=9, colspan=9))
    subplots.append(plt.subplot2grid((30, 30), (21, 20), rowspan=9, colspan=9))

    ax=subplots[0]

    # ax.figure.set_size_inches(20,30)
    ax.set_facecolor('#fffde3')
    ax.grid(color='gray')

    y_axis = -1
    ids = []
    xs = []
    ys = []
    cs = []
    ss = []

    algos=np.sort(np.unique(df_statistics["algo"]))
    for i_a, cur_alg in enumerate(algos):
        for i_d, cur_ds in enumerate(datasets):

            ax.plot([-0.5, 13.5], [(len(datasets))*i_a -0.5, (len(datasets))*i_a-0.5], linestyle=(0, (5, 10)), linewidth  =4, c='black')



            y_axis += 1
            id = "{}_{}".format(cur_alg, cur_ds)
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
                ss.append(np.max([np.log2(l_modules[cur][4]), 1]) * 75)

            txt = "EHR={}/{}".format(int(cur_series['emp_corrected']), int(cur_series['hg_corrected']))
            ax.scatter(13, y_axis, s=0)
            ax.annotate(txt, (10.3, y_axis), color='green', size=20)

    cs = cs
    sc = ax.scatter(xs, ys, s=ss, c=cs, cmap='bwr')

    ids = ids
    ax.set_yticks(np.arange(len(ids)))
    ax.set_yticklabels(tuple(ids), size=25)

    ax.set_xticks(np.arange(11))
    ax.set_xticklabels(np.append("", (np.arange(12) - 1)[1:]), size=25)

    ax.set_xlabel("module index", size=25)
    ax.set_ylabel("solution\n(algo-dataset combination)", size=25)
    cax = figure.add_axes([0.67, 0.01,0.01, 0.97])
    # aspect = 30
    # pad_fraction = 0.5
    # from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
    # divider = make_axes_locatable(ax_)
    # width = axes_size.AxesY(ax_, aspect=1. / aspect)
    # pad = axes_size.Fraction(pad_fraction, width)
    # cax = divider.append_axes("right", size=0.2, pad=0.4)
    plt.colorbar(mappable=sc, cax=cax)
    cax.tick_params(labelsize=25)
    cax.set_ylabel("EHR", size=25)

    plt.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "modules_statistics_plot.png"))

    max_n_modules_th =2
    summary_m_ehr = pd.DataFrame()
    summary_sig = pd.DataFrame()
    n_modules_fraction = []
    for n_modules_th in range(1, max_n_modules_th + 1):

        l_n_modules = []
        EHRs = pd.DataFrame()
        for cur_ds in datasets:
            df_ds = pd.DataFrame()
            for cur_alg in algos:
                y_axis += 1
                id = "{}_{}".format(cur_alg, cur_ds)
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

        summary_m_ehr[n_modules_th] = EHRs.groupby(by=['algo'])['ehr'].mean()
        print summary_m_ehr[n_modules_th]

        summary_sig[n_modules_th] = sig.groupby(by=['algo'])['qval_bool'].sum()
        print summary_sig[n_modules_th]

    summary_m_ehr = summary_m_ehr[np.sort(summary_m_ehr.columns.values)]
    summary_sig = summary_sig[np.sort(summary_m_ehr.columns.values)]



    # plt.clf()
    for i, cur_row in summary_m_ehr.iterrows():
        subplots[1].plot(np.arange(1, max_n_modules_th + 1), cur_row)

    subplots[1].set_xlabel("# top ranked EHR modules\n(modules head threshold)", fontsize=25)
    subplots[1].set_ylabel("average EHR", fontsize=25)
    # subplots[1].set_title("average EHR as function of modules' head threshold")
    subplots[1].set_facecolor('#fffde3')
    subplots[1].grid(color='gray')
    subplots[1].legend(fontsize=25)
    # subplots[0].savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,"summary_m_ehr.png"))

    # plt.clf()
    for i, cur_row in summary_sig.iterrows():
        subplots[2].plot(np.arange(1, max_n_modules_th + 1), cur_row)
    subplots[2].set_xlabel("# top ranked EHR modules\n(modules head threshold)", fontsize=25)
    subplots[2].set_ylabel("# significant datasets", fontsize=25)
    # subplots[2].set_title("# significant datasets as function of modules' head threshold")
    subplots[2].set_facecolor('#fffde3')
    subplots[2].grid(color='gray')
    subplots[2].legend(fontsize=25)
    # subplots[1].savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,"summary_sig.png"))

    # plt.clf()
    subplots[3].plot(np.arange(1, max_n_modules_th + 1), n_modules_fraction)
    subplots[3].set_xlabel("# top ranked EHR modules\n(modules head threshold)", fontsize=25)
    subplots[3].set_ylabel("fraction of solutions where\n# modules >= modules head threshold", fontsize=25)
    # subplots[3].set_title(
    #     "fraction of (modules >= modules head threshold) solutions\nas function of modules' head threshold")
    subplots[3].set_facecolor('#fffde3')
    subplots[3].grid(color='gray')
    subplots[3].legend(fontsize=25)

    plt.figtext(0.01, 0.99, "A:", weight='bold', size=25)
    plt.figtext(0.69, 0.99, "B:", weight='bold', size=25)
    plt.figtext(0.69, 0.33, "D:", weight='bold', size=25)
    plt.figtext(0.69, 0.66, "C:", weight='bold', size=25)


    plt.tight_layout()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "figure_5_{}.png".format(prefix)))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="Breast_Cancer.G50,Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50")  # TNFa_2,HC12,SHERA,SHEZH_1,ROR_1,ERS_1,IEM Breast_Cancer.G50,Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50
    parser.add_argument('--prefix', dest='prefix', default="PASCAL_SUM") # PASCAL_SUM   GE
    parser.add_argument('--base_folder_format', dest='base_folder_format', default="/home/hag007/Desktop/aggregate{}_report/oob")
    parser.add_argument('--terms_file_name_format', dest='terms_file_name_format', default="emp_diff_modules_{}_{}_passed_oob.tsv")
    parser.add_argument('--algos', dest='algos',
                        default="jactivemodules_greedy,jactivemodules_sa,netbox,bionet,dcem")  # ,keypathwayminer_INES_GREEDY,hotnet2,my_netbox_td

    args = parser.parse_args()

    datasets = args.datasets.split(",")
    algos = args.algos.split(",")
    prefix = args.prefix
    base_folder_format = args.base_folder_format
    terms_file_name_format = args.terms_file_name_format




    # prefix = "GE"
    # algos = ["dcem", "dcem3", "dcem4", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox", "keypathwayminer_INES_GREEDY",
    #          "hotnet2", "my_netbox_td"]
    # datasets = ["TNFa_2", "HC12", "SHERA", "SHEZH_1", "ROR_1", "ERS_1", "IEM"]
    # omic_type = ""
    # plot_modules_ehr_summary(prefix, datasets, algos, base_folder_format.format(omic_type), terms_file_name_format)

    prefix = "PASCAL_SUM"
    algos = ["dcem", "dcem3", "dcem4", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox", "keypathwayminer_INES_GREEDY",
             "hotnet2", "my_netbox_td"]
    datasets=["Breast_Cancer.G50", "Crohns_Disease.G50", "Schizophrenia.G50", "Triglycerides.G50", "Type_2_Diabetes.G50"]
    omic_type = "_gwas"
    plot_modules_ehr_summary(prefix, datasets, algos, base_folder_format.format(omic_type), terms_file_name_format)



