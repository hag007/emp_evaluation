import json
import matplotlib
from matplotlib_venn import venn2
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

from scipy.stats import hypergeom

from utils.revigo import sumbit_to_revigo

def load_gs():
    bg=set([a[3:].replace("_"," ") for a in np.char.lower(pd.read_csv('/home/hag007/Desktop/aggregate_gwas_report/gs_bc/GO_BP_200.tsv', sep='\t')['gene_set_name'].values.astype(np.str))])
    gs=set([a[3:].lower().strip().replace("_"," ") for a in open('/home/hag007/Desktop/aggregate_gwas_report/gs_bc/gold_standard_terms.tsv').readlines()])

    return bg,gs


def check_intersection(base_folder, algos, datasets, bg,gs):
    results=[]
    for cur_ds in datasets:
        for cur_alg in algos:

            df=pd.read_csv(
                    "/home/hag007/Desktop/aggregate_gwas_report/oob/emp_diff_modules_{}_{}_passed_oob_reduced.tsv".format(
                        cur_ds, cur_alg), sep='\t').dropna()
            alg_set=set([a.replace("_", " ") for a in np.char.lower(
                df.loc[df["passed_oob_permutation_test"].apply(
                    lambda a: np.any(np.array(a[1:-1].split(", ")) == "True")).values, :].sort_values(by=["hg_pval_max"], ascending=False)['GO name'].values.astype(np.str))])

            # sumbit_to_revigo(algorithm=cur_alg, dataset=cur_ds)
            # df=pd.read_csv(
            #         "/home/hag007/Desktop/aggregate_gwas_report/oob/REVIGO.csv".format(
            #             cur_ds, cur_alg), sep=',').sort_values(by=["log10 p-value"], ascending=True)
            # df=df.loc[df['eliminated']==0]
            # alg_set = set([a.replace("_", " ") for a in np.char.lower(
            #     df['description'].values.astype(np.str))])


            bg_overlap=alg_set.intersection(bg)
            gs_overlap=alg_set.intersection(gs)
            print "dataset: {}, algo: {}. OS: {}/{} (n={})".format(cur_ds, cur_alg, len(gs_overlap), len(bg_overlap), len(alg_set))
            print "fraction before filtering bg: {}".format(float(len(gs_overlap))/max(len(alg_set),1.0))
            print "fraction after filtering bg: {}".format(float(len(gs_overlap))/max(len(bg_overlap),1.0))
            results.append((bg_overlap, gs_overlap, alg_set))

    return results


def plot_retrieved_terms(datasets, algos, results):
    for i, cur_ds in enumerate(datasets):
        for j, cur_alg in enumerate(algos):
            cur_result=results[i*len(datasets)+j]
            counter=0
            xs=[]
            x=0
            ys = []
            auc=0.0
            for cur_set in cur_result[2]:
                if cur_set not in cur_result[0]: continue
                if cur_set in cur_result[1]:
                    counter +=1
                ys.append(counter)
                xs.append(x)
                x+=1
                auc += counter

            y_l=0
            x_l=0
            if len(ys)!=0:
               y_l=ys[-1]
            if len(xs)!=0:
                x_l=xs[-1]

            hg_score=hypergeom.sf(y_l, 3856, 289, x_l)  + hypergeom.pmf(y_l, 3856, 289, x_l)
            print hg_score


            plt.cla()
            plt.plot(xs,ys, label="aggregated gold standard terms")
            plt.plot([0,x_l],[0,y_l], 'r--', label="diagonal")
            plt.xlabel("ranked algo. terms (bg)")
            plt.ylabel("GS terms")
            plt.title("normalized AUC: {}\ntotal fraction of retrieved GS terms {}, HG pval: {}".format(round(auc/max((len(xs)*(len(xs)+1)/2.0),1.0),2), round(float(y_l)/max(x_l,1.0),2), '{:0.3e}'.format(hg_score)))
            plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,"retrieving_score_{}_{}.png".format(cur_ds,cur_alg)))
            plt.legend()





if __name__=="__main__":
    algos=["my_netbox_td"] # ,'keypathwayminer_INES_GREEDY','jactivemodules_greedy','jactivemodules_sa','netbox','bionet', "my_netbox_td"
    datasets=['Breast_Cancer2.G50']
    bg, gs=load_gs()

    results=check_intersection('/home/hag007/Desktop/aggregate_gwas_report/oob', algos, datasets, bg,gs)

    plot_retrieved_terms(datasets, algos, results)