import sys
sys.path.insert(0, '../')

import pandas as pd

from fastsemsim.SemSim import *

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as ml_colors
from matplotlib.lines import Line2D

from rpy2.robjects import pandas2ri
pandas2ri.activate()

import constants

import argparse

from utils.go_similarity import calc_intra_similarity

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="TNFa_2,HC12,ROR_1,SHEZH_1,ERS_1,IEM") #  SHERA "Breast_Cancer.G50,Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50 TNFa_2,HC12,SHERA,ROR_1,SHEZH_1,ERS_1,IEM"
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy,jactivemodules_sa,bionet,netbox,my_netbox_td,hotnet2,keypathwayminer_INES_GREEDY") #
    # parser.add_argument('--module_indices', dest='module_indices',
    #                     default="0,1,2")  # jactivemodules_greedy,jactivemodules_sa,bionet,netbox,my_netbox_td
    parser.add_argument('--pf', dest='pf', default=3)
    parser.add_argument('--base_folder', dest='base_folder', default=os.path.join(constants.OUTPUT_GLOBAL_DIR,"emp_fdr","MAX"))
    parser.add_argument('--sim_method', dest='sim_method', default="Resnik")
    parser.add_argument('--file_format', dest='file_format', default="emp_diff_modules_{}_{}_passed_oob.tsv")
    parser.add_argument('--cutoffs', dest='cutoffs', default="0.0, 1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0")

    args = parser.parse_args()

    prefix = args.prefix
    base_folder = args.base_folder
    sim_method= args.sim_method
    file_format=args.file_format
    # semsim = SetSemSim(ontology, ac, TSS=sim_method, MSS="BMA")
    datasets=["{}".format(x) for x in args.datasets.split(",")]
    # module_indices=args.module_indices.split(",")
    algos = args.algos.split(",")
    cutoffs = np.array(args.cutoffs.split(','),dtype=float)
    pf=int(args.pf)
    print "test"
    h_scores = pd.DataFrame()
    df_full_data = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "modules_full_data_{}.tsv".format(prefix)), index_col=0, sep='\t')

    for cutoff in cutoffs:
        suffix=str(cutoff)
        for cur_ds in datasets:
            for cur_alg in algos:

                ####

                try:

                    emp_results = pd.read_csv(
                        os.path.join(base_folder,
                                     file_format.format(cur_ds, cur_alg)), sep='\t', index_col=0)

                except:
                    print "could not find {}".format(os.path.join(base_folder,
                                                                  file_format.format(cur_ds, cur_alg)), cur_ds, cur_alg)
                    continue

                emp_results = emp_results.sort_values(by='emp_rank')
                emp_results_fdr = emp_results.dropna().loc[emp_results.dropna()["passed_oob_permutation_test"].apply(
                    lambda a: np.any(np.array(a[1:-1].split(", ")) == "True")).values, :]

                #####


                module_indices=[a.split("_")[-1] for a in df_full_data.loc[(df_full_data['algo']== cur_alg).values & (df_full_data['dataset']== cur_ds).values & (df_full_data['EHR']>0.2).values].index.values]
                for cur_module_index_0 in module_indices:
                    set_0=df_full_data.loc["{}_{}_module_{}".format(cur_ds,cur_alg,cur_module_index_0), "tp"]
                    if not type(set_0) is str:
                        set_0=[]
                    else:
                        set_0=[a.split(": ")[0] for a in df_full_data.loc["{}_{}_module_{}".format(cur_ds,cur_alg,cur_module_index_0), "tp"].split("\n")]
                    for cur_module_index_1 in module_indices:

                        if cur_module_index_0 > cur_module_index_1:
                            continue

                        set_1 = df_full_data.loc["{}_{}_module_{}".format(cur_ds, cur_alg, cur_module_index_1), "tp"]
                        if not type(set_1) is str:
                            set_1 = []
                        else:
                            set_1 = [a.split(": ")[0] for a in df_full_data.loc["{}_{}_module_{}".format(cur_ds, cur_alg, cur_module_index_1), "tp"].split("\n")]

                        print "current dataset: {} {} {} {}".format(cur_ds, cur_alg, cur_module_index_0, cur_module_index_1)
                        if cur_module_index_0<cur_module_index_1:
                            cache_file=os.path.join(constants.CACHE_GLOBAL_DIR, "similarity_cache_{}_{}_{}_{}_{}.npy".format(sim_method, cur_ds,cur_alg, cur_module_index_0, cur_module_index_1))
                        else:
                            cache_file=os.path.join(constants.CACHE_GLOBAL_DIR, "similarity_cache_{}_{}_{}_{}_{}.npy".format(sim_method, cur_ds, cur_alg, cur_module_index_1, cur_module_index_0))

                        # cc_to_v_ratio=cc_to_v_ratio_analysis(cache_file)
                        # # homogeneity_score = 1-np.mean(cc_to_v_ratio)
                        # factor=0.1
                        # limit=10
                        # homogeneity_score = 1 - np.sum([cc_to_v_ratio[-1-i]*factor+(cc_to_v_ratio[-2-i]-cc_to_v_ratio[-1-i])/2.0 for i, x in enumerate(np.linspace(0,limit,limit/factor + 1)[:-1])])/float(limit)

                        # adj=calc_similarity_matrix(set_0, set_1, pf=pf, cache_file=cache_file, sim_method=sim_method)
                        homogeneity_score=0
                        if cur_module_index_0 == cur_module_index_1:
                            all_go_terms_r, all_go_terms_o, adj = calc_intra_similarity(None, pf, emp_results_fdr, cache_file, sim_method, cutoff, set_0, set_1)
                            # homogeneity_score = (len(all_go_terms_o)-len(all_go_terms_r))/float(min(len(set_0), len(set_1)))
                            homogeneity_score = len(all_go_terms_r)
                            # homogeneity_score=np.mean([a for a in adj.values() if a >=0]) if len(adj) >0 else -1


                        h_scores.loc["{}_{}_modules_{}_{}".format(cur_ds, cur_alg, cur_module_index_0, cur_module_index_1), "homogeneity_score"]=homogeneity_score
                        h_scores.loc["{}_{}_modules_{}_{}".format(cur_ds, cur_alg, cur_module_index_0,
                                                                  cur_module_index_1), "dataset"] = cur_ds
                        h_scores.loc["{}_{}_modules_{}_{}".format(cur_ds, cur_alg, cur_module_index_0,
                                                                  cur_module_index_1), "algo"] = cur_alg
                        h_scores.loc["{}_{}_modules_{}_{}".format(cur_ds, cur_alg, cur_module_index_0,
                                                                  cur_module_index_1), "module_0"] = cur_module_index_0
                        h_scores.loc["{}_{}_modules_{}_{}".format(cur_ds, cur_alg, cur_module_index_0,
                                                                  cur_module_index_1), "module_1"] = cur_module_index_1




        h_scores.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "h_scores_{}.tsv".format(prefix)), sep='\t')

        h_scores=pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "h_scores_{}.tsv".format(prefix)), sep='\t',index_col=0)

        averaged_heterogeneities=[]
        averaged_homogeneities=[]
        cs=[]
        txts=[]
        fig, ax = plt.subplots()
        df_homogeneity_avg=pd.DataFrame()
        df_heterogeneity_avg = pd.DataFrame()
        for cur_ds in datasets:
            for cur_alg in algos:
                ds_filtered_scores=h_scores.loc[(h_scores['algo']== cur_alg).values & (h_scores['dataset']== cur_ds).values]
                averaged_heterogeneities.append(ds_filtered_scores[ds_filtered_scores["module_0"] !=ds_filtered_scores["module_1"]]['homogeneity_score'].values.mean())
                averaged_homogeneities.append(ds_filtered_scores[ds_filtered_scores["module_0"] == ds_filtered_scores["module_1"]]['homogeneity_score'].values.mean())
                df_homogeneity_avg.loc[cur_alg, cur_ds]=averaged_homogeneities[-1]
                df_heterogeneity_avg.loc[cur_alg, cur_ds] = averaged_heterogeneities[-1]
                cs.append(algos.index(cur_alg)/float(len(algos)))
                txts.append(cur_ds)
                # labels.append()

        df_homogeneity_avg.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "homogeneity_avg_matrix_{}_{}.tsv".format(prefix,suffix)), sep='\t')
        df_heterogeneity_avg.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "heterogeneity_avg_matrix_{}_{}.tsv".format(prefix,suffix)), sep='\t')

        averaged_homogeneities=np.array(averaged_homogeneities)
        averaged_heterogeneities= np.array(averaged_heterogeneities)
        cs = np.array(cs)
        txts=np.array(txts)
        # cs=cs[~np.isnan(averaged_heterogeneities)]
        # averaged_homogeneities=averaged_homogeneities[~np.isnan(averaged_heterogeneities)]
        # txts=txts[~np.isnan(averaged_heterogeneities)]
        # averaged_heterogeneities=averaged_heterogeneities[~np.isnan(averaged_heterogeneities)]
        DEFAULT_VAL=0
        ax.scatter([a if not np.isnan(a) else DEFAULT_VAL for a in averaged_heterogeneities], [a if not np.isnan(a) else DEFAULT_VAL for a in averaged_homogeneities], c=cs, cmap='jet')
        for i,data in enumerate(zip([a if not np.isnan(a) else DEFAULT_VAL for a in averaged_heterogeneities], [a if not np.isnan(a) else DEFAULT_VAL for a in averaged_homogeneities])):
            x, y = data
            ax.annotate(txts[i], (x,y))
        ax.set_xlabel("avg_heterogeneity (var={})".format(round(np.var([filter(lambda a : not np.isnan(a), averaged_heterogeneities)]),3)))
        ax.set_ylabel("avg_homogeneity (var={})".format(round(np.var([filter(lambda a : not np.isnan(a), averaged_homogeneities)]),3)))

        cmap = plt.cm.jet
        colorlist = [ml_colors.rgb2hex(cmap(a / float(np.size(algos) - 1))) for a in np.arange(np.size(algos))]
        patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                          markerfacecolor=c) for i, a, c in zip(list(range(len(algos))), algos, colorlist)]
        ax.legend(handles=list(reversed(patches)), loc='lower right', framealpha=0.5)


        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,"h_plot_{}_{}.png".format(prefix,suffix)))



