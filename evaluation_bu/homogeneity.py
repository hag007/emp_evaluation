import sys
sys.path.insert(0, '../')

import pandas as pd

from fastsemsim.SemSim import *

import constants
import argparse
import networkx as nx

from modules_report import summarize_modules_ehr
from utils.go_similarity import calc_intra_similarity

def calc_homogeneity(cache_file, dataset=None, algo=None , base_folder=None, cutoff=1, module=[], pf=3, file_format = None, sim_method="Resnik"):


    print "current cur_algo: {}".format(algo)

    try:
        emp_results = pd.read_csv(
            os.path.join(base_folder,
                         file_format.format(dataset, algo)), sep='\t', index_col=0)
    except Exception,e :
        print e
        print "could not find {}".format(os.path.join(base_folder, file_format.format(dataset, algo)), dataset, algo)
        return

    emp_results = emp_results.sort_values(by='emp_pval_max')
    emp_results_fdr = emp_results.dropna().loc[emp_results.dropna()["passed_oob_permutation_test"].apply(
    lambda a: np.any(np.array(a[1:-1].split(", ")) == "True")).values, :]

    all_go_terms = emp_results_fdr.index.values

    if len(all_go_terms) <=1:
        return np.nan

    vertices = []
    edges=[]
    G = nx.Graph()

    all_go_terms_r, all_go_terms_o, adj = calc_intra_similarity(all_go_terms, pf, emp_results_fdr, cache_file,
                                                                sim_method, reduce_list=False)

    for k, v in adj.iteritems():
        if float(v) > 0:
            edges.append((k.split("_")[0], k.split("_")[1], float(v)))
            vertices.append(k.split("_")[0])
            vertices.append(k.split("_")[1])

    vertices = np.unique(vertices)
    G.add_nodes_from(vertices)
    edges.sort(key=lambda a: a[2], reverse=True)

    i_edge = 0
    while i_edge < len(edges):
        if edges[i_edge][2] >= cutoff:
            G.add_edge(edges[i_edge][0], edges[i_edge][1], weight=edges[i_edge][2])
            i_edge += 1
        else:
            break

    if len(vertices) == 0:
        return np.nan

    module = list(set(module))
    module_r = list(set(module).difference(set(module).difference(G.nodes)))
    if len(module_r) < len(module):
        print "warning: {} terms of module are not included in graph".format(len(module) - len(module_r))

    global_fraction = len(G.edges) / float(len(G.nodes) * (len(G.nodes) - 1))
    intra_fraction = len(G.subgraph(module_r).edges) / max(float(len(module_r) * (len(module_r) - 1)), 10e-10)
    # print(intra_fraction)
    # print(global_fraction)
    return min(10e10, intra_fraction / max(global_fraction,10e-5))


def main(prefix, base_folder, sim_method, file_format, pf, datasets, algos, cutoffs, recalc_module_report):

    h_scores = pd.DataFrame()

    if recalc_module_report:
        summarize_modules_ehr(prefix, datasets, algos, base_folder)

    df_full_data = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "module_cache_files", "modules_full_data_{}.tsv".format(prefix)), index_col=0, sep='\t')

    for cutoff in cutoffs:
        for cur_ds in datasets:
            print "cur dataset: {}".format(cur_ds)
            for cur_alg in algos:
                print "cur cutoff: {}, dataset: {}, algo: {}".format(cutoff, cur_ds, cur_alg)
                # try:
                #     emp_results = pd.read_csv(os.path.join(base_folder, file_format.format(cur_ds, cur_alg)), sep='\t', index_col=0)
                #
                # except:
                #     print "could not find {}".format(os.path.join(base_folder, file_format.format(cur_ds, cur_alg)), cur_ds, cur_alg)
                #     continue
                #
                # emp_results = emp_results.sort_values(by='emp_rank')
                # emp_results_fdr = emp_results.dropna().loc[emp_results.dropna()["passed_oob_permutation_test"].apply(
                #     lambda a: np.any(np.array(a[1:-1].split(", ")) == "True")).values, :]

                module_indices=sorted([int(a.split("_")[-1]) for a in df_full_data.loc[(df_full_data['algo']== cur_alg).values & (df_full_data['dataset']== cur_ds).values & (df_full_data['EHR']>0.2).values].index.values])
                for cur_module_index in module_indices:
                    go_set=df_full_data.loc["{}_{}_module_{}".format(cur_ds,cur_alg,cur_module_index), "tp"]
                    if not type(go_set) is str:
                        go_set=[]
                    else:
                        go_set=[a.split(": ")[0] for a in df_full_data.loc["{}_{}_module_{}".format(cur_ds,cur_alg,cur_module_index), "tp"].split("\n")]

                    print "current module: {} {} {}".format(cur_ds, cur_alg, cur_module_index)
                    cache_file=os.path.join(constants.CACHE_GLOBAL_DIR, "similarity_cache_emp_{}_{}_{}.npy".format(cur_ds, cur_alg, sim_method))
                    score=calc_homogeneity(dataset=cur_ds, algo=cur_alg, base_folder=base_folder, cache_file=cache_file, cutoff=cutoff, module=go_set, file_format=file_format)

                    print "module score: {}".format(score)

                    h_scores.loc["{}_{}_modules_{}".format(cur_ds, cur_alg, cur_module_index), "homogeneity_score"]=score
                    h_scores.loc["{}_{}_modules_{}".format(cur_ds, cur_alg, cur_module_index), "dataset"] = cur_ds
                    h_scores.loc["{}_{}_modules_{}".format(cur_ds, cur_alg, cur_module_index), "algo"] = cur_alg
                    h_scores.loc["{}_{}_modules_{}".format(cur_ds, cur_alg, cur_module_index), "module"] = cur_module_index

        std_homogeneities=[]
        averaged_homogeneities=[]
        df_homogeneity_avg=pd.DataFrame()
        df_homogeneity_std= pd.DataFrame()
        for cur_ds in datasets:
            for cur_alg in algos:
                ds_filtered_scores=h_scores.loc[(h_scores['algo']== cur_alg).values & (h_scores['dataset']== cur_ds).values]
                std_homogeneities.append(ds_filtered_scores['homogeneity_score'].values.std())
                averaged_homogeneities.append(ds_filtered_scores['homogeneity_score'].values.mean())
                df_homogeneity_avg.loc[cur_alg, cur_ds]=averaged_homogeneities[-1]
                df_homogeneity_std.loc[cur_alg, cur_ds] = std_homogeneities[-1]

        h_scores.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation", "h_scores_{}_{}.tsv".format(prefix,str(cutoff))), sep='\t')
        df_homogeneity_avg.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation", "homogeneity_avg_matrix_{}_{}.tsv".format(prefix,str(cutoff))), sep='\t')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--prefix', dest='prefix', default="PASCAL_SUM")
    parser.add_argument('--datasets', dest='datasets', default="Breast_Cancer.G50,Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50,Coronary_Artery_Disease.G50,Bone_Mineral_Density.G50,Height1.G50,Age_Related_Macular_Degeneration.G50,Atrial_Fibrillation.G50") # Breast_Cancer.G50,Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50,Coronary_Artery_Disease.G50,Bone_Mineral_Density.G50,Height1.G50,Alzheimer.G50Age_Related_Macular_Degeneration.G50,Atrial_Fibrillation.G50 # TNFa_2,HC12,SHERA,ROR_1,SHEZH_1,ERS_1,IEM,APO,CBX,IFT
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy,jactivemodules_sa,bionet,netbox,domino_original")
    parser.add_argument('--pf', dest='pf', default=10)
    parser.add_argument('--base_folder', dest='base_folder', default=os.path.join(constants.OUTPUT_GLOBAL_DIR, "oob"))
    parser.add_argument('--sim_method', dest='sim_method', default="Resnik")
    parser.add_argument('--file_format', dest='file_format', default="emp_diff_modules_{}_{}_passed_oob.tsv")
    parser.add_argument('--cutoffs', dest='cutoffs', default="1.0,2.0,3.0,4.0") # 1.0,2.0,3.0,4.0  ,5.0,6.0,7.0,8.0,9.0,10.0,11.0
    parser.add_argument('--recalc_module_report', dest='recalc_module_report', default="true") # 1.0,2.0,3.0,4.0  ,5.0,6.0,7.0,8.0,9.0,10.0,11.0

    args = parser.parse_args()

    prefix = args.prefix
    base_folder = args.base_folder
    sim_method= args.sim_method
    file_format=args.file_format
    pf=int(args.pf)
    datasets=["{}".format(x) for x in args.datasets.split(",")]
    algos = args.algos.split(",")
    cutoffs = np.array(args.cutoffs.split(','),dtype=float)
    recalc_module_report=args.recalc_module_report=="true"

    datasets=["tnfa", "hc", "ror", "shera", "shezh", "ers", "iem", "apo", "cbx", "ift"] 
    algos=["DOMINO2", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"] 
    prefix="GE"
    main(prefix, base_folder, sim_method, file_format, pf, datasets, algos, cutoffs, recalc_module_report)

    datasets=["brca", "crh", "scz", "tri", "t2d", "cad", "bmd", "hgt", "amd", "af"]
    algos=["DOMINO2", "netbox", "jactivemodules_greedy", "jactivemodules_sa", "bionet", "keypathwayminer_INES_GREEDY", "hotnet2"] 
    prefix="PASCAL_SUM"
    main(prefix, base_folder, sim_method, file_format, pf, datasets, algos, cutoffs, recalc_module_report)


