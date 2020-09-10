import math
import random
import sys
sys.path.insert(0, '../')

import networkx as nx

import pandas as pd

from fastsemsim.SemSim import *
from fastsemsim.Ontology import ontologies
from fastsemsim.Ontology import AnnotationCorpus
from fastsemsim.SemSim.SetSemSim import SetSemSim

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as ml_colors
from matplotlib.lines import Line2D

from rpy2.robjects import pandas2ri
pandas2ri.activate()

import constants

import multiprocessing

from utils.daemon_multiprocessing import func_star

import argparse

# import utils.go_hierarcies as go_hierarcies

import simplejson as json

from networkx.algorithms.components import connected_component_subgraphs as connected_component_subgraphs

from networkx.algorithms.community.quality import performance, coverage

from utils.go_similarity import calc_intra_similarity

def cc_to_v_ratio_analysis(cache_file, intervals=np.linspace(0,10,101)):
    vertices = []
    edges=[]
    G = nx.Graph()
    adj = json.load(open(cache_file, 'r'))
    for k, v in adj.iteritems():
        if float(v) > 0:
            edges.append((k.split("_")[0], k.split("_")[1], float(v)))
            vertices.append(k.split("_")[0])
            vertices.append(k.split("_")[1])

    vertices = np.unique(vertices)
    G.add_nodes_from(vertices)
    edges.sort(key=lambda a: a[2], reverse=True)

    i_edge = 0
    cc_to_v_ratio = []
    for cur_th in np.flip(intervals):
        while i_edge < len(edges):
            if edges[i_edge][2] >= cur_th:
                G.add_edge(edges[i_edge][0], edges[i_edge][1], weight=edges[i_edge][2])
                i_edge += 1
            else:
                break

        if len(vertices) > 0:
            # cc_to_v_ratio.append(len(list(connected_component_subgraphs(G))) / float(len(vertices)))
            cc_to_v_ratio.append(len(list(connected_component_subgraphs(G))) / float(len(vertices)))
        else:
            cc_to_v_ratio.append(1)

    return cc_to_v_ratio


def community_performance_analysis(cache_file, intervals=np.linspace(1,10,11), modules=[], pf=3, dataset=None, algo=None, base_folder=None, file_format = None, sim_method="Resnik"):
    vertices = []
    edges=[]
    G = nx.Graph()

    print "current cur_algo: {}".format(algo)

    if algo is None or base_folder is None or file_format is None:

        adj = json.load(open(cache_file, 'r'))

    else:

        try:

            emp_results = pd.read_csv(
                os.path.join(base_folder,
                             file_format.format(dataset, algo)), sep='\t', index_col=0)

            emp_results = emp_results.sort_values(by='emp_rank')
            emp_results_fdr = emp_results.dropna().loc[emp_results.dropna()["passed_oob_permutation_test"].apply(
                lambda a: np.any(np.array(a[1:-1].split(", ")) == "True")).values, :]

            all_go_terms = emp_results_fdr.index.values

            all_go_terms_r, all_go_terms_o, adj = calc_intra_similarity(all_go_terms, pf, emp_results_fdr, cache_file,
                                                                        sim_method)
            adj=dict(adj)

        except Exception,e :
            print e
            print "could not find {}".format(os.path.join(base_folder,
                                                          file_format.format(dataset, algo)), dataset, algo)
            return

    for k, v in adj.iteritems():
        if float(v) > 0:
            edges.append((k.split("_")[0], k.split("_")[1], float(v)))
            vertices.append(k.split("_")[0])
            vertices.append(k.split("_")[1])

    vertices = np.unique(vertices)
    G.add_nodes_from(vertices)
    edges.sort(key=lambda a: a[2], reverse=True)

    i_edge = 0
    performaces = []
    for cur_th in np.flip(intervals):
        while i_edge < len(edges):
            if edges[i_edge][2] >= cur_th:
                G.add_edge(edges[i_edge][0], edges[i_edge][1], weight=edges[i_edge][2])
                i_edge += 1
            else:
                break

        if len(vertices) > 0:

            if len(G.edges) < len(G.nodes) : continue #

            are_same_modules= modules[0] == modules[1]
            modules[0] = list(set(modules[0]))
            excluded_graph = list(set(G.nodes).difference(modules[0]))
            module_r = list(set(modules[0]).difference(set(modules[0]).difference(G.nodes)))
            # module_r=module[0]
            if len(module_r) < len(modules[0]):
                print "warning: {} terms of module are not included in graph".format(len(modules[0]) - len(module_r))

            modules[0] = module_r

            global_fraction = len(G.edges) / float(len(G.nodes) * (len(G.nodes) - 1))
            if are_same_modules:
                    intra_fraction = len(G.subgraph(modules[0]).edges) / max(float(len(modules[0]) * (len(modules[0]) - 1)), 10e-10)
                    performaces.append(min(10e10, intra_fraction / global_fraction))
                    # performaces.append(coverage(G, [module[0], excluded_graph]))

            else:
                # modules_union=list(set(modules[0] + modules[1]))
                # m0_reduced=list(set(modules[0]).difference(modules[1]))
                # m1_reduced = list(set(modules[1]).difference(modules[0]))
                # m_intersection=set(modules[0]).intersection(modules[1])
                # between_edges_n=len(G.subgraph(modules_union).edges) - (len(G.subgraph(modules[0]).edges) + len(G.subgraph(modules[1]).edges) - len(G.subgraph(m_intersection).edges))
                #
                # inter_fraction=between_edges_n / max(1.0, float(len(m0_reduced)*len(m1_reduced)))
                #
                # performaces.append(((inter_fraction/global_fraction)
                #                 * (float(len(modules_union))/len(m0_reduced + m1_reduced))))

                modules[1]=excluded_graph
                between_edges_n=len(G.edges) - (len(G.subgraph(modules[0]).edges) + len(G.subgraph(modules[1]).edges))

                inter_fraction=between_edges_n / max(1.0, float(len(modules[0])*len(modules[1])))

                performaces.append(inter_fraction/global_fraction)


                # inter_fraction = between_edges_n / float(len(excluded_graph) * len(module))
                # performaces.append(inter_fraction / global_fraction)
                                    # * (float(len(modules_union)) / len(m0_reduced + m1_reduced))))



        # else:
        #     performaces.append(0)

    return performaces

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="TNFa_2,HC12,ROR_1,ERS_1,IEM") # "Breast_Cancer.G50,Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50 TNFa_2,HC12,SHERA,ROR_1,SHEZH_1,ERS_1,IEM"
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy,jactivemodules_sa,netbox,my_netbox_td,bionet,hotnet2,keypathwayminer_INES_GREEDY") # jactivemodules_greedy,jactivemodules_sa,bionet,netbox,my_netbox_td
    # parser.add_argument('--module_indices', dest='module_indices',
    #                     default="0,1,2")  # jactivemodules_greedy,jactivemodules_sa,bionet,netbox,my_netbox_td
    parser.add_argument('--pf', dest='pf', default=3)
    parser.add_argument('--base_folder', dest='base_folder', default=constants.OUTPUT_GLOBAL_DIR)
    parser.add_argument('--sim_method', dest='sim_method', default="Resnik")

    args = parser.parse_args()

    prefix = args.prefix
    base_folder = args.base_folder
    sim_method= args.sim_method
    datasets=["{}".format(x) for x in args.datasets.split(",")]
    # module_indices=args.module_indices.split(",")
    algos = args.algos.split(",")
    pf=int(args.pf)
    print "test"
    h_scores = pd.DataFrame()
    limit=11
    factor=0.1
    intervals=np.linspace(0,limit,limit/factor + 1)
    df_avg_cc=pd.DataFrame()
    for cur_ds in datasets:
        print "cur_ds: {}".format(cur_ds)
        for cur_alg in algos:
            print "cur_alg: {}".format(cur_alg)

            cache_file=os.path.join(constants.CACHE_GLOBAL_DIR, "similarity_cache_emp_{}_{}_{}.npy".format(cur_ds,cur_alg, sim_method))
            cc_to_v_ratio = cc_to_v_ratio_analysis(cache_file, intervals)

            agg_score=np.mean(cc_to_v_ratio[-31:-29])#np.sum([cc_to_v_ratio[-1-i]*factor+(cc_to_v_ratio[-2-i]-cc_to_v_ratio[-1-i])/2.0 for i, x in enumerate(np.linspace(0,limit,limit/factor + 1)[:-1])])/float(limit)#np.mean(cc_to_v_ratio)
            plt.plot(np.flip(intervals), cc_to_v_ratio, label="{} ({})".format(cur_alg, round(agg_score,2)))
            df_avg_cc.loc[cur_alg,cur_ds]=agg_score

        plt.legend()
        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "cc_2_v_ratio_{}_{}.png".format(cur_ds,"")))
        plt.clf()
        df_avg_cc.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "cc_2_v_ratio.tsv"), sep='\t')



