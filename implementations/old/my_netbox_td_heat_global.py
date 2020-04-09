import sys

sys.path.insert(0, "../")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pcst_fast
from networkx.algorithms.community.quality import modularity


from networkx.algorithms.community.centrality import girvan_newman

from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import fdrcorrection0

import os
import constants

from networkx.algorithms.components import connected_component_subgraphs as connected_component_subgraphs

from implementations.my_netbox_scoring_system import MyNetboxScoringSystems

from implementations.modules_reader import read_preprocessed_slices

from heat_diffusion_service import HeatDiffusionService

from scipy.stats import zscore

import random

import multiprocessing
from utils.daemon_multiprocessing import func_star


pertubed_node_score=1

def calc_emp_pval(cur_rv, cur_dist):
    cur_dist = np.array(cur_dist, np.float32)
    emp_pvals = []

    hg_pvals = np.array([cur_rv], dtype=np.float32)

    for hg_pval in hg_pvals:
        pos = np.size(cur_dist) - np.searchsorted(np.sort(cur_dist), hg_pval, side='left')
        emp_pval = pos / float(np.size(cur_dist))
        if emp_pval==0:
            emp_pval=1/float(len(cur_dist))
        emp_pvals.append(emp_pval)

    return emp_pvals


def extract_scores(scores_file):
    """"""
    scores = pd.read_csv(scores_file, sep='\t', index_col=0, header=None)
    if "pval" in scores.columns:
        scores["score"] = scores["pval"]
    else:
        scores["score"] = 1
    return scores


def adjust_scores(scores, nb_scoring_system):
    return nb_scoring_system.calculate_adjusted_scores(scores)


def build_network(network_file):
    """"""
    edges_dataset = pd.read_csv(network_file, sep='\t', header=None)
    edges = []
    for ind, row in edges_dataset.iterrows():
        # if row.iloc[0]!=row.iloc[2]:
        edges.append((row.iloc[0], row.iloc[2]))
    G = nx.Graph()
    G.add_edges_from(edges)
    nx.set_node_attributes(G, 0, 'adjusted_score')
    nx.set_node_attributes(G, 0, 'score')

    return G


def add_scores_to_nodes(G, scores):
    """"""
    inds = []
    for ind, row in scores.iterrows():
        if ind in G.nodes:
            inds.append(ind)
            G.nodes[ind]["score"] = row["score"]
            G.nodes[ind]["adjusted_score"] = row["adjusted_score"]

    print "inds: {}".format(len(set(inds)))
    return G


def mark_extracted_nodes(G, nb_scoring_system):
    pertubed_nodes = []
    for cur_node in G.nodes:
        if G.nodes[cur_node]["adjusted_score"] >= nb_scoring_system.get_threshold():
            G.nodes[cur_node]["extract_node"] = True
            G.nodes[cur_node]["pertubed_node"] = pertubed_node_score
            G.nodes[cur_node]["color"] = "red"
            pertubed_nodes.append(cur_node)
        else:
            G.nodes[cur_node]["pertubed_node"] = 0

    print "total # of pertubed genes: {}".format(len(pertubed_nodes))
    linker_pvals = {}
    for cur_node in G.nodes:
        if not G.nodes[cur_node]["pertubed_node"]:
            pertubed_neighbors = 0
            for cur_neighbor in G.neighbors(cur_node):
                if G.nodes[cur_neighbor]["pertubed_node"]:
                    pertubed_neighbors += 1
            q = 1
            if pertubed_neighbors > 1:
                linker_pvals[cur_node] = hypergeom.sf(pertubed_neighbors, len(G.nodes), len(pertubed_nodes),
                                                      len(list(G.neighbors(cur_node)))) \
                                         + hypergeom.pmf(pertubed_neighbors, len(G.nodes), len(pertubed_nodes),
                                                         len(list(G.neighbors(cur_node))))
                G.nodes[cur_node]["linker_pval"] = linker_pvals[cur_node]
            else:
                G.nodes[cur_node]["extract_node"] = False

    fdr_bh_results = fdrcorrection0(linker_pvals.values(), alpha=nb_scoring_system.linker_threshold, method='indep',
                                    is_sorted=False)

    print "total # of linkers: {}".format(sum(fdr_bh_results[0]))

    for i, cur_key in enumerate(linker_pvals.keys()):
        G.nodes[cur_key]["linker_qval"] = fdr_bh_results[1][i]
        G.nodes[cur_key]["passed_fdr"] = fdr_bh_results[0][i]
        G.nodes[cur_key]["extract_node"] = fdr_bh_results[0][i]
        G.nodes[cur_key]["color"] = "blue"

    included_edges = []
    for cur_edge in G.edges:
        if cur_edge[0] in G.nodes and G.nodes[cur_edge[0]]["extract_node"] and cur_edge[1] in G.nodes and \
                G.nodes[cur_edge[1]]["extract_node"]:
            included_edges.append(cur_edge)

    print "total # of edges: {}".format(len(included_edges))

    return G


def extract_modules(G):
    """"""
    G_extracted_modules = G.copy()

    print "total # of nodes: {}".format(len(G_extracted_modules.nodes))
    for cur_node in G.nodes:
        if not G.nodes[cur_node]["extract_node"]:
            G_extracted_modules.remove_node(cur_node)
            # print "total # of nodes after node removal: {}".format(len(G_extracted_modules.nodes))
            G.nodes[cur_node]["extract_node"] = False
        else:
            print "node {} is : {}".format(cur_node,
                                           "pertubed_node" if G.nodes[cur_node]["pertubed_node"] else "linker")

    print "total # of nodes after nodes removal: {}".format(len(G_extracted_modules.nodes))

    print "# of cc before modularity optimization: {}".format(
        len(list(connected_component_subgraphs(G_extracted_modules))))

    # plt.subplots(1,1,figsize=(50,50))
    # nx.draw_networkx(G_extracted_modules, pos=nx.spring_layout(G_extracted_modules))
    # plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,"cc_before_modularity.png"))

    return G_extracted_modules


def optimize_modularity(G, improvement_delta=0):
    """"""
    G_optimized = G.copy()

    G_optimized.remove_edges_from(list(nx.selfloop_edges(G_optimized)))
    G_optimized.remove_nodes_from(list(nx.isolates(G_optimized)))
    print "after optimizing for self loops - number of edges: {}, nodes: {}".format(len(G_optimized.edges),
                                                                                    len(G_optimized.nodes))

    reached_optimal_modularity = False
    optimal_modularity = -1

    pertubed_nodes = []
    for cur_node in G.nodes:
        if G.nodes[cur_node]["pertubed_node"]:
            pertubed_nodes.append(cur_node)

    cur_components = list(connected_component_subgraphs(G_optimized))
    while not reached_optimal_modularity:

        for cur_cc in cur_components:
            if len(cur_cc) < 4:
                G_optimized.remove_nodes_from(cur_cc)
            else:
                pertubed_nodes_in_cc = []
                for cur_node in cur_cc:
                    if G.nodes[cur_node]["pertubed_node"]:
                        pertubed_nodes_in_cc.append(cur_node)


                sig_score = hypergeom.sf(len(pertubed_nodes_in_cc), len(G.nodes), len(pertubed_nodes),
                                         len(cur_cc)) \
                            + hypergeom.pmf(len(pertubed_nodes_in_cc), len(G.nodes), len(pertubed_nodes),
                                            len(cur_cc))

                if len(cur_cc) < 8 and sig_score > 0.05:
                    G_optimized.remove_nodes_from(cur_cc)

        print "after optimizing for sig modules - number of edges: {}, nodes: {}".format(len(G_optimized.edges),
                                                                                         len(G_optimized.nodes))
        cur_components = list(connected_component_subgraphs(G_optimized))
        print "number of cc: ", len(cur_components)

        optimized_connected_components = girvan_newman(G_optimized)
        cur_components = sorted(next(optimized_connected_components))

        cur_modularity = modularity(G_optimized, cur_components, weight='weight')
        if cur_modularity > optimal_modularity + improvement_delta:
            optimal_modularity = cur_modularity
            optimal_components = cur_components

            edges_to_remove = []
            for cur_edge in G_optimized.edges:
                included = False
                for cur_cc in optimal_components:
                    if cur_edge[0] in cur_cc and cur_edge[1] in cur_cc:
                        included = True
                if not included:
                    edges_to_remove.append(cur_edge)

            G_optimized.remove_edges_from(edges_to_remove)

            print "new: ", len(cur_components), len([a for a in cur_components if len(a) > 3]), optimal_modularity



        else:
            reached_optimal_modularity = True

    nodes_to_remove = []
    for cur_node in G_optimized.nodes:
        if len(list(G_optimized.neighbors(cur_node))) == 0:
            nodes_to_remove.append(cur_node)

    G_optimized.remove_nodes_from(nodes_to_remove)

    # plt.subplots(1, 1, figsize=(50, 50))
    # nx.draw_networkx(G_optimized, pos=nx.spring_layout(G_optimized), node_color=[G_optimized.nodes[k]["color"] for k in G_optimized.nodes])
    # plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "cc_after_modularity.png"))
    print "modularity: ", modularity(G_optimized, list(connected_component_subgraphs(G_optimized)), weight='weight')
    return G_optimized, connected_component_subgraphs(G_optimized)


def filter_modules(G_modularity, connected_components):
    """"""
    pertubed_nodes = []
    passed_ccs=[]
    for cur_node in G_modularity.nodes():
        if G_modularity.nodes[cur_node]["pertubed_node"]:
            pertubed_nodes.append(cur_node)

    sig_scores = []
    for cur_cc in connected_components:
        pertubed_nodes_in_cc = []

        for cur_node in cur_cc:
            if G_modularity.nodes[cur_node]["pertubed_node"]:
                pertubed_nodes_in_cc.append(cur_node)

        sig_scores.append(hypergeom.sf(len(pertubed_nodes_in_cc), len(G_modularity.nodes), len(pertubed_nodes),
                                       len(cur_cc)) \
                          + hypergeom.pmf(len(pertubed_nodes_in_cc), len(G_modularity.nodes), len(pertubed_nodes),
                                          len(cur_cc)))

    print "sig_scores: {}".format(sig_scores)
    fdr_bh_results = fdrcorrection0(sig_scores, alpha=0.05, method='indep',
                                    is_sorted=False)
    print "adjust_sig_scores: {}\n(n={})".format(fdr_bh_results[1], sum(fdr_bh_results[0]))

    for cur_cc, sig_score in zip(connected_components, sig_scores):
        if sig_score <=0.05:
            passed_ccs.append(cur_cc)



def prune_network_by_modularity(G, modules):
    G_modularity=G.copy()
    edges_to_remove=[]
    for cur_edge in G_modularity.copy().edges:
        in_cc = False
        for cur_module in modules:
            if cur_edge[0] in cur_module and cur_edge[1] in cur_module:
                in_cc=True
        if not in_cc:
            edges_to_remove.append(cur_edge)

    G_modularity.remove_edges_from(edges_to_remove)

    return G_modularity

def leave_significant_modules(G_original, G_modularity, module_sig_th):

    pertubed_nodes = []
    for cur_node in G_modularity.nodes():
        if G_modularity.nodes[cur_node]["pertubed_node"]:
            pertubed_nodes.append(cur_node)

    connected_components=list(connected_component_subgraphs(G_modularity))
    edges_to_remove = []
    sig_scores=[]
    passed_modules=[]
    large_modules=[]
    G_optimized=G_modularity.copy()
    for cur_cc in connected_components:
        if len(cur_cc) < 4:
            G_optimized.remove_nodes_from(cur_cc)
        else:
            large_modules.append(cur_cc)
            pertubed_nodes_in_cc = []
            for cur_node in cur_cc:
                if G_modularity.nodes[cur_node]["pertubed_node"]: # or G_modularity.nodes[cur_node]["extract_node"]:
                    pertubed_nodes_in_cc.append(cur_node)


            sig_scores.append(hypergeom.sf(len(pertubed_nodes_in_cc), len(G_original.nodes), len(pertubed_nodes),
                                           len(cur_cc)) \
                        + hypergeom.pmf(len(pertubed_nodes_in_cc), len(G_original.nodes), len(pertubed_nodes),
                                        len(cur_cc)))

            # print "hypergeom params: {}, {}, {}, {} = {}".format(len(pertubed_nodes_in_cc), len(G_original.nodes),
            #                                                 len(pertubed_nodes),
            #                                                 len(cur_cc), sig_scores[-1])

    fdr_bh_results = fdrcorrection0(sig_scores, alpha=module_sig_th, method='indep',
                                    is_sorted=False)

    print "# of significant modules: ", np.sum(fdr_bh_results[0])

    i=0
    for cur_cc, is_passed_th, qval in zip(large_modules, fdr_bh_results[0], fdr_bh_results[1]):
        if is_passed_th:
            print "module #{} (n={}, qval={}):\n{}".format(i, len(cur_cc), qval, list(cur_cc.nodes))
            passed_modules.append(list(cur_cc.nodes))
        else:
            G_optimized.remove_nodes_from(cur_cc)

        i += 1

    return G_optimized, passed_modules, fdr_bh_results[1]



def retrieve_modules(G, optimized_connected_components):
    """"""

    print "# of cc after modularity optimization: {}".format(len(optimized_connected_components))
    for i, cur_cc in enumerate(optimized_connected_components):
        if len(cur_cc) > 3:
            print "cc #{}: n={}\n{}".format(i, len(cur_cc), cur_cc)

    optimized_connected_components



def main(dataset_name="PASCAL_SUM_Breast_Cancer2.G50", score_file_name="/media/hag007/Data/bnet/datasets/GE_HC12/output/ge_list.txt", network_file_name="dip.sif", modules_file_name=os.path.join(constants.NETWORKS_DIR, "dip_ng_modularity_components.txt"), module_sig_th=0.05):
    constants.update_dirs(DATASET_NAME_u=dataset_name)
    network_file_name = os.path.join(constants.NETWORKS_DIR, "dip.sif")

    # scores = extract_scores(score_file_name)
    # nb_ss = MyNetboxScoringSystems(scores)
    #
    # scores = adjust_scores(scores, nb_ss)
    # G = build_network(network_file_name)
    # G = add_scores_to_nodes(G, scores)
    # G = mark_extracted_nodes(G, nb_ss)
    # G_optimized, optimized_connected_components = optimize_modularity(G)
    # filter_modules(G_optimized, optimized_connected_components)
    # retrieve_modules(G_optimized, optimized_connected_components)

    scores = extract_scores(score_file_name)
    nb_ss = MyNetboxScoringSystems(scores)

    scores = adjust_scores(scores, nb_ss)
    G = build_network(network_file_name)
    G = add_scores_to_nodes(G, scores)
    G = mark_extracted_nodes(G, nb_ss)
    modularity_connected_components = read_preprocessed_slices(modules_file_name)
    G_modularity=prune_network_by_modularity(G, modularity_connected_components)

    G_modularity, optimized_connected_components, qvals = leave_significant_modules(G, G_modularity, module_sig_th)
    # retrieve_modules(G_modularity, optimized_connected_components)
    # p_cc = linear_threshold(G, [n for n in G.nodes if G.nodes[n]['pertubed_node'] > 0], steps=1)[0]
    hds = HeatDiffusionService()

    scores={}
    G_cc = nx.Graph(G)
    p_cc = hds.diffusion(G_cc, 'pertubed_node', "output", False, 1.0)
    # lt_gene_sets = linear_threshold(p_cc, [n for n in p_cc.nodes if p_cc.nodes[n]['pertubed_node'] > 0], steps=1)
    lt_genes=[n for n in p_cc.nodes] #  lt_gene_sets[0]+lt_gene_sets[1]
    for p_node in p_cc.nodes:
        if p_node in lt_genes:
            scores[p_node] = p_cc.nodes[p_node]['output_heat']

    nodes = list(G_cc.nodes)
    labels = {n: G_cc.nodes[n] for n in nodes}
    ###########

    root = -1
    num_clusters = 1
    pruning = 'strong'  # 'none'
    verbosity_level = 0
    edges_grid = []
    for cur_edge in G_cc.edges:
        edges_grid.append([nodes.index(cur_edge[0]), nodes.index(cur_edge[1])])

    ###########

    p = multiprocessing.Pool(5)
    manager = multiprocessing.Manager()
    permuted_scores=manager.dict()
    for cur_node in nodes:
        permuted_scores[cur_node]=manager.list()
    n_interations=50
    params=[]
    for i_iteration in np.arange(n_interations):
        params.append([create_bg_sample, [G_cc, hds, i_iteration, labels, nodes, permuted_scores]])


    p.map(func_star, params)
    p.close()

    pvals = {}
    for cur_node in nodes:
        if cur_node in lt_genes:
            pvals[cur_node]= calc_emp_pval(scores[cur_node], np.array(permuted_scores[cur_node]))[0]
            print cur_node, pvals[cur_node], scores[cur_node]


    print "pvals: {}".format(np.sort(pvals.values()))
    print "zscores: {}".format(zscore(np.sort(pvals.values())))
    print "var: {}".format(np.var(pvals.values()))

    fdr_results = fdrcorrection0(pvals.values(), alpha=0.3, method='indep', is_sorted=False)
    qvals = {}
    diffusion_cc = []
    for cur_gene, is_sig, cur_q in zip(pvals.keys(), fdr_results[0], fdr_results[1]):
        qvals[cur_gene] = cur_q
        if is_sig:
            diffusion_cc.append(cur_gene)


    print "qvals: {}".format(np.sort(fdr_results[1]))
    print "graph shrinkage factor: {} ({}/{})".format(sum(fdr_results[0])/float(len(fdr_results[0])), sum(fdr_results[0]), len(fdr_results[0]))

    # sg=nx.subgraph(G, diffusion_cc)

    ####
    scores={k: -np.log10(v) for k, v in pvals.iteritems()}


    vertices_prizes = []
    for cur_node in nodes:
        vertices_prizes.append(scores[cur_node])

    ALPHA_FACTOR=20.0
    alphas = [a for a in np.arange(int(ALPHA_FACTOR)+1)]
    COST_FACTOR=0.7
    ys = []
    y_plateaus = []
    y_plateau_values = []
    y_edge_costs = []
    plateau_value = 0
    plateau_index = 0
    pcst_results = {}
    for i_alpha, alpha in enumerate(alphas):
        # alpha = 0.1

        edges_costs=[np.percentile(scores.values(), 100) * (COST_FACTOR ** alpha) for cur_edge in edges_grid]

        vertices, edges = pcst_fast.pcst_fast(edges_grid, vertices_prizes, edges_costs, root, num_clusters, pruning,
                                              verbosity_level)
        pcst_results[alpha] = (vertices, edges)
        ys.append(len(vertices) / float(len(nodes)))
        if np.abs(ys[-1] - plateau_value) > 10 ** -1:
            if alpha >= 1 and alpha - plateau_index > 1:
                y_plateaus.append((plateau_index, alphas[i_alpha - 1]))
                y_plateau_values.append(plateau_value)
                y_edge_costs.append(edges_costs[0])
            plateau_index = alpha
            plateau_value = ys[-1]

    plt.plot(alphas, ys)
    plateau_max = (0, 0)
    plateau_max_interval = 0
    plateau_max_value = 0
    plateau_max_edge_cost = 0
    for cur_plateau, cur_y_plateau_value , y_edge_cost in zip(y_plateaus, y_plateau_values, y_edge_costs):
        if plateau_max_interval < cur_plateau[1] - cur_plateau[0]:
            plateau_max_interval = cur_plateau[1] - cur_plateau[0]
            plateau_max = cur_plateau
            plateau_max_value = cur_y_plateau_value
            plateau_max_edge_cost=y_edge_cost

        plt.plot(list(cur_plateau), [cur_y_plateau_value, cur_y_plateau_value])

    plt.title("plateau_max_value: {}\n plateau_max: {}\n plateau cost max: {}".format(plateau_max_value, plateau_max, plateau_max_edge_cost), fontsize=10)
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "alphas_m_graph.png"))

    vertices, edges = pcst_results[plateau_max[-1]]
    print "PCST vertices: {} ({})".format(vertices, len(vertices))

    ####

    G_pcst = nx.Graph()
    G_pcst.add_edges_from([(nodes[edges_grid[e][0]], nodes[edges_grid[e][1]]) for e in edges])

    nx.set_node_attributes(G_pcst, {n: labels[n] for n in  G_pcst.nodes})

    G_optimized, modules=optimize_modularity(G_pcst, improvement_delta=10**-2)

    module_nodes=[]
    for i_cur_module, cur_module in enumerate(modules):
        module_nodes.append(list(cur_module.nodes))
        print "module_agg_score", i_cur_module, sum([scores[n] for n in cur_module.nodes])-len(cur_module.edges)*plateau_max_edge_cost



    return module_nodes


def create_bg_sample(G_cc, hds, i_iteration, labels, nodes, permuted_scores):
    print "graph bg iteration {}".format(i_iteration)
    c_nodes = list(nodes)
    random.shuffle(c_nodes)
    mapping = {k: v for k, v in zip(nodes, c_nodes)}
    p_cc = nx.relabel_nodes(G_cc, mapping)
    nx.set_node_attributes(p_cc, labels)
    p_cc = hds.diffusion(p_cc, 'pertubed_node', "output", False, 1)
    for cur_node in p_cc.nodes:
        permuted_scores[cur_node]+= [p_cc.nodes[cur_node]["output_heat"]]


if __name__=="__main__":
    all_qvals=[]
    datasets=["PASCAL_SUM_Breast_Cancer.G50"] # "PASCAL_SUM_Breast_Cancer.G50","PASCAL_SUM_Crohns_Disease.G50","PASCAL_SUM_Schizophrenia.G50","PASCAL_SUM_Triglycerides.G50","PASCAL_SUM_Type_2_Diabetes.G50",
    for cur_ds in datasets:
        res=main(dataset_name=cur_ds, score_file_name="/media/hag007/Data/bnet/datasets/{}/output/ge_list.txt".format(cur_ds))
        # all_qvals=all_qvals+list(qvals)
    # sns.distplot([a for a in all_qvals if a <1.0], norm_hist=False, kde=False,bins=100)
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "qval_dist.png"))
    print all_qvals