import math
import random
import sys
sys.path.insert(0, '../')

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

import utils.go_hierarcies as go_hierarcies

import simplejson as json

ontology_type = 'GeneOntology'
ignore_parameters = {'ignore': {}}
source_type = 'obo'
source = os.path.join(os.path.join(constants.GO_DIR, constants.GO_FILE_NAME))

print "\n######################"
print "# Loading ontology... #"
print "######################\n"

ontology = ontologies.load(source=source, source_type=source_type, ontology_type=ontology_type,
                           parameters=ignore_parameters)

print "\n######################"
print "# Loading Annotation Corpus... #"
print "######################\n"
ac = AnnotationCorpus.AnnotationCorpus(ontology)
ac.parse(os.path.join(constants.GO_DIR, "goa_human.gaf"), "gaf-2.0")
ac.isConsistent()

print "\n#################################"
print "# Annotation corpus successfully loaded."
print "#################################\n"


ROOT='GO:0008150'
dict_result, go2geneids, geneids2go, entrez2ensembl = go_hierarcies.build_hierarcy(
    roots=[ROOT])
vertices = dict_result.values()[0]['vertices']

go_hierarcies
CUTOFF=7

def calc_similarity(mat_adj, i_x, i_y, x, y, semsim):
    key="{}_{}".format(x,y)
    key_inv="{}_{}".format(y,x)
    if mat_adj[key] != -200: return
    if(x==y):
       mat_adj[key]=-100
       return 

    mat_adj[key] = semsim.SemSim(x, y)

    if np.isnan(mat_adj[key]):
        mat_adj[key] = -100
    mat_adj[key_inv] = mat_adj[key]

def calc_similarity_matrix(set_0, set_1, pf, cache_file ,sim_method):
    cache_loaded = False
    if os.path.exists(cache_file):
        try:
            adj=json.load(open(cache_file,'r'))
            cache_loaded=True
        except:
            pass

    if not cache_loaded:
        semsim = SetSemSim(ontology, ac, TSS=sim_method, MSS="BMA")
        manager = multiprocessing.Manager()
        adj = manager.dict()
        for x in set_0:
            for y in set_1:
                adj["{}_{}".format(x, y)] = -200
        params = []
        for i_x, x in enumerate(set_0):
            for i_y, y in enumerate(set_1):
                calc_similarity(adj, i_x, i_y, x, y, semsim)
                # params.append([calc_similarity, [adj, i_x, i_y, x, y, semsim]])
        print "len(params): {}".format(len(params))


        p = multiprocessing.Pool(pf)
        p.map(func_star, params)
        p.close()
        p.join()

    open(cache_file, 'w+').write(json.dumps(dict(adj)))
    return adj

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="TNFa_2,HC12,SHERA,ROR_1,SHEZH_1,ERS_1,IEM") # "Breast_Cancer.G50,Crohns_Disease.G50,Schizophrenia.G50,Triglycerides.G50,Type_2_Diabetes.G50 TNFa_2,HC12,SHERA,ROR_1,SHEZH_1,ERS_1,IEM"
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy,netbox,my_netbox_td") # jactivemodules_greedy,jactivemodules_sa,bionet,netbox,my_netbox_td
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
    algos = list(np.sort(args.algos.split(",")))
    pf=int(args.pf)
    print "test"
    h_scores = pd.DataFrame()
    df_full_data = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, "modules_full_data_{}.tsv".format(prefix)), index_col=0, sep='\t')
    for cur_ds in datasets:
        for cur_alg in algos:
            module_indices=[a.split("_")[-1] for a in df_full_data.loc[(df_full_data['algo']== cur_alg).values & (df_full_data['dataset']== cur_ds).values & (df_full_data['EHR']>0.2).values].index.values]
            for cur_module_index_0 in module_indices:
                set_0=df_full_data.loc["{}_{}_module_{}".format(cur_ds,cur_alg,cur_module_index_0), "tp"]
                if not type(set_0) is str:
                    set_0=[]
                else:
                    set_0=[a.split(": ")[0] for a in df_full_data.loc["{}_{}_module_{}".format(cur_ds,cur_alg,cur_module_index_0), "tp"].split("\n")]
                for cur_module_index_1 in module_indices:
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
                    adj=calc_similarity_matrix(set_0, set_1, pf=pf, cache_file=cache_file, sim_method=sim_method)
                    h_scores.loc["{}_{}_modules_{}_{}".format(cur_ds, cur_alg, cur_module_index_0, cur_module_index_1), "homogeneity_score"]=np.mean([a for a in adj.values() if a >=0]) if len(adj) >0 else -1
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
    algos_in_scatter=[]
    for cur_ds in datasets:
        for cur_alg in algos:
            ds_filtered_scores=h_scores.loc[(h_scores['algo']== cur_alg).values & (h_scores['dataset']== cur_ds).values]
            averaged_heterogeneities.append(ds_filtered_scores[ds_filtered_scores["module_0"] !=ds_filtered_scores["module_1"]]['homogeneity_score'].values.mean())
            averaged_homogeneities.append(ds_filtered_scores[ds_filtered_scores["module_0"] == ds_filtered_scores["module_1"]]['homogeneity_score'].values.mean())
            cs.append(algos.index(cur_alg)/float(len(algos)-1))
            txts.append(cur_ds)
            algos_in_scatter.append(cur_alg)

    averaged_homogeneities=np.array(averaged_homogeneities)
    averaged_heterogeneities= np.array(averaged_heterogeneities)
    cs = np.array(cs)
    txts=np.array(txts)
    cs=cs[~np.isnan(averaged_heterogeneities)]
    averaged_homogeneities=averaged_homogeneities[~np.isnan(averaged_heterogeneities)]
    txts=txts[~np.isnan(averaged_heterogeneities)]
    algos_in_scatter=np.unique(np.array(algos_in_scatter)[~np.isnan(averaged_heterogeneities)])   
    averaged_heterogeneities=averaged_heterogeneities[~np.isnan(averaged_heterogeneities)]
    ax.scatter([a if not np.isnan(a) else 10 for a in averaged_heterogeneities], [a if not np.isnan(a) else 0 for a in averaged_homogeneities], c=cs, cmap='jet')
    for i,data in enumerate(zip([a if not np.isnan(a) else 10 for a in averaged_heterogeneities], [a if not np.isnan(a) else 0 for a in averaged_homogeneities])):
        x, y = data
        ax.annotate(txts[i], (x,y))
    ax.set_xlabel("avg_heterogeneity")
    ax.set_ylabel("avg_homogeneity")

    cmap = plt.cm.jet
    colorlist = [ml_colors.rgb2hex(cmap(a / float(np.size(algos_in_scatter)-1))) for a in np.arange(np.size(algos_in_scatter))]
    patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                      markerfacecolor=c) for i, a, c in zip(list(range(len(algos_in_scatter))), algos_in_scatter, colorlist)]
    ax.legend(handles=list(reversed(patches)), loc='lower right', framealpha=0.5)


    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,"h_plot_{}.png".format(prefix)))



