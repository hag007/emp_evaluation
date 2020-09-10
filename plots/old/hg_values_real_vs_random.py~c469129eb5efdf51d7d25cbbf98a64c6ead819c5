import sys
sys.path.insert(0, '../')
import utils.add_GO_terms_metadata_agg
import pandas as pd
import numpy as np
import os
import constants
import argparse
import utils.go_hierarcies as go_hierarcies
from statsmodels.sandbox.stats.multicomp import fdrcorrection0


dict_result, go2geneids, geneids2go, entrez2ensembl = go_hierarcies.build_hierarcy(
    roots=['GO:0008150'])
vertices = dict_result.values()[0]['vertices']

terms_to_genes={}

def get_all_genes_for_term(vertices, cur_root, term, in_subtree):
    if term in terms_to_genes:
        return terms_to_genes[term]
    in_subtree = in_subtree or term == cur_root
    all_genes = set()
    if in_subtree:
        all_genes.update(go2geneids[cur_root])

    for cur_child in vertices[cur_root]["obj"].children:
        all_genes.update(get_all_genes_for_term(vertices, cur_child.id, term, in_subtree))

    terms_to_genes[term]=all_genes
    return all_genes

kv={}
def get_permuted_values(row):
    kv[row.loc["GO id"]]=np.array(row.loc["dist_n_samples"][1:-1].split(", "), dtype=np.float32)

def main(dataset="SOC", algo="jactivemodules_sa", n_permutations=300,
         csv_file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr",
                                    "MAX/emp_diff_modules_{dataset}_{algo}.tsv"), prefix="GE"):
    global kv
    csv_file_name = csv_file_name.format(dataset=dataset, algo=algo)
    df = pd.read_csv(csv_file_name, sep='\t', index_col=0)
    df["GO id"]=df.index
    df = df.rename(columns={"filtered_pval": "hg_pval"})

    n_genes = [len(get_all_genes_for_term(vertices, cur_go_id, cur_go_id, cur_go_id == cur_go_id)) for i, cur_go_id in
               enumerate(df.index.values)]
    df["n_genes"] = pd.Series(n_genes, index=df.index)

    depth = [dict_result.values()[0]['vertices'][cur_go_id]['D'] for i, cur_go_id in enumerate(df.index.values)]
    df["depth"] = pd.Series(depth, index=df.index)

    kv={}
    df.loc[:,["GO id", "dist_n_samples"]].apply(get_permuted_values, axis=1)
    for i in np.arange(n_permutations):
        print i
        df["hg_pval_random_{}".format(i)]=[kv[a][i] if len(kv[a])>i else 0 for a in df.index]

    df["real_hg_pval"]=df["hg_pval"].apply(lambda a: a if type(a)!=str else max(np.array(a[1:-1].split(", "), dtype=np.float32)))

    df=df.drop(["hg_pval", "dist_n_samples"], axis=1)
    df.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"hg_real_random_{}_{}.tsv".format(dataset,algo)), sep='\t')



if __name__=="__main__":


    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--datasets', dest='datasets', default="SHERA")
    parser.add_argument('--prefix', dest='prefix', default="GE")
    parser.add_argument('--algos', dest='algos', default="jactivemodules_greedy")
    parser.add_argument('--n_permutations', dest='n_permutations', default="5")

    args = parser.parse_args()

    datasets=args.datasets.split(",")
    algos=args.algos.split(",")
    prefix = args.prefix
    n_permutations=int(args.n_permutations)

    n_terms = pd.DataFrame(index=algos, columns=datasets)
    hg_cutoffs = pd.DataFrame(index=algos, columns=datasets)
    emp_cutoffs = pd.DataFrame(index=algos, columns=datasets)
    for cur_ds in datasets:
        for cur_alg in algos:

            main(dataset=cur_ds, algo=cur_alg, n_permutations=n_permutations,
         csv_file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR, "emp_fdr",
                                    "MAX/emp_diff_modules_{dataset}_{algo}.tsv"), prefix="GE")





