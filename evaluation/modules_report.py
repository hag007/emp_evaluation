import sys

sys.path.insert(0, '../')

import pandas as pd

from fastsemsim.SemSim import *

import constants
from statsmodels.sandbox.stats.multicomp import fdrcorrection0


def summarize_modules_ehr(prefix, datasets, algos, base_folder,
                          terms_file_name_format="emp_diff_modules_{}_{}_passed_oob.tsv"):
    df_statistics = pd.DataFrame()
    df_full_data = pd.DataFrame()
    for cur_ds in datasets:
        for cur_alg in algos:
            terms_file_name = os.path.join(base_folder, terms_file_name_format.format(cur_ds, cur_alg))
            modules_file_name = os.path.join(constants.TRUE_SOLUTIONS_DIR, "{}_{}".format(cur_ds, cur_alg), "report",
                                             "modules_summary.tsv")

            res = modules_ehr_for_solution(cur_alg, cur_ds, prefix, terms_file_name=terms_file_name,
                                           modules_file_name=modules_file_name)
            if res is None: continue
            tps, fps, sig_hg_genes, sig_emp_genes, statistics, full_data = res
            statistics["algo"] = cur_alg
            statistics["dataset"] = cur_ds
            statistics["id"] = "{}_{}".format(cur_alg, cur_ds)
            df_statistics = df_statistics.append(statistics, ignore_index=True)

            full_data["algo"] = cur_alg
            full_data["dataset"] = cur_ds
            df_full_data = df_full_data.append(full_data)

    df_statistics = df_statistics.set_index("id")
    df_statistics.to_csv(
        os.path.join(constants.OUTPUT_GLOBAL_DIR, "module_cache_files", "modules_statistics_{}.tsv".format(prefix)),
        sep='\t')
    df_full_data.to_csv(
        os.path.join(constants.OUTPUT_GLOBAL_DIR, "module_cache_files", "modules_full_data_{}.tsv".format(prefix)),
        sep='\t')


def modules_ehr_for_solution(algo_sample=None, dataset_sample=None, prefix=None,
                             terms_file_name=os.path.join(constants.OUTPUT_GLOBAL_DIR, "oob", "emp_diff_{}_{}_oob.tsv"),
                             modules_file_name=os.path.join(constants.TRUE_SOLUTIONS_DIR,
                                                            "{}_{}/report/modules_summary.tsv"), emp_ratio_th=0.5):
    try:
        output_terms = pd.read_csv(terms_file_name.format(dataset_sample, algo_sample), sep='\t', index_col=0).dropna()
    except Exception as e:
        print("error: {}".format(e))
        return None

    output_terms = output_terms.sort_values(["hg_pval_max"], ascending=False).sort_values(["emp_pval_max"])
    output_terms = output_terms.rename(columns={"filtered_pval": "hg_pval_max"})
    filtered_genes = output_terms.loc[
        np.logical_and.reduce([output_terms["n_genes"].values > 5, output_terms["n_genes"].values < 500]), ["GO name",
                                                                                                            "hg_pval_max",
                                                                                                            "emp_pval_max",
                                                                                                            "passed_oob_permutation_test"]]

    try:
        output_modules = pd.read_csv(modules_file_name.format(dataset_sample, algo_sample), sep='\t',
                                     index_col=0).dropna()
    except Exception as e:
        print("error: {}".format(e))
        return None

    statistics = {}
    full_data = pd.DataFrame()

    ## correct HG scores
    sorted_genes_hg = filtered_genes.sort_values(by=['hg_pval_max'], ascending=False)
    sig_genes_hg_pval = np.append(sorted_genes_hg["hg_pval_max"].values,
                                  np.zeros(7435 - np.size(sorted_genes_hg["hg_pval_max"].values)))
    sig_genes_hg_pval = [10 ** (-x) for x in sig_genes_hg_pval]
    fdr_results = fdrcorrection0(sig_genes_hg_pval, alpha=0.05, method='indep', is_sorted=False)
    n_hg_true = len([cur for cur in fdr_results[0] if cur])
    sig_hg_genes = sorted_genes_hg.iloc[:n_hg_true, :] if n_hg_true > 0 else 0

    ## summarize HG score
    hg_cutoff = 0
    statistics["hg_cutoff"] = 0
    statistics["hg_corrected"] = 0
    if type(sig_hg_genes) != int:
        hg_cutoff = sig_hg_genes.iloc[- 1]["hg_pval_max"]
        statistics["hg_cutoff"] = round(hg_cutoff, 2)
        statistics["hg_corrected"] = len(sig_hg_genes.index)
        print "HG cutoff: {}, n={}".format(hg_cutoff, len(sig_hg_genes.index))

    ## correct EMP score
    sorted_genes_emp = filtered_genes.sort_values(by=['emp_pval_max'])
    sorted_genes_emp.loc[sorted_genes_emp['emp_pval_max'] == 0, 'emp_pval_max'] = 1.0 / 5000
    sig_genes_emp_pval = sorted_genes_emp["emp_pval_max"].values
    fdr_results = fdrcorrection0(sig_genes_emp_pval, alpha=0.05, method='indep', is_sorted=False)
    n_emp_true = len([cur for cur in fdr_results[0] if cur])
    sig_emp_genes = sorted_genes_emp.iloc[:n_emp_true, :]
    emp_cutoff = sig_emp_genes.iloc[- 1]["emp_pval_max"] if n_emp_true > 0 else 0

    ## summarize EMP scores
    statistics["emp_cutoff"] = emp_cutoff
    statistics["emp_corrected"] = len(sig_emp_genes.index)
    print "EMP cutoff: {}, n={}".format(emp_cutoff, len(sig_emp_genes.index))

    ## extract number of modules
    n_modules = 0
    if output_terms.shape[0] != 0:
        n_modules = output_terms.iloc[0]["hg_pval"].count(',') + 1
    statistics["n_modules"] = n_modules

    ## identify true positive and false positive terms
    tps = {}
    fps = {}
    for a in range(n_modules):
        tps[a] = []
        fps[a] = []
    for go_id, cur in output_terms.iterrows():
        hgs = list(np.array(cur["hg_pval"][1:-1].split(", "), dtype=np.float32))
        emps = list(np.array(cur["emp_pval"][1:-1].split(), dtype=np.float32))
        for i, v in enumerate(zip(hgs, emps)):
            hg, emp = v
            if hg >= hg_cutoff and emp <= emp_cutoff:
                tps[i].append("{}: {}".format(go_id, cur["GO name"]))
            elif hg >= hg_cutoff and emp > emp_cutoff:
                fps[i].append("{}: {}".format(go_id, cur["GO name"]))

    ## summarize figures per module

    real_modules_counter = 0
    for a in range(n_modules):
        statistics["module_{}_emp_ratio".format(a)] = round(float(len(tps[a])) / max(len(tps[a]) + len(fps[a]), 1), 2)
        statistics["module_{}_tp".format(a)] = len(tps[a])
        statistics["module_{}_fp".format(a)] = len(fps[a])
        statistics["module_{}_total".format(a)] = len(tps[a]) + len(fps[a])
        statistics["module_{}_size".format(a)] = output_modules.loc[a, '#_genes']

        full_data.loc["{}_{}_module_{}".format(dataset_sample, algo_sample, a), "EHR"] = round(
            float(len(tps[a])) / max(len(tps[a]) + len(fps[a]), 1), 2)
        full_data.loc["{}_{}_module_{}".format(dataset_sample, algo_sample, a), "tp"] = "\n".join(tps[a])
        full_data.loc["{}_{}_module_{}".format(dataset_sample, algo_sample, a), "fp"] = "\n".join(fps[a])
        full_data.loc["{}_{}_module_{}".format(dataset_sample, algo_sample, a), "module_size".format(a)] = \
        output_modules.loc[a, '#_genes']

        if statistics["module_{}_emp_ratio".format(a)] > emp_ratio_th:
            real_modules_counter += 1

    full_data.sort_index(inplace=True)
    statistics["real_modules_ratio"] = round(float(real_modules_counter) / max(n_modules, 1), 2)

    file_name = os.path.join(constants.OUTPUT_GLOBAL_DIR, "module_cache_files",
                             "full_modules_report_{}_{}.tsv".format(algo_sample, dataset_sample))
    full_data.to_csv(file_name, sep='\t')

    return tps, fps, sig_hg_genes, sig_emp_genes, statistics, full_data
