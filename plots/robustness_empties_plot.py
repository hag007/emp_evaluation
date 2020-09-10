import matplotlib.pyplot as plt
import pandas as pd
import os
import constants
from matplotlib.lines import Line2D


def main(algos):
    prefixes = ["GE", "PASCAL_SUM"]
    fig, axs = plt.subplots(1, 2, figsize=(20, 10))
    for prefix, ax, cur_algos in zip(prefixes, axs,
                                     [algos, algos]):
        ax.set_facecolor("#FFFFFF")
        ax.grid(color='gray')
        df_means = pd.DataFrame()
        empties_file_format = "robustness_{}_100_{{}}_matrix_empty.tsv".format(prefix)
        for ss_ratio in [0.4, 0.3, 0.2, 0.1]:
            df_means = pd.concat([df_means, pd.read_csv(
                os.path.join(constants.OUTPUT_GLOBAL_DIR, "evaluation", empties_file_format.format(ss_ratio)), sep='\t',
                index_col=0).mean(axis=1)], axis=1)

        df_means = df_means.loc[cur_algos, :]
        for idx, values in df_means.iterrows():
            ax.plot([0.6, 0.7, 0.8, 0.9], values, linestyle=(list(df_means.index).index(idx) * 4, (4, 6)),
                    linewidth=3.0, c=constants.COLORDICT[idx], label=constants.ALGOS_ACRONYM[idx])

        ax.set_title("GE" if prefix is "GE" else "GWAS", fontsize=22)
        ax.set_xlabel("subsample fraction")
        ax.set_ylabel("fraction of non-empty solutions")

        patches = [Line2D([0], [0], linestyle='--', color=constants.COLORDICT[a], label=constants.ALGOS_ACRONYM[a],
                          markersize=12, markerfacecolor=constants.COLORDICT[a], alpha=0.7) for a in algos]
        ax.legend(handles=patches, fontsize=22, loc=(0, 1.1), ncol=3, facecolor='#ffffff')
    plt.tight_layout()


if __name__=="__main__":
    algos=["jactivemodules_greedy", "jactivemodules_sa", "bionet", "netbox", "keypathwayminer_INES_GREEDY", "DOMINO", "hotnet2"]
    main(algos)
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "plots", "figure_19.png"))
