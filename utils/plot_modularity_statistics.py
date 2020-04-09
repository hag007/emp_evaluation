import sys

sys.path.insert(0, "../")

import numpy as np
import matplotlib.pyplot as plt

import os
import constants

with open(os.path.join(constants.OUTPUT_GLOBAL_DIR, "modularity_statistics.txt")) as f:
    n_modules = [int(a) for a in f.readline().strip()[1:-1].split(', ')]
    n_large_modules = [int(a) for a in f.readline().strip()[1:-1].split(', ')]
    modularity_scores = [float(a) for a in f.readline().strip()[1:-1].split(', ')]

    prev=0

    for i, cur_score in enumerate(modularity_scores):
        if cur_score<prev:
            stop=i

    plt.subplots(figsize=(20,17))
    plt.cla()
    plt.plot(np.arange(len(n_modules)), n_modules, label="n_modules")
    plt.plot(np.arange(len(n_large_modules)), n_large_modules, label="n_large_modules")
    plt.plot(np.arange(len(modularity_scores)), [a for a in modularity_scores], label="modularity_scores")



    for i, a in enumerate(modularity_scores):
        plt.annotate("{}/{}".format(n_large_modules[i], n_modules[i]), (i,a))

    plt.xlabel("iteration")
    plt.ylabel("modularity score")
    plt.title("Newman-Girvan modularity curve\nin DIP network")
    plt.legend()
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "modularity_statistics.png"))

