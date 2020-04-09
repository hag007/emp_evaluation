import sys
sys.path.insert(0, "../")

import pandas as pd
import numpy as np
import networkx as nx

BH_METRIC="bh_metric"
BINARY_METRIC="binary_metric"

class MyNetboxScoringSystems:

    def __init__(self, scores, metric=BINARY_METRIC, linker_threshold=0.05):
        self.scores=scores
        self.metric=metric
        self.linker_threshold=linker_threshold

    def calculate_adjusted_scores(self,scores):
        return getattr(self, "{}_{}".format(self.metric,"scores"))(scores)


    def calculate_adjusted_score(self, score):
        return getattr(self, "{}_{}".format(self.metric,"score"))(score)

    def get_threshold(self):
        return getattr(self, "get_{}_threshold".format(self.metric))()


    def binary_metric_scores(self, scores):
        scores["adjusted_score"]=scores["score"]
        return scores

    def linker_threshold(self):
        return self.linker_threshold



    def binary_metric_score(self, score):
        return score

    def get_binary_metric_threshold(self):
        return 1

