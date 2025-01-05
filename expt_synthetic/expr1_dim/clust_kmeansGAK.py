# This script loads datasets and run k-shape from tslearn for clustering.

import time
import numpy as np
from tslearn.utils import save_time_series_txt, load_time_series_txt
from tslearn.preprocessing import TimeSeriesScalerMeanVariance
from tslearn.clustering import *

import clustMetrics

# paths
datapath = "../data_expr1/"
respath  = "../result_expr1/"

# dataset spec.
mList = (2, 4, 8)
nPara = len(mList)
nExpr = 40
K = 8

# load datasets
for iPara in range(nPara):
    m = mList[iPara];
    print("benchmark for m=%i: start" % m)

    labelM = None   # 2d array, each row is clustering result for a dataset
    for iExpr in range(nExpr):
        fname = datapath + "m" + str(m) + "_dataset" + str(iExpr+1)
        data = load_time_series_txt(fname  + ".txt")
        labels_gt = np.loadtxt(fname + "_labels.txt", delimiter=' ')
        labels_gt[labels_gt == K] = 0

        # kmeans+DTW clustering
        sTime = time.time()
        clust = KernelKMeans(n_clusters=K, kernel="gak",
                             kernel_params={"sigma": "auto"},
                             random_state=0).fit(data)
        eTime = time.time() - sTime
        prec,_ = clustMetrics.evalPrec(clust.labels_, labels_gt, K)
        print("... dataset #%i: %.4f sec, %.4f%%" % (iExpr, eTime, prec*100))

        if labelM is None:
            labelM = clust.labels_
        else:
            labelM = np.vstack((labelM, clust.labels_))

    # export clustering results for "m" testing
    fname = respath + "kMeansGAK_m" + str(m) + "_result.txt"
    np.savetxt(fname, labelM, fmt='%i', delimiter=',')
    print("benchmark for m=%i: done\n" % m)
