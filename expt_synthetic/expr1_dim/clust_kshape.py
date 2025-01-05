# This script loads datasets and run k-shape from tslearn for clustering.

import time
import numpy as np
from tslearn.utils import save_time_series_txt, load_time_series_txt
from tslearn.preprocessing import TimeSeriesScalerMeanVariance
from tslearn.clustering import KShape
from tslearn.clustering import TimeSeriesKMeans

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

        sTime = time.time()
        # k-shape clustering
        X = TimeSeriesScalerMeanVariance(mu=0., std=1.).fit_transform(data)
        clust = KShape(n_clusters=K, n_init=1, random_state=0,
                    max_iter=100).fit(X)
        # clust.labels_
        eTime = time.time() - sTime
        prec,_ = clustMetrics.evalPrec(clust.labels_, labels_gt, K)
        print("... dataset #%i: %.4f sec, %.4f%%" % (iExpr, eTime, prec*100))

        if labelM is None:
            labelM = clust.labels_
        else:
            labelM = np.vstack((labelM, clust.labels_))

    # export clustering results for "m" testing
    fname = respath + "kShape_m" + str(m) + "_result.txt"
    np.savetxt(fname, labelM, fmt='%i', delimiter=',')
    print("benchmark for m=%i: done\n" % m)
