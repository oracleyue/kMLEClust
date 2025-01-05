# This script loads datasets and run k-shape from tslearn for clustering.

import time
import numpy as np
from tslearn.utils import save_time_series_txt, load_time_series_txt
from tslearn.preprocessing import TimeSeriesScalerMeanVariance
from tslearn.clustering import KShape
from tslearn.clustering import TimeSeriesKMeans

# paths
datapath = "../data_expr3/"
respath  = "../result_expr3/"

# dataset spec.
SNRVec = (0, 2, 4, 8, 10) #(30, 20, 10, 0)
nPara = len(SNRVec)
nExpr = 40
K = 8

# load datasets
for iPara in range(nPara):
    snr = SNRVec[iPara]
    print("benchmark for SNR=%i dB: start" % snr)

    labelM = None   # 2d array, each row is clustering result for a dataset
    for iExpr in range(nExpr):
        fname = datapath + "SNR" + str(snr) + "_dataset" + str(iExpr+1) + ".txt"
        data = load_time_series_txt(fname)

        sTime = time.time()
        # k-shape clustering
        X = TimeSeriesScalerMeanVariance(mu=0., std=1.).fit_transform(data)
        clust = KShape(n_clusters=K, n_init=1, random_state=0,
                    max_iter=100).fit(X)
        # clust.labels_
        eTime = time.time() - sTime
        print("... dataset #%i: %.4f sec" % (iExpr, eTime))

        if labelM is None:
            labelM = clust.labels_
        else:
            labelM = np.vstack((labelM, clust.labels_))

    # export clustering results for "m" testing
    fname = respath + "kShape_SNR" + str(snr) + "_result.txt"
    np.savetxt(fname, labelM, fmt='%i', delimiter=',')
    print("benchmark for SNR=%i dB: done\n" % snr)
