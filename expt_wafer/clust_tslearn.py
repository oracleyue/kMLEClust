# This script loads datasets and run k-shape from tslearn for clustering.

import time
import numpy as np

from tslearn.utils import save_time_series_txt, load_time_series_txt
from tslearn.preprocessing import TimeSeriesScalerMeanVariance
from tslearn.clustering import *

# paths
datapath = "data/"
respath = "results/"

# choose dataset
name = "Wafer"; K = 2

# load datasets
data = load_time_series_txt(datapath + name + '_data.txt')
# labels_gt = np.loadtxt(datapath + name + "_labels.txt", delimiter=' ')

# # kmeans+DTW clustering
# sTime = time.time()
# clust = TimeSeriesKMeans(n_clusters=K, metric="dtw",
#                          max_iter=100, random_state=0).fit(data)
# eTime = time.time() - sTime

# # save results
# np.savetxt(respath + name + '_kMeansDTW' + '_labels.txt', clust.labels_, fmt='%i', delimiter=',')
# print("... kMeansDTW: %.4f sec" % (eTime,))

# # k-Shape clustering
# sTime = time.time()
# X = TimeSeriesScalerMeanVariance(mu=0., std=1.).fit_transform(data)
# clust = KShape(n_clusters=K, n_init=1, random_state=0,
#                max_iter=100).fit(X)
# eTime = time.time() - sTime

# Kernel (GAK) k-means clustering
sTime = time.time()
X = TimeSeriesScalerMeanVariance(mu=0., std=1.).fit_transform(data)
clust = KernelKMeans(n_clusters=2, n_init=20,
                     kernel="gak", kernel_params={"sigma": "auto"},
                     random_state=0).fit(X)
eTime = time.time() - sTime

# save results
np.savetxt(respath + name + '_kGAK' + '_labels.txt', clust.labels_, fmt='%i', delimiter=',')
print("... tslearn: %.4f sec" % (eTime,))
