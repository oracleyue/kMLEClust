# k-VARs: simulation and benchmark

## Folders

`expr1_dim`: generating data for *Experiment 1*, which benchmarks
our algorithms against the baselines for the scalability of clustering.

`expr2_BIC`: uses the extended BIC method to choose model order and the
number of clusters.

`expr3_SNR`: compares the performance of methods for time series generated under different SNRs.

`data_expr1`, `data_expr3`: keep data for the corresponding experiments.

`result_expr1`, `result_expr2`, `result_expr3`: save results for the corresponding experiments.


## Scripts

**Experiment 1:**

- `expr1_datagen.m`: generates synthetic datasets for benchmark.

- `clust_kVARs.m`: runs Vector k-ARs or MVAR on the synthetic
  dataset.

- `clust_kshape.py`: runs k-Shape method from `tslearn`.

- `clust_kmeansDTW.py`: runs k-DBA method (k-means + DTW barycenter
  average) from `tslearn`.

- `clust_kmeansGAK.py`: runs kernel k-means method from `tslearn`.

- `clust_kSC.m`: runs k-SC method from the third-party library in `ROOT/toolkits/k-SC/`.

- `expr1_boxplot.m`: after all methods being applied and generating
  results, it runs to draw the boxplot of clustering performance given
  in our manuscript.

(Misc.: `clustMetrics.py` is a script which is used in Python files to temporarily check the precision.)

**Experiment 2:**

`expr2_BIC_heatmap.m` runs the extended BIC method to select model order and the number of clusters, which exports a heatmap figure in your manuscript.

**Experiment 3:**

- `expr3_datagen.m`: generates synthetic datasets for benchmark.

- `clust_kVARs.m`: runs k-VARs clustering for Experiment 3.

- `clust_kshape.py`: runs k-Shape method from `tslearn` clustering for Experiment 3.

- `clust_kmeansDTW.py`: runs k-DBA method (k-means + DTW barycenter
  average) from `tslearn` clustering for Experiment 3.

- `clust_kmeansGAK.py`: runs kernel k-means method from `tslearn` clustering for Experiment 3.

- `clust_kSC.m`: runs k-SC method from the third-party library in `ROOT/toolkits/k-SC/` clustering for Experiment 3.

- `expr3_boxplot.m`: after all methods being applied and generating
  results, it runs to draw the boxplot of clustering performance given
  in our manuscript.

- `simVARs_SNRs.m`: the additional VAR simulator that generates the synthetic datasets with given SNRs. (Misc.: `demo_SNRs.m` is a script which tests/demonstrates how this simulator can be applied.)


## Experiment reproduction

To reproduce the results in our manuscript, in each experiment, you may need to first run the corresponding `expr?_datagen.m` script to generate data, then apply each clustering script with a name prefixed by `clust_????.m/py`, and lastly run `expr?_boxplot.m` to process the results and generate images.

## PDF files

You may or may not see PDF files that are the raw figures, and which should be self-explained by in our manuscript.


<br/>Last modified on 04 Jun 2021.
