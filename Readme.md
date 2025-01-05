# A k-MLE framework with its k-VARs algorithm for multivariate time-series clustering

## Folders

`/` or `ROOT`: the project root directory.

`expt_synthetic/`: contains the scripts that simulate data and perform
benchmarks that appear in our manuscript.

`expt_wafer/`: contains the scripts that applies our algorithms on the
WAFER dataset, together with the scripts in Python/MATLAB that uses the
third-party libraries in the comparative study.

`toolkits/`: collects the third-party libraries of methods that compared
in our manuscript. The scripts directly under this folder are those the
author wrote as wrapper functions.

## Scripts

Here we list the main scripts/functions that supports the proposed
method. For other scripts for specific experiments, refer to Readme.md file in
each sub-folder.

- `demo.m`: a minimal example to apply k-VARs for time-series clustering.

Scripts in `functions/` for clustering algorithms and simulator:

- `kVARs.m`: implements the k-VARs algorithm for clustering.

- `simVARs.m`: random generation of stable VAR models and simulate to
  gather a collection of multivariate time series for later clustering
  tasks.

- `calcBIC.m`: computes the extended BIC scores to choose model order or the number of clusters, e.g., K, p1,...,pK.

Scripts in `measures/` for clustering performance measures:

- `perfRI.m`: implements Rand Index (RI).

- `perfARI.m`: implements Adjusted Rand Index (ARI).

- `perfNMI.m`: implements Normalised Mutual Information (NMI).

- `perfNID.m`: implements the Normalized Information Distance (NID).

## Experiments in our manuscript

- synthetic-data experiments: read `Readme.md` in the folder `expt_synthetic/`.

- real-data experiments: read `Readme.md` in the folder `expt_wafer/`.


<br/>Last modified on 05 Jan 2025.
