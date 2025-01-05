# k-VARs clustering for WAFER dataset

Data Name: `WAFER`

Source: `http://www. mustafabaydogan.com`

## Folders

`.` (current folder): keeps all scripts to reproduce results in our NIPS submission.

`data/`: contains the original WAFER dataset, and re-formated data files for `tslearn` tools.

`results/`: where to save the processing results.

`data_plots/`: contains scripts to draw figures for Wafer dataset.

## Scripts

**Clustering**:

- `clust_kVARs.m`: applies k-VARs method;

- `clust_kSC.m`: applies k-SC method;

- `clust_tslearn.py`: runs methods (`kmeans+DTW`, `k-Shape`, `KernelKMeans`) from `tslearn` toolbox.

**Preparing data**:

- `data_summary.m`: loads the dataset `Wafer.mat` and print is basic information;

- `data_mat2txt.m`: converts `wafer.mat` into the `txt` format for `tslearn`, which includes `Wafer_data.txt` and `Wafer_labels.txt` saved in `data/`.

And basic statistics for Wafer data:

- `ts_plots.m`: plots time series from WAFER;

- `acf_plots.m`: generates ACF/XCF plots, including all `acf_plot_*.pdf*` files;

- `pacf_plots.m`: generates PACF plots, including `pacf_plots.pdf` file.


<br/>Last modified on 05 Jan 2025.
