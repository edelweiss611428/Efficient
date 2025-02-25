# Clustering Algorithms for Optimising the Average Silhouette Width

### Description
This package implements clustering methods for optimising the Average Silhouette Width (ASW), including:
- The PAMSil algorithm (Van der Laan & Pollard 2003), a k-medoid clustering algorithm for optimising the ASW.
- The Efficient Optimum Silhouette algorithm (effOSil), which performs the exact Optimum Silhoutte clustering (Batool & Hennig 2021, OSil) but O(N) times faster, where N is the number of observations in the dataset.
- The Scalable Optimum Silhouette algorithm (scalOSil), which performs the exact Fast Optimum Silhouette clustering (Batool & Hennig 2021, FOSil) but O(n) times faster, where n is the sub-sample size.
- Automatic selection of the optimal clustering solution.
- Automatic selection of the optimal number of clusters.

![10c4ad28-6fbc-4351-9a8c-00a3ada882ff](https://github.com/user-attachments/assets/0e3b6ecb-1d0d-4793-ab2e-8ef53c5708ae)

### Installation

To install the package, Rtools and Rcpp are required. However, we have been unable to compile the ASW package under the latest version of Rcpp (Rcpp v. 1.0.13). Users should use older Rcpp versions (e.g., Rcpp v.1.0.12) in the meantime.

```
install.packages("devtools")
pkg = "https://cran.r-project.org/src/contrib/Archive/Rcpp/Rcpp_1.0.12.tar.gz"
install.packages(pkg)
devtools::install_github("edelweiss611428/ASW")
```
To check the current Rcpp version, use
```
Rcpp::getRcppVersion()
```

### Contact

To report bugs or seek help with installation or running the package, please contact edelweiss611428@gmail.com.

### References

[Batool, F. and Hennig, C., 2021. Clustering with the average silhouette width. Computational Statistics & Data Analysis, 158, p.107190.](https://www.sciencedirect.com/science/article/abs/pii/S0167947321000244)

[Van der Laan, M., Pollard, K. and Bryan, J., 2003. A new partitioning around medoids algorithm. Journal of Statistical Computation and Simulation, 73(8), pp.575-584.](https://www.tandfonline.com/doi/abs/10.1080/0094965031000136012)

[Batool, F., 2019. Initialization methods for optimum average silhouette width clustering. arXiv preprint arXiv:1910.08644.](https://arxiv.org/abs/1910.08644)

[Rousseeuw, P.J., 1987. Silhouettes: a graphical aid to the interpretation and validation of cluster analysis. Journal of computational and applied mathematics, 20, pp.53-65.](https://www.sciencedirect.com/science/article/pii/0377042787901257)

