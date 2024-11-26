#' @name effOSil
#' @title The Efficient Optimum Silhouette algorithm
#'
#' @description  This function implements the Efficient Optimum Silhouette (effOSil) algorithm.
#'
#' @usage effOSil(dx, K, initMethod, variant)
#'
#' @param dx A dist object, which can be computed using the stats::dist() function.
#' @param K An integer vector specifying the number of clusters. By default, K = 2:12.
#' @param initMethod A character vector specifying initialisation methods. By default,
#' initMethod = "average"; however, to achieve the best initialisation in terms of the ASW,
#' various initialisation methods should be used (e.g., initMethod = c("single", "average", "complete", "pam")).
#' See ?Init for more details.
#' @param variant An algorithmic variant. Options include "efficient" and "original". By default, variant = "efficient", indicating that effOSil is used.
#' If variant = "original", the original, computationally expensive OSil algorithm is used.
#'
#' @return
#' \describe{
#' \item{best_clustering}{The clustering achieving the highest ASW value.}
#' \item{best_asw}{The highest ASW value.}
#' \item{k}{The estimated number of clusters.}
#' \item{clusterings}{The effOSil clusterings for all k in K.}
#' \item{asw}{The ASW values associated with the clusterings.}
#' \item{nIter}{The numbers of iterations needed for convergence.}
#' }
#'
#' @details
#' This function implements the Efficient Optimum Silhouette (effOSil) algorithm, an O(N) runtime improvement of
#' the original, computationally expensive Fast OSil (FOSil) algorithm proposed by Batool & Hennig (2021) where N is
#' the number of observations. This function also implements the OSil algorithm for comparision purporses.
#'
#'
#' @examples
#' x = scale(faithful)
#' dx = dist(x)
#' effOSil_clustering = effOSil(dx = dx, K = 2:12)
#' par(mfrow = c(1,2))
#' plot(faithful, col = effOSil_clustering$best_clustering, pch = effOSil_clustering$best_clustering)
#' plot(2:12, effOSil_clustering$asw, type = "l", xlab = "k", ylab = "ASW")
#' par(mfrow = c(1,1))
#'
#' @references
#' Batool, F. and Hennig, C., 2021. Clustering with the average silhouette width. Computational Statistics & Data Analysis, 158, p.107190.
#'
#' @importFrom cluster pam
#' @importFrom stats dist
#'
#' @author Minh Long Nguyen \email{edelweiss611428@gmail.com}
#' @export

effOSil = function(dx, K = 2:12, initMethod = "average", variant = "efficient"){

  if(inherits(dx, "dist") == TRUE){
    N = attr(dx, "Size")
  } else{
    stop("effOSil only inputs a distance matrix of class 'dist'!")
  }

  nK = length(K)

  if((!is.numeric(K)) | (nK == 0)){
    stop("K must be an integer vector!")
  }

  K = as.integer(K)
  nuniqueK = length(unique(K))
  minK = min(K)
  maxK = max(K)

  if(nuniqueK != nK){
    stop("Duplicated number of clusters!")
  } else if(minK <= 1){
    stop("The number of clusters must be larger than 1!")
  } else if(maxK > N){
    stop("The number of clusters cannot be larger than the number of observations!")
  }

  if(length(variant) != 1){
    stop("Only ONE variant could be specified!")
  }

  clusterings = matrix(integer(N*nK), nrow = N)
  nIter = integer(nK)
  asw = numeric(nK)
  colnames(clusterings) = K
  names(asw) = K
  names(nIter) = K

  for(i in 1:nK){
    init = Init(dx, K[i], initMethod)$clustering - 1L

    if(variant == "efficient"){
      OSilres = .effOSilCpp(dx, init, N, K[i])
    } else{
      OSilres = .OSilCpp(dx, init, N, K[i])
    }

    clusterings[,i] = OSilres$Clustering
    asw[i] = OSilres$ASW
    nIter[i] = OSilres$nIter

  }

  idx_max = which.max(asw)
  best_asw = asw[idx_max]
  best_clustering = clusterings[,idx_max]
  k = K[idx_max]

  return(list(best_clustering = best_clustering, best_asw = best_asw, k = k,
              clusterings = clusterings, asw = asw, nIter = nIter))

}



