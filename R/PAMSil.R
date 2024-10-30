#' @name PAMSil
#' @title The PAMSil algorithm
#'
#' @description  This function implements the PAMSil algorithm.
#'
#' @usage PAMSil(dx, K)
#'
#' @param dx A dist object, which can be computed using the stats::dist() function.
#' @param K An integer vector specifying the number of clusters. By default, K = 2:12.
#'
#' @return
#' \describe{
#' \item{best_clustering}{The clustering achieving the highest ASW value.}
#' \item{best_asw}{The highest ASW value.}
#' \item{best_medoids}{The medoids associated with the clustering maximizing the ASW.}
#' \item{k}{The estimated number of clusters.}
#' \item{clusterings}{The PAMSil clusterings for all k in K.}
#' \item{asw}{The ASW values associated with the clusterings.}
#' \item{medoids}{The medoids associated with the clustering solutions.}
#' \item{nIter}{The numbers of iterations needed for convergence.}
#' }
#'
#' @details
#' This function implements the PAMSil algorithm proposed by Van der Laan & Pollard (2003),
#' a k-medoid clustering algorithm whose objective function is the ASW.
#'
#' @examples
#' library("cluster")
#' x = scale(faithful)
#' dx = dist(x)
#' PAMSil_clustering = PAMSil(dx = dx, K = 2:12)
#' par(mfrow = c(1,2))
#' plot(faithful, col = PAMSil_clustering$best_clustering, pch = PAMSil_clustering$best_clustering)
#' plot(2:12, PAMSil_clustering$asw, type = "l", xlab = "k", ylab = "ASW")
#' par(mfrow = c(1,1))
#'
#' @references
#' Van der Laan, M., Pollard, K. and Bryan, J., 2003. A new partitioning around medoids algorithm. Journal of Statistical Computation and Simulation, 73(8), pp.575-584.
#'
#' @importFrom cluster pam
#' @importFrom stats dist
#'
#' @author Minh Long Nguyen \email{edelweiss611428@gmail.com}
#' @export

PAMSil = function(dx, K = 2:12){

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

  clusterings = matrix(integer(N*nK), nrow = N)
  nIter = integer(nK)
  asw = numeric(nK)
  medoids = vector("list", length = nK)
  colnames(clusterings) = K
  names(asw) = K
  names(nIter) = K
  names(medoids) = K

  for(i in 1:nK){

    PAM = pam(dx, K[i])
    PAMSilres = .PAMSilCpp(dx, PAM$clustering-1L, PAM$id.med-1L, N, K[i])
    clusterings[,i] = PAMSilres$Clustering
    asw[i] = PAMSilres$ASW
    nIter[i] = PAMSilres$nIter
    medoids[[i]] = PAMSilres$medoids

  }

  idx_max = which.max(asw)
  best_asw = asw[idx_max]
  best_clustering = clusterings[,idx_max]
  k = K[idx_max]
  best_medoids = medoids[[idx_max]]

  return(list(best_clustering = best_clustering, best_asw = best_asw, best_medoids = best_medoids, k = k,
              clusterings = clusterings, asw = asw, medoids = medoids, nIter = nIter))

}





