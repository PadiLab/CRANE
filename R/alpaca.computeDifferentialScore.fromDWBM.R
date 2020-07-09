#' Compute Differential modularity score from differential modularity matrix
#'
#' This functions takes the precomputed differential modularity matrix and the genLouvain membership to compute the differential modularity score.
#' @param dwbm differential modularity matrix
#' @param louv.memb louvain community membership
#'
#' @return Vector of differntial modularity score
#' @export
#'
#' @examples
alpaca.computeDifferentialScore.fromDWBM = function(dwbm,louv.memb){
  louv.Ascores <- NULL
  louv.Bscores <- NULL
  for (i in 1:max(louv.memb)){
    this.comm <- names(louv.memb)[louv.memb==i]
    this.tfs <- this.comm[grep("_A",this.comm)]
    this.genes <- this.comm[grep("_B",this.comm)]
    if (length(this.tfs)>1){
      tf.sums <- apply(as.matrix(dwbm[this.tfs,this.genes]),1,sum)
      gene.sums <- apply(as.matrix(dwbm[this.tfs,this.genes]),2,sum)
    } else {
      tf.sums <- sum(dwbm[this.tfs,this.genes])
      gene.sums <- dwbm[this.tfs,this.genes]
    }
      this.denom <- sum(dwbm[this.tfs,this.genes])
      louv.Ascores <- c(louv.Ascores,tf.sums/this.denom)
      louv.Bscores <- c(louv.Bscores,gene.sums/this.denom)
  }

  louv.scores <- c(louv.Ascores,louv.Bscores)

  return(louv.scores)
}
