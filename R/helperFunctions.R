#' Compute Differential modularity score from differential modularity matrix
#'
#' This functions takes the precomputed differential modularity matrix and the genLouvain membership to compute the differential modularity score.
#' @param dwbm differential modularity matrix
#' @param louv.memb louvain community membership
#'
#' @return Vector of differntial modularity score
#'
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

#' get the member vector from alpaca object
#'
#' @param alp alpaca object
#' @param target tf, gene, or all
#'
#' @return member vector
#'
alpaca.getMember = function(alp, target="all") {
  alp2 = alp[[1]]
  memb_vector = as.vector(alp2)
  names(memb_vector) = names(alp2)
  memb_vector=memb_vector[!is.na(memb_vector)]
  if (target =="tf"){
    return(memb_vector[grep("_A",names(memb_vector))])
  }else if (target == "gene"){
    return(memb_vector[grep("_B",names(memb_vector))])
  }else{
    return(memb_vector)
  }
}

#' creates condor object
#'
#' @param elist edge list
#'
#' @return condor object
#' @import condor
condor.createObject = function(elist){
  elist=elist.addTags(elist)
  elist=elist[elist[,3]>0,]
  elist=data.frame(red=elist[,2],blue=elist[,1],weights=elist[,3])
  elist=create.condor.object(elist,return.gcc=T)
  return(elist)
}

#' Run CONDOR clustering
#'
#' @param elist edge list
#' @param qscore TRUE = output qscore / FALSE = do not output qscore
#'
#' @import condor
#' @return condor object
#'
condor.run = function(elist,qscore=F){
  cond.object=condor.createObject(elist)
  cond.object=condor.cluster(cond.object)
  if (qscore){
    cond.object=condor.qscore(cond.object)
  }
  return(cond.object)
}

#' Adds "_A" to first column and "_B" to second column
#'
#' @param elist edge list
#'
#' @return edge list
#'
elist.addTags = function(elist){
  elist[,1]=paste(elist[,1],"A",sep="_")
  elist[,2]=paste(elist[,2],"B",sep="_")
  return(elist)
}

#' check if first two columns are identical
#'
#' @param elist1 edge list
#' @param elist2 edge list
#'
#' @return boolean
#'
elist.isEdgeOrderEqual=function(elist1,elist2){
  check1=identical(as.character(elist1[,1]),as.character(elist2[,1]))
  check2=identical(as.character(elist1[,2]),as.character(elist2[,2]))
  if (check1 & check2){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

#' undo elist.addTags
#'
#' @param elist edge list
#'
#' @return edge list
#'
elist.removeTags = function(elist){
  elist[,1]=gsub("_A","",elist[,1])
  elist[,2]=gsub("_B","",elist[,2])
  return(elist)
}

#' Sorts the edge list based on first two columns in alphabetical order
#'
#' @param elist edge list
#' @return edge list
#'
elist.sort=function(elist){
  A=elistToAdjMat(elist,isBipartite=T)
  ridx=order(rownames(A))
  cidx=order(colnames(A))
  elist=adjMatToElist(A[ridx,cidx])
  return(elist)
}

#' Converts edge list to adjacency matrix
#'
#' @param elist edge list
#' @param isBipartite TRUE = for bipartite / FALSE = for unipartite
#'
#' @import igraph
#' @return Adjcency Matrix
#' @export
elistToAdjMat= function(elist,isBipartite=F){
  if (isTRUE(isBipartite)){
    elist=elist.addTags(elist)
  }
  colnames(elist)=c('from','to','weight')
  A_up = as.matrix(get.adjacency(graph.data.frame(elist),attr = 'weight')) # upper triangular matrix
  A = A_up[rownames(A_up) %in% elist[,1], colnames(A_up) %in% elist[,2]]
  return(A)
}

#' Converts alpaca output into list of data frames
#'
#' @param alp alpaca object
#'
#' @return list of data frames
#'
alpacaObjectToDfList= function(alp){
  alp.out=list()
  if (length(alp)<3){
    tf=alpaca.getMember(alp,target="tf")
    df=data.frame(tf,alp[[2]][names(tf)])
    colnames(df)=c("Membership","Differential Modularity Score")
    rownames(df)=gsub("_A","",rownames(df))
    alp.out[["TF"]]=df
    gene=alpaca.getMember(alp,target="gene")
    df=data.frame(gene,alp[[2]][names(gene)])
    colnames(df)=c("Membership","Differential Modularity Score")
    rownames(df)=gsub("_B","",rownames(df))
    alp.out[["Gene"]]=df
  }else{
    tf=alpaca.getMember(alp,target="tf")
    df=data.frame(tf,alp[[2]][names(tf)],alp[[3]][names(tf)],alp[[4]][names(tf)])
    colnames(df)=c("Membership","Differential Modularity Score","Pvalue","t-statistic")
    rownames(df)=gsub("_A","",rownames(df))
    alp.out[["TF"]]=df
    gene=alpaca.getMember(alp,target="gene")
    df=data.frame(gene,alp[[2]][names(gene)],alp[[3]][names(gene)],alp[[4]][names(gene)])
    colnames(df)=c("Membership","Differential Modularity Score","Pvalue","t-statistic")
    rownames(df)=gsub("_B","",rownames(df))
    alp.out[["Gene"]]=df
  }
  return(alp.out)
}

#' converts adjacency matrix to edge list
#'
#' @param adj_mat adjacency matrix
#'
#' @return edge list
#'
adjMatToElist = function(adj_mat){
  from=rep(row.names(adj_mat), ncol(adj_mat))
  to=rep(colnames(adj_mat), each=nrow(adj_mat))
  weight=as.numeric(unlist(c(adj_mat),use.names = F))
  elist=data.frame(from,to,weight)
  elist=elist.removeTags(elist)
  return(elist)
}

#' CRANE Beta perturbation function. This function will add noice to the node strength sequence.
#'
#' @param nodeD Vector of node strength
#' @param beta beta
#' @param beta_slope TRUE=use predetermined slope to add noise / FALSE = use constant value for noise
#' @return vector with new strength distribution
#'
jutterDegree=function(nodeD,beta,beta_slope=T){
  if (beta==0){
    print("No BETA")
    return(nodeD)
  }
  print(paste("Applying Beta =",beta))
  if (beta<0){
    beta=abs(beta)
    correction=0.06
    const.sig=1
    yint=const.sig*beta
  }else{
    # Correction needed based on number of genes
    correction=0.025
    const.sig=400
    yint=const.sig*beta
  }
  if(beta_slope){
    beta.pos=beta*correction
    beta.neg=beta*correction*0.5
    print(paste("applying beta slope:",beta.pos))
    temp.sigma=nodeD
    temp.sigma[temp.sigma>=0]=abs(temp.sigma[temp.sigma>=0]*beta.pos)+yint
    temp.sigma[temp.sigma<0]=abs(temp.sigma[temp.sigma<0]*beta.neg)+yint
    temp.sigma=abs(temp.sigma)
  }else{
    print(paste("NO Slope Applied, Gene perturbation:", const.sig*beta))
    temp.sigma=const.sig*beta
  }
  nodeD=nodeD+rnorm(length(nodeD),0,temp.sigma)

  return(nodeD)
}

#' Check if data frame is an edge list
#'
#' @param df some data frame
#'
#' @return Boolean
#'
isElist = function(df){
  return(ncol(df)==3 & (is.factor(df[,1])|is.character(df[,1]))&(is.factor(df[,2])|is.character(df[,2]))&is.numeric(df[,3]) )
}

