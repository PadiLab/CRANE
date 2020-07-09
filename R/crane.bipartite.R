#' Pertrubs the bipartite network with fixed node strength
#'
#' @param df Adjacency Matrix or Edge list
#' @param alpha alpha paramter perturbs each edge weights
#' @param beta beta parameter perturbs the strength of each node. Set this to 0 if you want nodes to have node strength identical to the orignal network.
#' @param getAdj TRUE = this will return adjacency matrix instead of edge list
#' @param randomStart FALSE = initialize the first row with completely random edge weights.
#'
#' @return edge list
#' @export
#' @examples
#' \dontrun{
#'
#' # Using Edge list as input
#' elist=crane.bipartite(nonAng)
#' elist=crane.bipartite(nonAng,alpha=0.3)
#'
#' # Using Edge list as input and Adjcency Matrix as output
#' adjMatrix=crane.bipartite(nonAng,alpha=0.1,getAdj=T)
#'
#' # Using Edge list as input and Adjcency Matrix as output
#' A=elistToAdjMat(nonAng)
#' elist=crane.bipartite(A)
#'
#' }

crane.bipartite= function(df,alpha=0.1,beta=0,getAdj=F,randomStart=F){
  if (isElist(df)){
    print("Converting Edgelist to Adj Matrix")
    A=elistToAdjMat(df,isBipartite=T)
  }else{
    A=df
  }

  print(paste("Applying Alpha =",alpha))
  rown=nrow(A)
  coln=ncol(A)

  # randomize nodes
  ridx=sample(1:rown,rown,replace = F)
  cidx=sample(1:coln,coln,replace = F)
  A=A[ridx,cidx]

  rlt_A=array(0,dim=dim(A))
  colnames(rlt_A)=colnames(A)
  rownames(rlt_A)=rownames(A)

  # beta implimentation: beta = 1 models subsample data
  if (beta!=0){
    print("Beta on TFs")
    tfD=rowSums(A)
    tfD.jit=jutterDegree(tfD,beta,beta_slope = T)
    tfD.delta=tfD.jit-tfD
    correction=tfD.delta/coln
    A=A+correction
    print("Beta on Genes")
    A=t(A)
    tfD=rowSums(A)
    tfD.jit=jutterDegree(tfD,beta,beta_slope = F)
    tfD.delta=tfD.jit-tfD
    correction=tfD.delta/rown
    A=t(A)
  }

  tfD=rowSums(A)
  gD=colSums(A)

  omin=min(A)
  omax=max(A)

  if (alpha >0){
    rlt_A[1,]=A[1,]
    if (randomStart){
      rlt_A[1,]=rnorm(length(A[1,]),mean(A[1,]),2)
    }
    print("Constructing Iteratuvely Perturbed Network")
    current_degree=rep(0,coln)
    for (i in 2:rown){
      if (i %% 100 ==0){
        print(i)
      }
      #Randomly add delta to previous the Rows by alpha
      temp=rlt_A[(i-1),]
      n=length(temp)
      avg=mean(temp)
      std=sd(temp)
      cur_min=min(temp)
      cur_max=max(temp)
      delta=rnorm(n,mean=0,sd=std)
      temp=temp+alpha*delta
      temp=(temp/sum(temp))*tfD[(i-1)]
      goal_degree=colSums(A[1:i,])
      k=0
      # Backtrack to Control min max
      while(1){
        if (k >150){
          break
        }
        temp_degree=current_degree+temp
        y=goal_degree-temp_degree

        lidx=which(y<omin)
        ridx=which(y>omax)
        idx=c(lidx,ridx)
        lcor=y[lidx]-omin
        rcor=y[ridx]-omax
        if(length(idx)==0){
          break
        }
        if (length(lidx)>0){
          temp[lidx]=temp[lidx]-abs(lcor)
        }
        if (length(ridx)>0){
          temp[ridx]=temp[ridx]+abs(rcor)
        }
        temp=(temp/sum(temp))*tfD[(i-1)]
        k=k+1
      }
      #Update Final Edges
      rlt_A[(i-1),]=temp
      current_degree=current_degree+temp
      y=goal_degree-current_degree
      rlt_A[i,]=y
    }
  }

  if(!isTRUE(all.equal(rowSums(rlt_A),tfD,check.attributes = F, use.names = F))){
    print("ROW ERROR Alpha limit has been reached")
    return(NULL)
  }
  if(!isTRUE(all.equal(colSums(rlt_A),gD,check.attributes = F, use.names = F))){
    print("COLUMN ERROR Alpha limit has been reached")
    return(NULL)
  }

  print("Sorting Nodes")
  ridx=order(rownames(rlt_A))
  cidx=order(colnames(rlt_A))
  if (getAdj){
    return(rlt_A[ridx,cidx])
  }
  elist=adjMatToElist(rlt_A[ridx,cidx])
  elist=elist.removeTags(elist)
  return(elist)
}
