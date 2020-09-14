#' Pertrubs the unipartite network with fixed node strength from adjacency matrix
#'
#' @param A Adjacency Matrix
#' @param alpha alpha paramter perturbs each edge weights
#' @param isSelfLoop TRUE/FALSE if self loop exists. co-expression matrix will have a self-loop of 1. Thus TRUE
#'
#' @return adjacency matrix
#' @export
#'
crane.unipartite = function(A,alpha=0.1,isSelfLoop=F){
  print(paste("Applying Alpha =",alpha))

  rown=nrow(A)
  coln=ncol(A)

  # randomize nodes
  ridx=sample(1:rown,rown,replace = F)
  A=A[ridx,ridx]

  rlt_A=array(0,dim=dim(A))
  colnames(rlt_A)=colnames(A)
  rownames(rlt_A)=rownames(A)

  Aorig=A
  oavg=mean(Aorig)
  diag(Aorig)=mean(Aorig)
  A[lower.tri(A)]=0

  rowD=rowSums(A)
  colD=colSums(A)
  gD=colSums(A)

  omin=min(A)
  oomin=min(A[A>0])
  omax=max(A)
  randn=rown-2
  if (alpha >0){
    rlt_A[1,]=A[1,]
    print("Constructing Iteratuvely Perturbed Network")
    current_degree=rep(0,coln)
    for (i in 2:(rown-1)){
      if (i %% 100 ==0){
        print(i)
      }
      #Randomly add delta to previous the Rows by alpha
      temp=rlt_A[(i-1),]
      n=randn
      avg=mean(temp[(i+1):rown])
      std=sd(Aorig[(i-1),])
      cur_min=min(temp[(i+1):rown])
      cur_max=max(temp[(i+1):rown])
      delta=rnorm(n,mean=0,sd=std)
      temp[(i+1):rown]=temp[(i+1):rown]+alpha*delta
      temp[(i+1):rown][is.na(temp[(i+1):rown])]=oomin
      if (any(temp[(i+1):rown]<omin)){
        temp[(i+1):rown]=temp[(i+1):rown]+abs(min(temp[(i+1):rown]))
      }
      temp[(i+1):rown]=(temp[(i+1):rown]/sum(temp[(i+1):rown]))*(rowD[(i-1)]-sum(temp[1:i]))
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
        temp[(i+1):rown][is.na(temp[(i+1):rown])]=oomin
        if (any(temp[(i+1):rown]<omin)){
          temp[(i+1):rown]=temp[(i+1):rown]+abs(min(temp[(i+1):rown]))+omin
        }
        temp[(i+1):rown]=(temp[(i+1):rown]/sum(temp[(i+1):rown]))*(rowD[(i-1)]-sum(temp[1:i]))
        k=k+1
      }
      #Final Update Last Edges
      rlt_A[(i-1),]=temp
      current_degree=current_degree+temp
      y=goal_degree-current_degree
      rlt_A[i,]=y
      randn=randn-1
    }
  }

  if (isSelfLoop){
    rlt_A[nrow(rlt_A),ncol(rlt_A)]=A[nrow(rlt_A),ncol(rlt_A)]
  }

  # Non-Iterative version of the algorithm (not recommended)
  if(!isTRUE(all.equal(rowSums(rlt_A),rowD,check.attributes = F, use.names = F))){
    print("ROW ERROR Alpha limit has been reached")
    return(NULL)
  }
  if(!isTRUE(all.equal(colSums(rlt_A),colD,check.attributes = F, use.names = F))){
    print("COLUMN ERROR Alpha limit has been reached")
    return(NULL)
  }

  for (i in 1:nrow(rlt_A)){
    for (j in 1:ncol(rlt_A)){
      rlt_A[j,i]=rlt_A[i,j]
    }
  }

  print("Sorting Nodes")
  ridx=order(rownames(rlt_A))
  cidx=order(colnames(rlt_A))

  return(rlt_A[ridx,cidx])
}
