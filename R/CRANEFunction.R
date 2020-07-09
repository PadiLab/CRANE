library(ALPACA)
library(condor)
library(igraph)
library(foreach)
library(doParallel)


alpaca.computeDifferentialScore.fromDWBM = function(dwbm,louv.memb,isNormal=T){
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
    if (isNormal){
      this.denom <- sum(dwbm[this.tfs,this.genes])
      louv.Ascores <- c(louv.Ascores,tf.sums/this.denom)
      louv.Bscores <- c(louv.Bscores,gene.sums/this.denom)
    }else{
      louv.Ascores <- c(louv.Ascores,tf.sums)
      louv.Bscores <- c(louv.Bscores,gene.sums)
    }
  }

  louv.scores <- c(louv.Ascores,louv.Bscores)


  #   louv.scores[louv.scores<0]=0
  return(louv.scores)
}
alpaca.crane = function(input,alp,a=0.1,b=0,iteration=30,isParallel=F){
  elist1=input[,1:3]
  tfcheck=paste(unique(elist1[,1]),"A",sep="_")
  gcheck=paste(unique(elist1[,2]),"B",sep="_")
  member=alpaca.getMember(alp)

  print("Converting to Adjecency Matrix...")
  A=elistToAdjMat(elist1,isBipartite=T)
  ridx=order(rownames(A))
  cidx=order(colnames(A))
  elist1=adjMatToElist(A[ridx,cidx])

  print("Generating CRANE Networks, this may take a long time (upto 10-20 minutes for Networks with 19000+ nodes) ...")
  if (isParallel){
    cores=detectCores()
    cl <- makeCluster(cores[1]-1)
    registerDoParallel(cl)
    exportFunctions=c("crane.bipartite.fromAdjMatrix","elist.removeTags","adjMatToElist","jutterDegree")
    elist.df=foreach(i=1:iteration,.export = exportFunctions, .combine='cbind', .packages=c("condor","ALPACA"),.verbose = F) %dopar% {
      crane.bipartite.fromAdjMatrix(A,alpha=a,beta=b)[,3]
    }
    stopCluster(cl)
  }else{
    elist.df=NULL
    for (i in 1:iteration){
      elist.df=cbind(elist.df,crane.bipartite.fromAdjMatrix(A,alpha=a,beta=b,gamma=g,beta_slope=T,iterative = T,adjSign =F)[,3])
    }
  }
  gc()

  print("Creating Test Distribution...")
  elist1=elist.addTags(elist1)

  print("Detecting communities in control network...")
  ctrl.pos = elist1[elist1[,3]>=0,1:3]
  ctrl.elist = data.frame(red=ctrl.pos[,2],blue=ctrl.pos[,1],weights=ctrl.pos[,3])
  ctrl.condor = create.condor.object(ctrl.elist)
  ctrl.condor = condor.cluster(ctrl.condor,project=F)
  ctrl.memb = c(ctrl.condor$red.memb[,2],ctrl.condor$blue.memb[,2])
  names(ctrl.memb) = c(as.character(ctrl.condor$red.memb[,1]),as.character(ctrl.condor$blue.memb[,1]))


  print("Comuting Differential Scores...")
  df=NULL
  for (i in 1:iteration){
    net.table=cbind(elist1,elist.df[,i])
    pos.table = net.table[intersect(which(net.table[,3]>=0),which(net.table[,4]>=0)),]
    pos.graph = graph.edgelist(as.matrix(pos.table[,1:2]),directed=T)

    if (length(setdiff(V(pos.graph)$name,names(ctrl.memb)))>0){
      uncounted <- setdiff(V(pos.graph)$name,names(ctrl.memb))
      unc.memb <- sample(1:max(ctrl.memb),length(uncounted),replace=T)
      names(unc.memb) <- uncounted
      ctrl.memb <- c(ctrl.memb,unc.memb)
    }
    if (i %% round(iteration/5)==0){
      print(paste("Computing Differential Score",(i/iteration)*100,"%"))
    }
    dwbm = computeDWBMmat.mscale(pos.table,ctrl.memb[V(pos.graph)$name])

    miss_tf=setdiff(tfcheck,rownames(dwbm))
    miss_g=setdiff(gcheck,colnames(dwbm))
    if (length(miss_tf)>0){
      temp=array(0,dim = c(length(miss_tf),ncol(dwbm)))
      rownames(temp)=miss_tf
      dwbm=rbind(dwbm,temp)
    }
    if (length(miss_g)>0){
      temp=array(0,dim = c(nrow(dwbm),length(miss_g)))
      colnames(temp)=miss_g
      dwbm=cbind(dwbm,temp)
    }

    temp=alpaca.computeDifferentialScore.fromDWBM(dwbm,member)

    if (i>1 & !identical(rownames(df),names(temp))){
      print("error")
      break
    }else{
      df=cbind(df,temp)
    }
  }

  isSig=c()
  tstat=c()
  if (identical(names(alp[[2]]),rownames(df))){
    for (i in 1:nrow(df)){
      if (rownames(df)[i]==""){
        isSig=c(isSig,NA)
        next
      }
      dist=df[i,]
      if (length(dist)>1){
        temp=t.test(dist, mu=alp[[2]][i], alternative = "less",conf.level=0.95)
        isSig=c(isSig,temp$p.value)
        temp2=temp$statistic[[1]]
        names(temp2)=rownames(df)[i]
        tstat=c(tstat,temp2)
      }else{
        temp.v=1
        names(temp.v)=rownames(df)[i]
        isSig=c(isSig,temp.v)
        temp.v[]=10e+200
        tstat=c(tstat,temp.v)
      }
    }
  }
  alp[[3]]=isSig
  alp[[4]]=tstat

  alp.out=alpacaObjectToDfList(alp)

  return(alp.out)
}
alpaca.getMember = function(alp, target="all", remove_na=T) {
  alp2 = alp[[1]]
  memb_vector = as.vector(alp2)
  names(memb_vector) = names(alp2)
  if (isTRUE(remove_na)){
    memb_vector=memb_vector[!is.na(memb_vector)]
  }
  if (target =="tf"){
    return(memb_vector[grep("_A",names(memb_vector))])
  }else if (target == "gene"){
    return(memb_vector[grep("_B",names(memb_vector))])
  }else if (target=="gtf"){
    x=alp.removeTag(memb_vector[grep("_A",names(memb_vector))])
    y=alp.removeTag(memb_vector[grep("_B",names(memb_vector))])
    z=intersect(names(y),names(x))
    rlt=y[z]
    names(rlt)=paste(names(rlt),"B",sep="_")
    return(rlt)
  } else{
    return(memb_vector)
  }
}
condor.createObject = function(elist){
  elist=elist.addTags(elist)
  elist=elist[elist[,3]>0,]
  elist=data.frame(red=elist[,2],blue=elist[,1],weights=elist[,3])
  elist=create.condor.object(elist,return.gcc=T)
  return(elist)
}
condor.run = function(elist,qscore=F){
  cond.object=condor.createObject(elist)
  cond.object=condor.cluster(cond.object)
  if (qscore){
    cond.object=condor.qscore(cond.object)
  }
  return(cond.object)
}
crane.bipartite.fromElist = function(edge_list,alpha=0.1,beta=0,iterative=T,randomStart=F,isPos=F){
  print("Converting Edgelist to Adj Matrix")
  A=elistToAdjMat(edge_list,isBipartite=T)
  if (isPos){
    elist=crane.bipartite.fromAdjMatrix.pos(A,alpha=alpha,beta=beta,randomStart=randomStart)
  }else{
    elist=crane.bipartite.fromAdjMatrix(A,alpha=alpha,beta=beta,randomStart=randomStart)
  }
  return(elist)
}
crane.bipartite.fromAdjMatrix = function(A,alpha=0.1,beta=0,getAdj=F,randomStart=F){
  print(paste("NORMAL PERTURBATION"))
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
    tfD.jit=jutterDegree(tfD,beta,beta_slope = T,tfCount=rown)
    tfD.delta=tfD.jit-tfD
    correction=tfD.delta/coln
    A=A+correction
    print("Beta on Genes")
    A=t(A)
    tfD=rowSums(A)
    tfD.jit=jutterDegree(tfD,beta,beta_slope = F,tfCount=rown)
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
      #Final Update Last Edges
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
crane.unipartite.fromAdjMatrix = function(A,alpha=0.1){
  print(paste("NORMAL PERTURBATION"))
  print(paste("Applying Alpha =",alpha))

  rown=nrow(A)
  coln=ncol(A)

  ridx=sample(1:rown,rown,replace = F)
  A=A[ridx,ridx]


  rlt_A=array(0,dim=dim(A))
  colnames(rlt_A)=colnames(A)
  rownames(rlt_A)=rownames(A)

  A[lower.tri(A)]=0

  rowD=rowSums(A)
  colD=colSums(A)
  gD=colSums(A)

  omin=min(A)
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
      avg=mean(temp)
      std=sd(temp)
      cur_min=min(temp)
      cur_max=max(temp)
      delta=rnorm(n,mean=0,sd=std)
      temp[(i+1):rown]=temp[(i+1):rown]+alpha*delta
      if (any(temp[(i+1):rown]<omin)){
        temp[(i+1):rown]=temp[(i+1):rown]+abs(min(temp[(i+1):rown]))+omax*alpha
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
        if (any(temp[(i+1):rown]<omin)){
          temp[(i+1):rown]=temp[(i+1):rown]+abs(min(temp[(i+1):rown]))+omax*alpha
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
elist.addTags = function(edge_list){
  edge_list[,1]=paste(edge_list[,1],"A",sep="_")
  edge_list[,2]=paste(edge_list[,2],"B",sep="_")
  return(edge_list)
}
elist.isEdgeOrderEqual=function(elist1,elist2){
  check1=identical(as.character(elist1[,1]),as.character(elist2[,1]))
  check2=identical(as.character(elist1[,2]),as.character(elist2[,2]))
  if (check1 & check2){
    return(TRUE)
  }else{
    return(FALSE)
  }
}
elist.removeTags = function(edge_list){
  edge_list[,1]=gsub("_A","",edge_list[,1])
  edge_list[,2]=gsub("_B","",edge_list[,2])
  return(edge_list)
}
elist.sort=function(edge_list){
  A=elistToAdjMat(edge_list,isBipartite=T)
  ridx=order(rownames(A))
  cidx=order(colnames(A))
  elist=adjMatToElist(A[ridx,cidx])
  return(elist)
}
elistToAdjMat= function(edge_list,isBipartite=F){
  if (isTRUE(isBipartite)){
    edge_list=elist.addTags(edge_list)
  }
  colnames(edge_list)=c('from','to','weight')
  A_up = as.matrix(get.adjacency(graph.data.frame(edge_list),attr = 'weight')) # upper triangular matrix
  A = A_up[rownames(A_up) %in% edge_list[,1], colnames(A_up) %in% edge_list[,2]]
  return(A)
}
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
adjMatToElist = function(adj_mat){
  from=rep(row.names(adj_mat), ncol(adj_mat))
  to=rep(colnames(adj_mat), each=nrow(adj_mat))
  weight=as.numeric(unlist(c(adj_mat),use.names = F))
  elist=data.frame(from,to,weight)
  elist=elist.removeTags(elist)
  return(elist)
}
jutterDegree=function(tfD,beta,beta_slope=T,tfCount){
  if (beta==0){
    print("No BETA")
    return(tfD)
  }
  print(paste("Applying Beta =",beta))
  if (beta<0){
    beta=abs(beta)
    # correction=0.000625
    correction=0.06
    # const.sig=10
    const.sig=1
    yint=const.sig*beta
  }else{
    # Correction needed based on number of genes
    correction=0.025
    const.sig=400
    yint=const.sig*beta
    # yint=abs(50*beta)
  }
  if(beta_slope){
    beta.pos=beta*correction
    beta.neg=beta*correction*0.5
    print(paste("applying beta slope:",beta.pos))
    temp.sigma=tfD
    temp.sigma[temp.sigma>=0]=abs(temp.sigma[temp.sigma>=0]*beta.pos)+yint
    temp.sigma[temp.sigma<0]=abs(temp.sigma[temp.sigma<0]*beta.neg)+yint
    temp.sigma=abs(temp.sigma)
  }else{
    print(paste("NO Slope Applied, Gene perturbation:", const.sig*beta))
    temp.sigma=const.sig*beta
  }
  tfD=tfD+rnorm(length(tfD),0,temp.sigma)
  # for (i in 1:length(tfD)){
  #   delta=rnorm(1,0,temp.sigma[i])
  #   tfD[i]=tfD[i]+delta
  # }
  return(tfD)
}





