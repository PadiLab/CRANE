#' Find the robust nodes in ALPACA community using CRANE
#'
#' @param input same input for alpaca: first column TF, second column Genes, third column edge weights from baseline condition, fourth column edge weights from disease condition.
#' @param alp alpca object in list format (output from alpaca package)
#' @param alpha alpha paramter perturbs each edge weights
#' @param beta beta parameter perturbs the strength of each node. Set this to 0 if you want nodes to have node strength identical to the orignal network.
#' @param iteration Number of CRANE distributions to create. Higher value leads to better ranking but longer runtime.
#' @param isParallel TRUE = use Multithread / FALSE = do not use Multithread
#'
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import ALPACA
#' @import condor
#'
#' @return list of data frames
#' @export
#' @examples
#' \dontrun{
#'
#' input=cbind(nonAng,ang[,3])
#' alp=alpaca(input,NULL,verbose = F)
#' alpListObject=alpaca.crane(input, alp, isParallel = T)
#'
#' }


alpaca.crane = function(input,alp,alpha=0.1,beta=0,iteration=30,isParallel=F){
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
    exportFunctions=c("crane.bipartite","elist.removeTags","adjMatToElist","jutterDegree")
    elist.df=foreach(i=1:iteration,.export = exportFunctions, .combine='cbind', .packages=c("condor","ALPACA"),.verbose = F) %dopar% {
      crane.bipartite(A,alpha=alpha,beta=beta)[,3]
    }
    stopCluster(cl)
  }else{
    elist.df=NULL
    for (i in 1:iteration){
      elist.df=cbind(elist.df,crane.bipartite(A,alpha=alpha,beta=beta)[,3])
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
