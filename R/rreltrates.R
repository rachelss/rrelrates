#' Get the distance between 2 nodes / tips.
#' @export
#' @param node1 a tip or node 
#' @param node2 a tip or node 
#' @param t tree containing the nodes 
allnodedist2<-function(node1,node2,t){
  return(dist.nodes(t)[node1,node2])
}

#' Get the relative rates of two decendant branches.
#' @export
#' @param bnode a node in the tree with two decendants
#' @param t the tree
#' @param finalrates ancestral relative rates, which are used as multipliers
#' Method:
#' Get the average distance from bnode to all decendant tips (avgdist).
#' For each of bnode's direct child nodes (c1,c2):
#' Get the total distance from this node to each tip (c1dist,c2dist).
#' Divide this total by avgdist (rate1,rate2).
#' If bnode is not the root node:
#' Multiply each rate by the ancestor relative rate.
#' @references Tamura et al. 2012 PNAS
decendDist<-function(bnode,t,finalrates){
  d<-Descendants(t,bnode,type="tips")[1]
  sumdist<-sum(lapply(d,allnodedist2,bnode,t)[[1]])
  avgdist<-sumdist/length(d[[1]])
  
  c<-Descendants(t,bnode,type="children")
  c1d<-Descendants(t,c[[1]],type="tips")
  c1dist<-sum(lapply(c1d,allnodedist2,bnode,t)[[1]])
  c1distAvg<-c1dist/avgdist
  rate1<-c1distAvg/length(c1d[[1]])
  
  c2d<-Descendants(t,c[[2]],type="tips")
  c2dist<-sum(lapply(c2d,allnodedist2,bnode,t)[[1]])
  c2distAvg<-c2dist/avgdist
  rate2<-c2distAvg/length(c2d[[1]])
  
  if(avgdist==0){
    rate1=0
    rate2=0
  }
  
  a<-Ancestors(t, bnode, type=c("parent"))
  if(a>0){
    r<-paste(a,bnode,sep=",")
    rate1<-rate1*finalrates[r][[1]]
    rate2<-rate2*finalrates[r][[1]]
  }
  rates<-c(rate1,rate2)
  names(rates)<-c(paste(toString(bnode),toString(c[[1]]),sep=","),paste(toString(bnode),toString(c[[2]]),sep=","))
  
  return(rates)
}

#' Checks if value is outside confidence interval.
#' @export
#' @param rate a value of interest
#' @param lb lower bound of confidence interval
#' @param upper bound of confidence interval
outsideCI<-function(rate,lb,ub){
  if(is.numeric(rate[[1]])){
    if(rate[[1]]<lb || rate[[1]]>ub){
      return(rate)
    }
  }
}

#' Gets name of tip.
#' @export
#' @param name tip number
#' @param t tree
rename_branch<-function(name,t){
  branches<-strsplit(name, ",")
  tlabel<-as.numeric(branches[[1]][2])
  return(t$tip.label[tlabel])
}

#' Get the primate tree for a gene from ensembl based on multiz alignment using rensembl.
#' @export
#' @gene name of gene 
#' @plottree specify as false to skip plotting the tree
gettree<-function(gene, plottree=TRUE){
  print(gene)
  t<-tryCatch(primate_tree(gene),error=function(e) NULL)
  if(! is.null(t)){
    if(plottree){
      plot(t)
      nodelabels()
    }
  }
  return(t)
}

#' Adjust a rate
#' 
adjust_br<-function(p,r,gentimes){  
  a<-sapply(names(r),strsplit,",")
  b<-mapply(match,as.character(p),a)
  c<-match(2,b)
  d<-r[[c]]*gentimes[[p]]
  return(d)
}

#' Adjust rates -> gen_time * rel rate.  #TIPS ONLY
#' 
#' Takes a list
#' Returns adjusted list.
adjust_tree<-function(r,gentimes){  #assumes gentimes sp names = tip numbers
  p<-length(gentimes)
#  gentimes<- append(gentimes,rep("NA", length(r)-length(gentimes)))
  r<-sapply(seq(p),adjust_br,r,gentimes)
  return(r)
}

#' Get the relative rates on a tree.
#' @export
#' @param t tree
#' Each rate is for a branch, which is labeled by the two nodes it connects
#' @references Tamura et al. 2012 PNAS
geneinfo<-function(t){
  if(! is.null(t)){
    Z = as.list(t$tip.label)
    X<-t$edge
    X[X[,2]%in%1:length(t$tip),2]<-t$tip[X[X[,2]%in%1:length(t$tip),2]]
    names(t$edge.length)<-paste(X[,1],X[,2],sep=",")
    
    numnodes<-seq(length(t$tip.label)+1,length(t$tip.label)+t$Nnode)
    finalrates<-list()
    for (i in 1:length(numnodes) ) {
      finalrates<-append(finalrates,decendDist(numnodes[i],t,finalrates))
    }
    return(finalrates)
  }
}

#' Get the mean, var, and values outside the confidence interval (mean +- 2 std dev)
#' @export
#' @param finalrates a list of values
#' @param a tree
suminfo<-function(finalrates,t){
  if(! is.null(t)){
    sd<-sqrt(var(as.numeric(finalrates)))
    lb<-mean(as.numeric(finalrates))-2*sd[[1]]
    ub<-mean(as.numeric(finalrates))+2*sd[[1]]
    highrates<-sapply(finalrates,outsideCI,lb,ub,USE.NAMES = TRUE)
    h<-Filter(Negate(is.null), highrates)
    names(h)<-lapply(names(h),rename_branch,t)
    m<-mean(as.numeric(finalrates))
    names(m)<-"mean"
    v<-var(as.numeric(finalrates))
    names(v)<-"var"
    h<-append(h,m)
    h<-append(h,v)
    return(h)
  }
}

