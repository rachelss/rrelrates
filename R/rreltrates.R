require(ape)
require(phytools)
require(phangorn)
require(plotrix)
require(rensembl)

#' Get the distance between 2 nodes / tips.
#' 
#' Takes node1,node2,tree
allnodedist2<-function(node1,node2,t){
  return(dist.nodes(t)[node1,node2])
}

#' Get the relative rates of two decendant branches.
#' 
#' See method for reltime in Tamura et al 2012 PNAS
#' Rate is adjusted by ancestor relative rates.
#' Takes a node in the tree (bnode), the tree (t), and ancestral relative rates, which are used as multipliers.
#' Method:
#' Get the average distance from bnode to all decendant tips (avgdist).
#' For each of bnode's direct child nodes (c1,c2):
#' Get the total distance from this node to each tip (c1dist,c2dist).
#' Divide this total by avgdist (rate1,rate2).
#' If bnode is not the root node:
#' Multiply each rate by the ancestor relative rate.
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
#' 
#' Takes a value, upper bound (mean+2sd), and lower bound (mean-2sd)
outsideCI<-function(rate,lb,ub){
  if(rate[[1]]<lb || rate[[1]]>ub){
    return(rate)
  }
}

#' Gets name of tip.
#' 
#' Takes a tip number and tree
rename_branch<-function(name,t){
  branches<-strsplit(name, ",")
  tlabel<-as.numeric(branches[[1]][2])
  return(t$tip.label[tlabel])
}

#' Get the mean and var of relative rates on a tree.
#' Also returns branches with rates outside the CI
#' 
#' Takes a gene name.
#' Gets the primate tree for this gene using rensembl.
#' Gets all relative rates.
geneinfo<-function(gene){
  print(gene)
  t<-primate_tree(gene)
  plot(t)
  nodelabels()
  
  Z = as.list(t$tip.label)
  X<-t$edge
  X[X[,2]%in%1:length(t$tip),2]<-t$tip[X[X[,2]%in%1:length(t$tip),2]]
  names(t$edge.length)<-paste(X[,1],X[,2],sep=",")
  
  numnodes<-seq(length(t$tip.label)+1,length(t$tip.label)+t$Nnode)
  finalrates<-list()
  for (i in 1:length(numnodes) ) {
    finalrates<-append(finalrates,decendDist(numnodes[i],t,finalrates))
  }
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

