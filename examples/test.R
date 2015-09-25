
#require(devtools)
#devtools::install_github("dwinter/rensembl")
#devtools::install_github("rachelss/rrelrates")
require(rensembl)
require(rrelrates)

#one gene
tree<-gettree("CCND1")
finalrates<-geneinfo(tree)
suminfo(finalrates,tree)

#multiple genes
genes<-c("CCND1","ATP6")
trees<-lapply(genes,gettree)
finalrates<-sapply(trees,geneinfo)
names(finalrates)<-genes
mapply(suminfo,finalrates,trees)