library(rrelrates)

a<-"(a:0.1,b:0.2);"
tree<-read.tree(a)
geneinfo(tree)

############

genes<-c("CCND1","ATP6")
trees<-lapply(genes,gettree)
finalrates<-sapply(trees,geneinfo)
names(finalrates)<-genes
mapply(suminfo,finalrates,trees,SIMPLIFY = FALSE)

#############

gos<-get_go_ids("mammary")
genenames<-get_gene_names(gos)

trees<-lapply(genenames[[1]],gettree)
finalrates<-lapply(trees,geneinfo)
names(finalrates)<-genenames[[1]]
results<-mapply(suminfo,finalrates,trees,SIMPLIFY = FALSE)
unlist(results)
results

last4 <- function(x){
  substr(x, nchar(x)-3, nchar(x))
}

spp<-unlist(lapply(lapply(results, `[`, 1),names))
names(spp)<-NULL
remove<-c("mean","var",NA)
extra<-setdiff(spp, remove)
table(unlist(lapply(extra,last4)))
table(unlist(lapply(unlist(lapply(trees, function(x) x$tip.label)),last4)))

################
primates <- c('papio_anubis', 'homo_sapiens', 'pongo_abelii' , 'pan_troglodytes', 'chlorocebus_sabaeus', 'gorilla_gorilla', 'macaca_mulatta', 'callithrix_jacchus')
seqs4<-get_seqs("HOXA9",primates)

#realign
seqs5<-as.DNAbin.list(seqs4)
seqs_aligned<-mafft(seqs5)
seqs_align_char=sapply(seqs_aligned, paste, collapse="")

#build phylogeny
p<-floor(runif(1, min=0, max=10001))
x<-floor(runif(1, min=0, max=10001))
rax_tree<-raxml(seqs_aligned, m = "GTRGAMMA", f = "a", N = 10, p = p, x = x, exec = '/Users/rachelschwartz/partitionfinder/programs/raxml')
plot.phylo(root(rax_tree$bestTree,"callithrix_jacchus"))

finalrates<-geneinfo(drop.tip(root(rax_tree$bestTree,"callithrix_jacchus"),"callithrix_jacchus"))
suminfo(finalrates,root(rax_tree$bestTree,"callithrix_jacchus"))

###################
t<-gettree("AREG")
finalrates<-geneinfo(t)
suminfo(finalrates,t)