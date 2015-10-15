#' Get GO id's using biomart - assumes human
#' @export
#' @param keyword term of interest e.g. "mammary"
#' @return list of GO id's related to keyword
get_go_ids<-function(keyword){
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
#  filters = listFilters(ensembl)
#  attributes = listAttributes(ensembl)
  
  goterms<-Term(GOTERM)
  sp_go<-grep(pattern = keyword,x = goterms,ignore.case = TRUE,value=TRUE)
  sp_go_id<-names(sp_go)
  return(sp_go_id)
}

#' Get gene names given GO id's
#' @export
#' @param go_ids vector of GO id - possibly from get_go_ids
#' @return list of genes names
get_gene_names<-function(go_ids){
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  results <- getBM(attributes = c("external_gene_name"), filters = c("go_id"), values = go_ids, mart = ensembl)
  return(results)
}

#' Get a set of sequences for a gene from the multiz alignment given a list of species
#' sequences are in a list with names = species, sequence as a vector of sites
#' @param symbol gene symbol
#' @param sp_list list of species
#' @export
get_seqs<-function(symbol,sp_list){
  symbol_info<-lookup_symbol(symbol)
  reg<-paste(symbol_info$seq_region_name,":",symbol_info$start,"-",symbol_info$end,sep="")
  a<-alignment_region(reg,species="homo_sapiens")
  seqs<-lapply(seq(length(a[[1]]$alignments)),function(x) a[[1]]$alignments[[x]]$seq)
  a_species<-lapply(seq(length(a[[1]]$alignments)),function(x) a[[1]]$alignments[[x]]$species)
  seqs2<-sapply(seqs,tolower)
  seqs3<-sapply(seqs2,strsplit,"")
  names(seqs3)<-unlist(a_species)
  which_species<-lapply(names(seqs3),is.element,sp_list)
  seqs4<-seqs3[unlist(which_species)]
  return(seqs4)
}
