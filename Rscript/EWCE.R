ewce <- function(gene_set,name=NULL,number_of_iteration=99999,gene_set_species=NULL) {
set.seed(52)

#For human - HGNC to ENTREZ
if(gene_set_species=="human"){
  species2entrez <- read.table("/Users/julienbryios/Documents/Data/Mouse2human_Orthologues/Pat_m2h_1to1/mus2human.mapping.all.1to1.txt",header=T)
}  else {
  stop("Please set the species the data is coming from")
}
if(length(name)==0){
  stop("Please provide a name to your gene set")
}

#Get proportion in each cell type
load("/Users/julienbryios/Documents/Data/Projects/Best_tissue_SCZ/Nathan3/Main2/celltype_data_allKImouse_wtHypo_MergedStriatal_1to1only_level1_thresh0_trim0.rda")
proportion <- as.data.frame(celltype_data[[1]]$cell_dists)
colnames(proportion) <- make.names(colnames(proportion))
proportion$musName <- rownames(proportion)

#Keep mouse genes with a 1to1 orthoglog and change mouse Name to Human Name
proportion <- merge(proportion,species2entrez,by="musName")
rownames(proportion) <- proportion$geneName
proportion <- proportion[c(2:(ncol(proportion)-5))]

colnames(gene_set)[1] <- "HGNC.symbol"
proportion_gene_set <- proportion[rownames(proportion)%in%gene_set$HGNC.symbol,]

#Print all genes with a 1 to 1 ortholog in mouse
cat(c("Number of genes with a 1to1 ortholog",nrow(proportion),"\n"))

#Print name of gene set
cat(c("Gene set:\t",name,"\n"))

#Print length of gene set
cat(c("Length of gene set original",nrow(gene_set),"\n"))

#Print number of genes in the gene set and also with a 1to1 ortholog
cat(c("Length of gene set with a 1to1 ortholog",nrow(proportion_gene_set),"\n"))

#Average proportion in each cell type
mean_proportion <- apply(proportion_gene_set,2,mean)

######
#Permutation
######

mean_proportion_bootstrap_df <- matrix(ncol=ncol(proportion_gene_set),nrow=number_of_iteration)
for (i in 1:number_of_iteration){
  mean_proportion_bootstrap_df[i,] <- apply(proportion[sample(nrow(proportion),nrow(proportion_gene_set),replace=F),],2,mean)
  if(i%%10000==0){
  cat(i)
  cat("\n")
  }
}

#Get Pvalue

pvalues <- vector("numeric", ncol(proportion_gene_set))
for (i in 1:ncol(mean_proportion_bootstrap_df)){
  number_null_more_extreme <- length(which(mean_proportion_bootstrap_df[,i] >= mean_proportion[i]))
  pvalues[i] <- (number_null_more_extreme+1)/(nrow(mean_proportion_bootstrap_df)+1)
}
names(pvalues) <- colnames(proportion_gene_set)

#Get Z-score

sd_boot <- apply(mean_proportion_bootstrap_df,2,sd)
mean_boot <- apply(mean_proportion_bootstrap_df,2,mean)
z_scores <- (mean_proportion-mean_boot)/sd_boot

results <- as.data.frame(t(rbind(pvalues,z_scores)))
results <- cbind(cell_type=gsub("\\."," ",rownames(results)),name,results)
rownames(results) <- NULL
return(results)
}

#a <- read.table("/Users/julienbryios/Documents/Data/Gene_sets/Antipsychotic_targets/antipsychotic.txt",header=F)
#test <- ewce(a,name="Antipsychotics",gene_set_species="human")