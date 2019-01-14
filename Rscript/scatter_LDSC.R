scatter_ldsc <- function(f,errorbar=TRUE,method="Bonferroni",plot_extended=FALSE) {
library(ggplot2)
library(ggrepel)

GWAS <- read.table(f, header = T)
GWAS$Category <- gsub("L2_0","",as.character(GWAS$Category))
GWAS$Category <- gsub("_"," ",GWAS$Category)
GWAS$Category <- gsub("\\."," ",GWAS$Category)
GWAS$Category <- gsub("bed","",GWAS$Category)
GWAS$Category <- gsub(" $","",GWAS$Category)

GWAS <- GWAS[!GWAS$Category == "base",]

#Get -log10p
GWAS$log10p <- -log10(GWAS$Enrichment_p)

if(plot_extended==FALSE){
  GWAS <- GWAS[!grepl("500",GWAS$Category),]
}

if(method == "Bonferroni") {
  threshold = 0.05/nrow(GWAS)
  file_name = paste0(f,".scatter.Bonf.pdf")
}

if(method == "5FDR") {
  threshold = max(GWAS$Enrichment_p[p.adjust(GWAS$Enrichment_p,method="fdr")<0.05])
  file_name = paste0(f,".scatter.5FDR.pdf")
}

#Add significant labels only
lab_to_add <- GWAS$Category[GWAS$Enrichment_p <= threshold]

GWAS <- GWAS[order(GWAS$Enrichment),]
GWAS$Category <- factor(GWAS$Category, levels = GWAS$Category)

limits <- aes(xmax = GWAS$Enrichment + GWAS$Enrichment_std_error, xmin = GWAS$Enrichment - GWAS$Enrichment_std_error)
p1 <- ggplot(GWAS, aes(x = Enrichment, y = log10p, colour = Category, label = Category)) + geom_point()
p1 <- p1  + geom_hline(yintercept = -log10(threshold),linetype = "longdash")
p1 <- p1 + geom_vline(xintercept = 1,linetype = "longdash")
p1 <- p1 + geom_label_repel(aes(label = ifelse(GWAS$Category %in% lab_to_add,as.character(GWAS$Category),'')), size = 3) 
p1 <- p1 + theme_classic() + theme(legend.position = "none")
p1 <- p1 + ylab(expression('-log'[10]*'(pvalue)')) + xlab(expression('h'^2*' Enrichment')) 
if (errorbar == TRUE){
  p1 <- p1 + geom_errorbarh(limits) 
}
ggsave(p1,filename = file_name)
}

#Run function
#scatter_LDSC("LDSC.results")
#scatter_ldsc("/Users/julienbryios/Documents/Data/Projects/AN_OCD_meta/Data/Partition_h2/AN_OCD_Meta_LDSC_ATAC_Neanderthal.results")
#scatter_ldsc("/Users/julienbryios/Documents/Data/Projects/ATAC-seq_acc/SCZ2_ATACacc_int.results",FALSE)
#scatter_ldsc("/Users/julienbryios/Documents/Data/Projects/STARseq/SCZ2_WGSS_int.results",FALSE)
#scatter_ldsc("/Users/julienbryios/Documents/Data/Projects/SCZ3/v1_5_10_2017/partitioned_h2/SCZ3_ATAC_seq_DLPFC.bed_dir.results",FALSE)
#scatter_ldsc("/Users/julienbryios/Documents/Data/Projects/SCZ3/v1_5_10_2017/partitioned_h2/SCZ3_ATAC_seq_DLPFC.bed_dir.results",FALSE,"5FDR")
#scatter_ldsc("/Users/julienbryios/Documents/Data/Projects/AN_OCD_meta/Results/Partition_h2/AN_OCD_Meta_LDSC_53.results",FALSE,"5FDR")
#scatter_ldsc("/Users/julienbryios/Documents/Data/Projects/SCZ2_asian/Report/LDSC_output/SCZ_EAS.ATAC.Neanderthal.results",FALSE,"5FDR",plot_extended=FALSE)
#scatter_ldsc("/Users/julienbryios/Documents/Data/Projects/SCZ2_asian/Report/LDSC_output/SCZ_EUR.ATAC.Neanderthal.results",FALSE,"5FDR",plot_extended=FALSE)
#scatter_ldsc("/Users/julienbryios/Documents/Data/Projects/MDD/Hardac_plot_final/MDD_no_int.txt",TRUE,"Bonferroni",plot_extended=TRUE)
#scatter_ldsc("/Users/julienbryios/Documents/Data/Projects/ATAC-seq/Revision/LD_score/SCZ2_ATAC_300bp.bed2_dir2.results",FALSE,"Bonferroni",plot_extended=FALSE)
scatter_ldsc("/Users/julienbryios/Documents/Data/Projects/SCZ3/Functional_Genomics/v2_12_2_2018/partitioned_h2/Results/SCZ3_Finucane_no_annot.results",FALSE,"Bonferroni",plot_extended=FALSE)
