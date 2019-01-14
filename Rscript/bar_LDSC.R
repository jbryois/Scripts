bar_ldsc <- function(f,plot_extended=FALSE) {
  library(ggplot2)
  library(gridExtra)
  library(grid)

  ldsc_results <- f
  GWAS <- read.table(ldsc_results,header  =  T)
  GWAS$Category <- gsub("L2_0","",as.character(GWAS$Category))
  GWAS$Category <- gsub("_"," ",GWAS$Category)
  GWAS$Category <- gsub("\\."," ",GWAS$Category)
  GWAS$Category <- gsub("bed","",GWAS$Category)
  GWAS$Category <- gsub(" $","",GWAS$Category)
  GWAS$Category <- gsub(" 2$","",GWAS$Category)
  
  GWAS$log10p <- -log10(GWAS$Enrichment_p)

  GWAS <- GWAS[!GWAS$Category  == "base",]

  if(plot_extended==FALSE){
    GWAS <- GWAS[!grepl("500",GWAS$Category),]
  }
  
  GWAS <- GWAS[order(GWAS$Enrichment),]
  GWAS$Category <- factor(GWAS$Category,levels  =  GWAS$Category)

  #5%FDR threshold
  #threshold = max(GWAS$Enrichment_p[p.adjust(GWAS$Enrichment_p,method="fdr")<=0.05])
  
  dodge <- position_dodge(width = 0.9)
  limits <- aes(ymax  =  GWAS$Enrichment + GWAS$Enrichment_std_error, ymin = GWAS$Enrichment - GWAS$Enrichment_std_error)
  p1 <- ggplot(GWAS, aes(y = Enrichment, x = Category,fill=Category)) 
  p1 <- p1 + geom_bar(position = "dodge", stat = "identity", color = "black") +xlab("")
  p1 <- p1 + geom_errorbar(limits, position = dodge, width = 0.25) + geom_hline(yintercept = 1) + coord_flip()
  p1 <- p1 + theme_bw() + theme(axis.text.y = element_text(angle = 0,hjust = 1,vjust = 0.5,size = 8),legend.position = "none") 
  p2 <- ggplot(GWAS, aes(y = log10p, x = Category,fill=Category)) + coord_flip() 
  p2 <- p2 + geom_bar(position = "dodge", stat = "identity", color = "black") 
  p2 <- p2 + ylab(expression('-log'[10]*'(pvalue)')) 
  p2 <- p2 + geom_hline(yintercept = -log10(0.05/nrow(GWAS)))
  #p2 <- p2 + geom_hline(yintercept = -log10(threshold),linetype=1)
  p2 <- p2 + theme_bw()
  p2 <- p2 + theme(axis.text.y = element_text(angle = 0,hjust = 1,vjust = 0.5,size = 8),legend.position="none")
  #p2 <- p2 + theme(axis.title.y=element_blank(),
  #                 axis.text.y=element_blank(),
  #                 axis.ticks.y=element_blank(),legend.position = "none")
  p <- arrangeGrob(p1,p2,ncol = 2)
  ggsave(p,filename  =  paste0(f,".bar.pdf"))
}

bar_ldsc_z <- function(f,plot_extended=FALSE) {
  library(ggplot2)
  library(gridExtra)
  library(grid)
  
  ldsc_results <- f
  GWAS <- read.table(ldsc_results,header  =  T)
  GWAS$Category <- gsub("L2_0","",as.character(GWAS$Category))
  GWAS$Category <- gsub("_"," ",GWAS$Category)
  GWAS$Category <- gsub("\\."," ",GWAS$Category)
  GWAS$Category <- gsub("bed","",GWAS$Category)
  GWAS$Category <- gsub(" $","",GWAS$Category)
  GWAS$Category <- gsub(" 2$","",GWAS$Category)
  
  GWAS$log10p <- -log10(1-pnorm(GWAS$Coefficient_z.score))
  
  GWAS <- GWAS[!GWAS$Category  == "base",]
  
  if(plot_extended==FALSE){
    GWAS <- GWAS[!grepl("500",GWAS$Category),]
  }
  
  GWAS <- GWAS[order(GWAS$log10p),]
  GWAS$Category <- factor(GWAS$Category,levels  =  GWAS$Category)
  
  #5%FDR threshold
  #threshold = max(GWAS$Enrichment_p[p.adjust(GWAS$Enrichment_p,method="fdr")<=0.05])
  
  p2 <- ggplot(GWAS, aes(y = log10p, x = Category,fill=Category)) + coord_flip() 
  p2 <- p2 + geom_bar(position = "dodge", stat = "identity", color = "black") 
  p2 <- p2 + ylab(expression('-log'[10]*'(pvalue)')) 
  p2 <- p2 + geom_hline(yintercept = -log10(0.05/nrow(GWAS)))
  #p2 <- p2 + geom_hline(yintercept = -log10(threshold),linetype=1)
  p2 <- p2 + theme_bw()
  p2 <- p2 + theme(axis.text.y = element_text(angle = 0,hjust = 1,vjust = 0.5,size = 8),legend.position="none")
  #p2 <- p2 + theme(axis.title.y=element_blank(),
  #                 axis.text.y=element_blank(),
  #                 axis.ticks.y=element_blank(),legend.position = "none")
  ggsave(p2,filename  =  paste0(f,".bar.Z.pdf"))
}


#f <- "~/Documents/Data/Projects/Best_tissue_SCZ/Differential_expression/All_cells/"
#f <- "~/Documents/Data/Projects/SCZ2_asian/Report/LDSC_output/SCZ2_all_tag_snps.EUR.merged.ALL.0.3_R2_cluster.1KG_phase3_essentials.median_af.bed.extended_LD.bed_dir.results"
#bar_ldsc ("/Users/julienbryios/Documents/Data/Projects/SCZ2_asian/Report/LDSC_output/SCZ_EAS.ATAC.Neanderthal.results")
#bar_ldsc ("/Users/julienbryios/Documents/Data/Projects/AN_OCD_meta/Data/Partition_h2/AN_OCD_Meta_LDSC_ATAC_Neanderthal.results")
#bar_ldsc("/Users/julienbryios/Documents/Data/Projects/ATAC-seq_acc/SCZ2_ATACacc_int.results")
#bar_ldsc("/Users/julienbryios/Documents/Data/Projects/STARseq/SCZ2_WGSS_int.results")
#bar_ldsc("/Users/julienbryios/Desktop/SCZ2_ATAC_nolimit.bed3_dir2.results")
#bar_ldsc("/Users/julienbryios/Documents/Data/Projects/ATAC-seq/Revision/LD_score/SCZ2_ATAC_300bp.bed2_dir2.results")
#bar_ldsc("/Users/julienbryios/Documents/Data/Projects/AN_OCD_meta/Results/Partition_h2/AN_OCD_Meta_LDSC_53.results")
bar_ldsc("/Users/julienbryios/Documents/Data/Projects/SCZ3/Functional_Genomics/v2_12_2_2018/partitioned_h2/Results/SCZ3_Finucane_no_annot.results")
bar_ldsc_z("/Users/julienbryios/Documents/Data/Projects/SCZ3/Functional_Genomics/v2_12_2_2018/partitioned_h2/Results/SCZ3_Finucane_no_annot.results")
bar_ldsc("/Users/julienbryios/Documents/Data/Projects/SCZ3/Functional_Genomics/v2_12_2_2018/partitioned_h2/Results/SCZ3_Gazal_no_annot.results")
bar_ldsc_z("/Users/julienbryios/Documents/Data/Projects/SCZ3/Functional_Genomics/v2_12_2_2018/partitioned_h2/Results/SCZ3_Gazal_no_annot.results")

files <- list.files("/Users/julienbryios/Documents/Data/Projects/ATAC-seq_acc/LDSC/Differentially_accessible/Results/",".results",full.names = T)
files <- list.files("/Users/julienbryios/Documents/Data/Projects/ATAC-seq_acc/LDSC/Differentially_accessible/Diff_accessible_Tissue_nuclei_v2/Results/",".results",full.names = T)

map(files,bar_ldsc)
map(files,bar_ldsc_z)
