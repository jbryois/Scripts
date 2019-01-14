library(edgeR)
library(tidyverse)

#Load exon file
d <- read_tsv("~/Downloads/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_exon_reads.txt")
d <- as.data.frame(d)

#Set rownames as the exon names
rownames(d) <- d$exon_id

#Create a DGEList object for edgeR
d2 <- DGEList(counts=d[-1])

#Get counts per millions for filtering. For each sample: read counts*10^6/(sum Reads counts)
#These are not normalised for exon size. Not good to compare expression between different exons. 
cpm_only <- cpm(d2, normalized.lib.sizes=FALSE)

#Optional: filter out exon with low expression accross tissues/samples
#Keep only exon with at least 10 cpm  in at least 50 samples
keep <- rowSums(cpm_only>10)>= 50
d2 <- d2[keep,]

#Get the TMM normalisation factors
d2 <- calcNormFactors(d2)

#Get final normalised values
cpm_TMM <- cpm(d2, normalized.lib.sizes=TRUE) 

#Write normalized read counts
write.table(cpm_TMM,"~/Downloads/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_exon_reads.filtered.TMM.txt",quote=F,sep="\t",col.names=T,row.names=T)
