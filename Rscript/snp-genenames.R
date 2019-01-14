########################################
####### get SNPs gene names  ###########
########################################
########################################
######## by Kaarina Kowalec 1/2019 ####
########################################

# rsIDs to Gene Names
# BiocManager::install("biomaRt", version = "3.8")

library(biomaRt)
#Mart used to map SNPs to Ensembl Gene IDs
grch37.snp = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", 
                     path="/biomart/martservice",dataset="hsapiens_snp")
#Mart used to map Ensembl Gene IDs to Gene name
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", 
                 path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

table1 <- getBM(attributes = c("refsnp_id", "ensembl_gene_stable_id"), 
                filters = "snp_filter", 
                values = yourSNPlist, 
                mart = grch37.snp)
table2 <- getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                               "external_gene_source"),
                filters = "ensembl_gene_id", 
                values =  table1$ensembl_gene_stable_id, 
                mart = grch37)

results <- merge(table1,table2, by.x = "ensembl_gene_stable_id", 
                 by.y="ensembl_gene_id", all.x=T)