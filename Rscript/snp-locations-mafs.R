########################################
####### get SNPs locations + MAFs ######
########################################
########################################
######## by Kaarina Kowalec 12/2018 ####
########################################

# load editing packages
library(dplyr)
library(data.table)
library(purrr)

# load in rsIDs
rsids <- fread("insertnamehere.txt") %>%
  as_tibble()

#########################################
##### GRCH37 ############################
#########################################

# Locations: dbSNP  
library(SNPlocs.Hsapiens.dbSNP144.GRCh37) # load package
snp_g37 <- SNPlocs.Hsapiens.dbSNP144.GRCh37 # Define formal class object
rsids_g37 <- snpsById(snp_g37, rsids, ifnotfound="drop") # grab locations

# MAF: GnomDB
library(MafDb.gnomAD.r2.1.hs37d5) # load package
mafdb.gnomAD.37 <- MafDb.gnomAD.r2.1.hs37d5 # define maf database
maf.gnomAD.g37 <- as.data.frame(gscores(mafdb.gnomAD.37, rsids_g37)) # define maf based on gnomAD
names(maf.gnomAD.g37)[6]<-"gnomadAF"

# MAF: TOPMed
library(MafDb.TOPMed.freeze5.hg19) # load package
mafdb.topmed.g37 <- MafDb.TOPMed.freeze5.hg19 # define maf database
seqlevelsStyle(cyp_g37) <- seqlevelsStyle(mafdb.topmed.g37) 
maf.topmed.g37 <- as.data.frame(gscores(mafdb.topmed.g37, rsids_g37)) # define maf based on TOPMed
maf.topmed.g37 <- select(maf.topmed.g37,-seqnames)
names(maf.topmed.g37)[5]<-"topmedAF"

# Merge all together
maf.g37 <- inner_join(maf.gnomAD.g37,maf.topmed.g37)

# remove un-necessary dataframes
rm(maf.gnomAD.g37,maf.topmed.g37, mafdb.gnomAD.37,mafdb.topmed.g37,snp_g37, rsids_g37)

#########################################
##### GRCH38 ############################
#########################################

# Locations: dbSNP
library(SNPlocs.Hsapiens.dbSNP151.GRCh38) # load package
snp_g38 <- SNPlocs.Hsapiens.dbSNP151.GRCh38 # Define formal class object
rsids_g38 <- snpsById(snp_g38, rsids, ifnotfound="drop") # grab locationss

# MAF: GnomDB 
library(MafDb.gnomAD.r2.0.1.GRCh38) # load package
mafdb.gnomAD.hg38 <- MafDb.gnomAD.r2.0.1.GRCh38 # define maf database
maf.gnomAD.g38 <- as.data.frame(gscores(mafdb.gnomAD.hg38, rsids_g38)) # define maf based on gnomAD
names(maf.gnomAD.g38)[6]<-"gnomadAF"

# MAF: TOPMed
library(MafDb.TOPMed.freeze5.hg38) # load package
mafdb.topmed.hg38 <- MafDb.TOPMed.freeze5.hg38 # define maf database
seqlevelsStyle(cyp_g38) <- seqlevelsStyle(mafdb.topmed.hg38) 
maf.topmed.g38 <- as.data.frame(gscores(mafdb.topmed.hg38, rsids_g38)) # define maf based on TOPMed
maf.topmed.g38 <- select(maf.topmed.g38,-seqnames)
names(maf.topmed.g38)[5]<-"topmedAF"

# Merge all together
maf.g38 <- inner_join(maf.gnomAD.g38,maf.topmed.g38)

# remove un-necessary dataframes
rm(maf.topmed.g38, maf.gnomAD.g38, mafdb.gnomAD.hg38, mafdb.topmed.hg38, rsids_g38, snp_g38, rsids)