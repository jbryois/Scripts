################################################################################
# Script to extract variances values for Heritability analysis from Eurobats
#   expression data Freezev2 January 2013.
# Modify from Kerrin Small.
# Author: Ana ViÒuela (ana.vinuela@kcl.ac.uk)
# Date: 13/01/2013
################################################################################
#read in the arguments
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1) {
    print("Error, this script requires 1 arguments:")
    print("START")
    stop("Exiting script")
}
library(Matrix)
library(lme4)
library(GenABEL)

start <- arguments[1]

data_file <- "/data/jbryois/projects/systemotwin/Solar_modeling/Heritability/prepare_exp_gen_files/t1.GC_mean.corrected.exp.csv"
pedigree <- "/data/jbryois/projects/systemotwin/Solar_modeling/Heritability/prepare_exp_gen_files/pedigree_file_t1.csv"


    reads = read.csv(data_file, header=T)
   	reads <- t(reads)
   	colnames(reads) <- reads[1,]
   	reads <- reads[-1,]
   	reads <- as.data.frame(reads)
    covs = read.csv(pedigree,header=T)
    covs$hhid <- paste("HH-",substr(covs$id,1,nchar(covs$id)-1),sep="")
	covs$mztwin <- ifelse(substr(covs$mztwin,1,1)=="M",paste("MZ-",substr(covs$id,1,nchar(covs$id)-1),sep=""),paste("DZ-",substr(covs$id,1,nchar(covs$id)),sep=""))
	unique_family <- names(which(table(covs$hhid)==2))
	unique_family_of_single_individuals <- names(which(table(covs$hhid)==1))
	bootstrap_unique_family <- sample(unique_family,length(unique_family),replace=T)
	family_to_keep_for_bootstrap <- as.data.frame(c(unique_family_of_single_individuals,bootstrap_unique_family))
	colnames(family_to_keep_for_bootstrap) <- "hhid"
	
	####BOOTSTRAPED COVS and READS
	covs <- merge(covs,family_to_keep_for_bootstrap,by="hhid")
	reads <- reads[as.character(covs$id)]
	normal_transformed_expression <- apply(reads,1, function(x) rntransform(as.numeric(x)))
	normal_transformed_expression <- t(normal_transformed_expression)
	#rntransform(as.numeric(reads[2312,]))==normal_transformed_expression[2312,]
	reads <- normal_transformed_expression
    
    #### CHANGE NAMES OF THE COVS THAT ARE DUPLICATED
	covs_temp <- covs
	id_list <- NULL
	for(i in 1:nrow(covs)){
		id_list <- c(id_list,covs$id[i])	
		number_of_occurence_in_list <- table(id_list)[as.character(covs$id[i])]
		covs_temp$hhid[i] <- paste(covs$hhid[i],number_of_occurence_in_list,sep="_")
		covs_temp$mztwin[i] <- paste(covs$mztwin[i],number_of_occurence_in_list,sep="_")
	}
	covs <- covs_temp

    DEL <- as.factor(as.matrix(covs['mztwin'])) # Same number for MZ twins (family ID), different number for DZ twins (twin ID + 00000)
    DZ <- as.factor(as.matrix(covs['hhid'])) # This is coded with the family ID.
    n_peaks <- nrow(reads)

	gene_name <- character()
	h2_vector <- numeric()
	c2_vector <- numeric()
	for(j in 1:n_peaks){
		gene_name <- c(gene_name,rownames(reads[j,]))
		lmer_totv <- lmer(as.matrix(as.double(reads[j,])) ~ 1 + (1 | DZ) + (1 | DEL))
		variance_E <- attr(VarCorr(lmer_totv),"sc")^2
		variance_DEL <- attr(VarCorr(lmer_totv)$DEL,"stddev")^2
		variance_DZ <- attr(VarCorr(lmer_totv)$DZ,"stddev")^2
		h2_fly <- 2*variance_DEL/(variance_DEL+variance_E+variance_DZ)
		c2_fly <- (variance_DZ-variance_DEL)/(variance_DEL+variance_E+variance_DZ)
		h2_vector <- c(h2_vector,h2_fly)
		c2_vector <- c(c2_vector,c2_fly)
		cat(j,"\t",rownames(reads[j,]),"\t",h2_fly,"\t",c2_fly,"\n")	
	}

heritability <- as.data.frame(cbind(rownames(reads),as.numeric(h2_vector),as.numeric(c2_vector)))
colnames(heritability) <- c("name","h2","c2")
write.table(heritability,paste("/data/jbryois/projects/systemotwin/Solar_modeling/Heritability_Ana/Bootstrap_results/t1.GC_mean.corrected.h2r_",start,sep=""), sep="\t",quote=F, row.names=F,col.names=T)
















