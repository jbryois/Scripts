library(Matrix)
library(lme4)
library(GenABEL)

data_file <- "/data/2backup/common/data/funpopgen_projects/systemotwin/data/expression_OURF/systemoTwin_232_v15_gene.coverage_scaled_10per.new_quantification.paper.txt"

#data_file <- "/data/2backup/common/data/funpopgen_projects/systemotwin/data/expression_OURF/systemoTwin_232_v15_gene.coverage_scaled_10per.new_quantification.paper.txt.GC_mean_lane_LibraryPrepDate_INSERT_SIZE_MODE.txt"

pedigree <- "/data/jbryois/projects/systemotwin/Solar_modeling/Heritability_from_scratch/pedigree_t1_no_outliers.csv"
pedigree_t2 <- "/data/jbryois/projects/systemotwin/Solar_modeling/Heritability_from_scratch/pedigree_t2_no_outliers.csv"

    reads = read.table(data_file, header=T)
    tmp <- names(reads)[-c(1:4)]
   	
   	rownames(reads) <- reads$TargetID
   	reads <- reads[-c(1:4)]
	
    covs = read.csv(pedigree,header=T,stringsAsFactors=F)
    covs$id_parsed <- ifelse(substr(covs$id,1,1)=="B", substr(covs$id,2,nchar(covs$id)-1),ifelse(substr(covs$id,1,1)=="L",substr(covs$id,4,nchar(covs$id)-1),substr(covs$id,5,nchar(covs$id)-1)))
    covs$hhid <- paste("HH-",substr(covs$fa,3,nchar(covs$fa)),sep="")
	covs$mztwin <- ifelse(substr(covs$mztwin,1,1)=="M",paste("MZ-",substr(covs$fa,3,nchar(covs$fa)),sep=""),paste("DZ-",covs$id_parsed,sep=""))

	
    covs_t2 = read.csv(pedigree_t2,header=T,stringsAsFactors=F)
    covs_t2$id_parsed <- ifelse(substr(covs_t2$id,1,1)=="B", substr(covs_t2$id,2,nchar(covs_t2$id)-1),ifelse(substr(covs_t2$id,1,1)=="L",substr(covs_t2$id,4,nchar(covs_t2$id)-1),substr(covs_t2$id,5,nchar(covs_t2$id)-1)))
    covs_t2$hhid <- paste("HH-",substr(covs_t2$fa,3,nchar(covs_t2$fa)),sep="")
	covs_t2$mztwin <- ifelse(substr(covs_t2$mztwin,1,1)=="M",paste("MZ-",substr(covs_t2$fa,3,nchar(covs_t2$fa)),sep=""),paste("DZ-",covs_t2$id_parsed,sep=""))

	covs <- rbind(covs,covs_t2)
	
	reads <- reads[as.character(covs$id)]
	normal_transformed_expression <- apply(reads,1, function(x) rntransform(as.numeric(x)))
	normal_transformed_expression <- t(normal_transformed_expression)
	#rntransform(as.numeric(reads[2312,]))==normal_transformed_expression[2312,]
	reads <- normal_transformed_expression
	
    n_peaks <- nrow(reads)
    
    DEL <- as.factor(as.matrix(covs['mztwin'])) # Same number for MZ twins (family ID), different number for DZ twins (twin ID + 00000)
    DZ <- as.factor(as.matrix(covs['hhid'])) # This is coded with the family ID.
	TP_random <- as.factor(as.matrix(covs['id_parsed'])) # This is coded with the family ID.

	info <- read.table("/data/2backup/common/data/funpopgen_projects/systemotwin/data/systemoTwin_TP_info_final.txt",header=T)
	info_ordered <- info[match(as.character(covs$id),info$sampleID),]
	#info_ordered$sampleID==as.character(covs$id)
	TP_continuous <- ifelse(info_ordered$Time_point==1,0,info_ordered$Difference_in_visits_.years.)
	Age <- info_ordered$Age_first_visit
	
	covariates <- read.table("/data/2backup/common/data/funpopgen_projects/systemotwin/data/systemoTwin_masterQC_final.txt",header=T,stringsAsFactors=F)
	covariates$sample[covariates$sample=="LE_30262B"] <- "LE_32062B"
	covariates$sample[covariates$sample=="LE_73681B"] <- "temp"
	covariates$sample[covariates$sample=="B73681A"] <- "LE_73681B"
	covariates$sample[covariates$sample=="temp"] <- "B73681A"
	covariates$sample[covariates$sample=="LE_73682B"] <- "temp"
	covariates$sample[covariates$sample=="B73682A"] <- "LE_73682B"
	covariates$sample[covariates$sample=="temp"] <- "B73682A"

	covariates_ordered <- covariates[na.omit(match(info_ordered$sampleID,covariates$sample)),]
	covariates_ordered$sample==info_ordered$sampleID
	covariates_ordered$sample==as.character(covs$id)
	#write.table(covariates_ordered,"/data/2backup/common/data/funpopgen_projects/systemotwin/data/systemoTwin_masterQC_final_sample_mislabel_corrected.txt",quote=F,sep="\t",col.names=T,row.names=F)
	
	
	GC <- covariates_ordered$GC_mean 
	lane <- covariates_ordered$lane 
	libprep <- covariates_ordered$LibraryPrepDate 
	insertsize <- covariates_ordered$INSERT_SIZE_MODE 

	gene_name <- character()
	h2_vector <- numeric()
	c2_vector <- numeric()
	beta_tp_vector <- numeric()
	pvalue_tp_vector <- numeric()
	for(j in 1:n_peaks){
		gene_name <- c(gene_name,rownames(reads)[j])
		lmer_totv <- lmer(as.matrix(as.double(reads[j,])) ~ TP_continuous + Age + GC + (1 | DZ) + (1 | DEL) + (1 | TP_random) + (1 | lane) + (1 | libprep) + (1 | insertsize))
		lmer_limited <- lmer(as.matrix(as.double(reads[j,])) ~ 1 + Age + GC + (1 | DZ) + (1 | DEL) + (1 | TP_random) + (1 | lane) + (1 | libprep) + (1 | insertsize))
		beta_tp <- fixef(lmer_totv)[2]
		pvalue_tp <- anova(lmer_totv,lmer_limited)[8][2,]
		variance_E <- attr(VarCorr(lmer_totv),"sc")^2
		variance_DEL <- attr(VarCorr(lmer_totv)$DEL,"stddev")^2
		variance_DZ <- attr(VarCorr(lmer_totv)$DZ,"stddev")^2
		h2_fly <- 2*variance_DEL/(variance_DEL+variance_E+variance_DZ)
		c2_fly <- (variance_DZ-variance_DEL)/(variance_DEL+variance_E+variance_DZ)
		h2_vector <- c(h2_vector,h2_fly)
		c2_vector <- c(c2_vector,c2_fly)
		beta_tp_vector <- c(beta_tp_vector,beta_tp)
		pvalue_tp_vector <- c(pvalue_tp_vector, pvalue_tp)
		cat(j,"\t",rownames(reads)[j],"\t",h2_fly,"\t",c2_fly,"\t",beta_tp,"\t",pvalue_tp,"\n")	
	}

heritability <- as.data.frame(cbind(gene_name,as.numeric(h2_vector),as.numeric(c2_vector),as.numeric(beta_tp_vector),as.numeric(pvalue_tp_vector)))
colnames(heritability) <- c("name","h2","c2","beta_tp","pvalue_tp")
write.table(heritability,"/data/jbryois/projects/systemotwin/Paper_analysis/Differentially_expressed_mixed_model/DE_genes_h2.txt", sep="\t",quote=F, row.names=F,col.names=T)

#################ANALYSE RESULTS

library(qvalue)
heritability <- read.table("/data/jbryois/projects/systemotwin/Paper_analysis/Differentially_expressed_mixed_model/DE_genes_h2.txt", header=T)
qvalue_h2 <- qvalue(as.numeric(as.matrix(heritability$pvalue_tp)),pi0.method="bootstrap")
1-qvalue_h2$pi0
#0.4616966

threshold <- max(qvalue_h2$pvalues[qvalue_h2$qvalues <= 0.05]) 
DE_5_FDR <- heritability[as.numeric(as.matrix(heritability$pvalue_tp))< threshold,]
DE_5_FDR_ordered <- DE_5_FDR[order(as.numeric(as.matrix(DE_5_FDR$pvalue_tp))),]

DE_5_FDR_up <- DE_5_FDR_ordered[as.numeric(as.matrix(DE_5_FDR_ordered$beta_tp))>0,]
DE_5_FDR_down <- DE_5_FDR_ordered[as.numeric(as.matrix(DE_5_FDR_ordered$beta_tp))<0,]

write.table(DE_5_FDR,"/data/jbryois/projects/systemotwin/Paper_analysis/Differentially_expressed_mixed_model/DE_5FDR_combined_t1_t2.GC_mean_lane_LibraryPrepDate_INSERT_SIZE_MODE.h2r",col.names=T,row.names=F,sep="\t",quote=F)
write.table(DE_5_FDR_up,"/data/jbryois/projects/systemotwin/Paper_analysis/Differentially_expressed_mixed_model/DE_5FDR_combined_t1_t2.GC_mean_lane_LibraryPrepDate_INSERT_SIZE_MODE.h2r_up",col.names=T,row.names=F,sep="\t",quote=F)
write.table(DE_5_FDR_down,"/data/jbryois/projects/systemotwin/Paper_analysis/Differentially_expressed_mixed_model/DE_5FDR_combined_t1_t2.GC_mean_lane_LibraryPrepDate_INSERT_SIZE_MODE.h2r_down",col.names=T,row.names=F,sep="\t",quote=F)

threshold <- max(qvalue_h2$pvalues[qvalue_h2$qvalues <= 0.5])

DE_80_FDR_more <- heritability[as.numeric(as.matrix(heritability$pvalue_tp))> threshold,]
write.table(DE_80_FDR_more,"/data/jbryois/projects/systemotwin/Paper_analysis/Differentially_expressed_mixed_model/DE_more_50FDR_combined_t1_t2.GC_mean_lane_LibraryPrepDate_INSERT_SIZE_MODE.h2r_down",col.names=T,row.names=F,sep="\t",quote=F)

threshold <- max(qvalue_h2$pvalues[qvalue_h2$qvalues <= 0.01]) 
DE_5_FDR <- heritability[as.numeric(as.matrix(heritability$pvalue_tp))< threshold,]
DE_5_FDR_ordered <- DE_5_FDR[order(as.numeric(as.matrix(DE_5_FDR$pvalue_tp))),]

DE_5_FDR_up <- DE_5_FDR_ordered[as.numeric(as.matrix(DE_5_FDR_ordered$beta_tp))>0,]
DE_5_FDR_down <- DE_5_FDR_ordered[as.numeric(as.matrix(DE_5_FDR_ordered$beta_tp))<0,]

write.table(DE_5_FDR,"/data/jbryois/projects/systemotwin/Paper_analysis/Differentially_expressed_mixed_model/DE_1FDR_combined_t1_t2.GC_mean_lane_LibraryPrepDate_INSERT_SIZE_MODE.h2r",col.names=T,row.names=F,sep="\t",quote=F)
write.table(DE_5_FDR_up,"/data/jbryois/projects/systemotwin/Paper_analysis/Differentially_expressed_mixed_model/DE_1FDR_combined_t1_t2.GC_mean_lane_LibraryPrepDate_INSERT_SIZE_MODE.h2r_up",col.names=T,row.names=F,sep="\t",quote=F)
write.table(DE_5_FDR_down,"/data/jbryois/projects/systemotwin/Paper_analysis/Differentially_expressed_mixed_model/DE_1FDR_combined_t1_t2.GC_mean_lane_LibraryPrepDate_INSERT_SIZE_MODE.h2r_down",col.names=T,row.names=F,sep="\t",quote=F)

threshold <- max(qvalue_h2$pvalues[qvalue_h2$qvalues <= 0.4])

DE_80_FDR_more <- heritability[as.numeric(as.matrix(heritability$pvalue_tp))> threshold,]
write.table(DE_80_FDR_more,"/data/jbryois/projects/systemotwin/Paper_analysis/Differentially_expressed_mixed_model/DE_more_40FDR_combined_t1_t2.GC_mean_lane_LibraryPrepDate_INSERT_SIZE_MODE.h2r_down",col.names=T,row.names=F,sep="\t",quote=F)
