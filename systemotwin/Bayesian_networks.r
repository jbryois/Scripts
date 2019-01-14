library(GenABEL)
library(bnlearn)

#read in the arguments
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2) {
    print("Error, this script requires 2 arguments:")
    print("Number_of_covariate")
    stop("Exiting script")
}

start <- arguments[1]
end <- arguments[2]

expression <- read.table("/data/2backup/common/data/funpopgen_projects/systemotwin/data/expression_OURF/systemoTwin_232_v15_gene.coverage_scaled_10per.new_quantification.paper.txt.GC_mean_lane_LibraryPrepDate_INSERT_SIZE_MODE.TP2corrected_for_delta_tp.paper.txt",header=T)
expression_names <- colnames(expression)
normal_transformed_expression <- apply(expression[-c(1:4)],1, function(x) rntransform(as.numeric(x)))
normal_transformed_expression <- t(normal_transformed_expression)
#rntransform(as.numeric(expression[-c(1:4)][2312,]))==normal_transformed_expression[2312,]
expression <- as.data.frame(cbind(expression[c(1:4)],normal_transformed_expression))
colnames(expression) <- expression_names

DE_genes_all_DE_up <- read.table("/data/jbryois/projects/systemotwin/Paper_analysis/Differentially_expressed_mixed_model/DE_1FDR_combined_t1_t2.GC_mean_lane_LibraryPrepDate_INSERT_SIZE_MODE.h2r_up",header=T)
DE_genes_all_DE_down <- read.table("/data/jbryois/projects/systemotwin/Paper_analysis/Differentially_expressed_mixed_model/DE_1FDR_combined_t1_t2.GC_mean_lane_LibraryPrepDate_INSERT_SIZE_MODE.h2r_down",header=T)

DE_genes_all_DE_up_expr <- expression[expression$TargetID%in%as.character(DE_genes_all_DE_up$name),]
DE_genes_all_DE_down_expr <- expression[expression$TargetID%in%as.character(DE_genes_all_DE_down$name),]

info <- read.table("/data/2backup/common/data/funpopgen_projects/systemotwin/data/systemoTwin_TP_info_final.txt",header=T)
info_ordered <- info[match(colnames(expression)[-c(1:4)],info$sampleID),]
info_ordered$sampleID==colnames(expression)[-c(1:4)]
TP_continuous <- ifelse(info_ordered$Time_point==1,0,info_ordered$Difference_in_visits_.years.)
LL <- TP_continuous

######COMPARE MODELS FUNCTION#######
## LL=Time point GG= gene upregulated TT= gene downregulated
#3 models: linearGE, linearM, indep
compareModels <- function(dat){

  AIC <- numeric()
  BIC <- numeric()
  LOG <- numeric()

  linearGE <- empty.graph(names(dat))
  modelstring(linearGE) <- "[LL][GG|LL][TT|GG]"
  BIC <- c(BIC,score(linearGE, dat, type = "bic-g"))
  AIC <- c(AIC,score(linearGE, dat, type = "aic-g"))
  LOG <- c(LOG,score(linearGE, dat, type = "loglik-g"))
  
  linearM <- empty.graph(names(dat))
  modelstring(linearM) <- "[LL][TT|LL][GG|TT]"
  BIC <- c(BIC,score(linearM, dat, type = "bic-g"))
  AIC <- c(AIC,score(linearM, dat, type = "aic-g"))
  LOG <- c(LOG,score(linearM, dat, type = "loglik-g"))
  
  indep <- empty.graph(names(dat))
  modelstring(indep) <- "[LL][TT|LL][GG|LL]"
  BIC <- c(BIC,score(indep, dat, type = "bic-g"))
  AIC <- c(AIC,score(indep, dat, type = "aic-g"))
  LOG <- c(LOG,score(indep, dat, type = "loglik-g"))
  
  #indepGE <- empty.graph(names(dat))
  #modelstring(indepGE) <- "[LL][GG][TT|LL:GG]"
  #BIC <- c(BIC,score(indepGE, dat, type = "bic-g"))
  #AIC <- c(AIC,score(indepGE, dat, type = "aic-g"))
  #LOG <- c(LOG,score(indepGE, dat, type = "loglik-g"))

  #indepM <- empty.graph(names(dat))
  #modelstring(indepM) <- "[LL][TT][GG|LL:TT]"
  #BIC <- c(BIC,score(indepM, dat, type = "bic-g"))
  #AIC <- c(AIC,score(indepM, dat, type = "aic-g"))
  #LOG <- c(LOG,score(indepM, dat, type = "loglik-g"))

  res <- as.data.frame(t(data.frame(bic=BIC, aic=AIC, loglike=LOG)))
  names(res) <- c("linearGE","linearM","indep")
  win <- apply(res,1,function(x) names(res)[which(x==max(x))])
  tmp <- res[row.names(res)=="aic", colnames(res)!=win[2]]
  relativeLikelihood <- exp(-(res[row.names(res)=="aic",win[2]]-tmp)/2)
  names(relativeLikelihood) <- names(tmp)
  relativeLikelihood[[win[2]]] <- 1
  rL <- as.numeric(relativeLikelihood) 
  names(rL) <- names(relativeLikelihood)
  rL <- sort(rL, decreasing=TRUE)
  #print(win)
  rL
}

cm <- NULL
for(i in start:end){
	gene_1_name <- as.character(DE_genes_all_DE_up_expr$TargetID[i])
	GG <- as.numeric(DE_genes_all_DE_up_expr[-c(1:4)][DE_genes_all_DE_up_expr$TargetID==gene_1_name,])
	for(j in 1:nrow(DE_genes_all_DE_down_expr)){
		gene_2_name <- as.character(DE_genes_all_DE_down_expr$TargetID[j])
		TT <- as.numeric(DE_genes_all_DE_down_expr[-c(1:4)][DE_genes_all_DE_down_expr$TargetID==gene_2_name,])
		dat <- as.data.frame(cbind(LL, GG, TT))	
		compMod <- compareModels(dat)
		cm <- rbind(cm,c(gene_1_name,gene_2_name, cor(LL,TT),cor(LL,GG), cor(TT,GG), compMod["linearGE"], compMod["linearM"], compMod["indep"]))
	}}	
cm <- as.data.frame(cm)
#colnames(cm) <- c("gene_up","gene_down", "cor_time_down_gene","cor_time_up_gene", "cor_up_down_gene", "time_UP_DOWN", "time_DOWN_UP","INDEP")

write.table(cm, paste("/data/jbryois/projects/systemotwin/Paper_analysis/Bayesian_network_for_differentially_expressed_genes/all_models_DE_1FDR.txt_",start,"_",end,sep=""), col.names=F,row.names=F,quote=F,sep="\t")