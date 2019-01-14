load("/data/jbryois/projects/Alspac/Interaction/GWAS-cis/cis_eQTLs_genotype.Robject")
genotype_cis<- genotype_all
genotype_cis <- unique(genotype_cis)

#expression_all <- read.table("/data/jbryois/projects/Alspac/genotype_expression_without_outliers/869expression_without_outliers.exp", header=T)
expression_all <- read.table("/data/jbryois/projects/Alspac/genotype_expression_without_outliers/869expression_without_outliers_normal_transformed.exp", header=T)
trans_effects_of_cis_eqtls <- read.table("/data/jbryois/projects/Alspac/Revisions/trans_analysis_of_cis_eQTLs/Trans_effect_of_cis_eqtls_5FDR.txt_filtered_col_5_final_table_arranged_for_paper_with_rep_MUTHER_GEUVADIS_for_paper_1MAF_1KG.txt.txt", header=T)

######COMPARE MODELS FUNCTION#######
## LL=genotype GG=trans exp TT=cis exp
#3 models: linearGE, linearM, indep
compareModels <- function(dat){
  library(bnlearn)

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

############## CIT
CausalityTestJM <- function(LL,GG,TT){
 	no.bootstrap = 50
 	# remove missing values
  	sel = (!is.na(LL)) & (!is.na(GG)) & (!is.na(TT))
  	
   dat_f = as.data.frame(cbind(LL,GG,TT),stringsAsFactors=FALSE)
   dat_f = dat_f[sel,]
   names(dat_f) = c("L","G","T")
   Lf = as.factor(dat_f$L)
   dat_f$L = as.integer(Lf) - 1
   llevels = as.integer(levels(as.factor(dat_f$L)))
   dfL = length(llevels) - 1

   #### Causal p-value
                
   pvec = rep(NA,4)
	if(dfL > 2){
	   return(c("NA","NA","NA"))

	}
   if(dfL == 2){
	dat_f$L1 = ifelse(dat_f$L == 1,1,0)
	dat_f$L2 = ifelse(dat_f$L == 2,1,0)

	fit0 = lm(T ~ 1,data=dat_f)
	fit1 = lm(T ~ L1 + L2,data=dat_f)
	fit2 = lm(G ~ T,data=dat_f)
	fit3 = lm(T ~ G,data=dat_f)
	fit4 = lm(G ~ T + L1 + L2,data=dat_f)
	fit5 = lm(T ~ G + L1 + L2,data=dat_f)

	pvec[1] = anova(fit0,fit1)$"Pr(>F)"[2]
	pvec[2] = anova(fit2,fit4)$"Pr(>F)"[2]
	pvec[3] = anova(fit1,fit5)$"Pr(>F)"[2]

	f_ = anova(fit3,fit5)$F[2]

	fit1G = lm(G ~ L1 + L2,data=dat_f)

	alg = summary(fit1G)$coefficients["(Intercept)",1]
	blg1 = summary(fit1G)$coefficients["L1",1]
	blg2 = summary(fit1G)$coefficients["L2",1]

	alt = summary(fit1)$coefficients["(Intercept)",1]
	blt1 = summary(fit1)$coefficients["L1",1]
	blt2 = summary(fit1)$coefficients["L2",1]

	dat_f$eG = resid(fit1G)
	dat_f$eT = resid(fit1)       

	ss = dim(dat_f)[1]
	fvecr = rep(NA,no.bootstrap)
	fvecr_r = rep(NA,no.bootstrap)

	for(i in 1:no.bootstrap){
		nni <- trunc(1 + ss*runif(ss, 0, 1)) ;
		dat_f$G_ = alg + blg1*dat_f$L1 + blg2*dat_f$L2 + dat_f$eG[nni]

		fit_0 = lm(T ~ G_,data=dat_f)
		fit_1 = lm(T ~ G_ + L1 + L2,data=dat_f)
		fvecr[i] = anova(fit_0,fit_1)$F[2]

		dat_f$T_ = alt + blt1*dat_f$L1 + blt2*dat_f$L2 + dat_f$eT[nni]

		fit_0 = lm(G ~ T_,data=dat_f)
		fit_1 = lm(G ~ T_ + L1 + L2,data=dat_f)
		fvecr_r[i] = anova(fit_0,fit_1)$F[2]
	}
   }#End dfL == 2
   if(dfL == 1){

	dat_f$L1 = ifelse(dat_f$L == 1,1,0)

	fit0 = lm(T ~ 1,data=dat_f)
	fit1 = lm(T ~ L1,data=dat_f)
	fit2 = lm(G ~ T,data=dat_f)
	fit3 = lm(T ~ G,data=dat_f)
	fit4 = lm(G ~ T + L1,data=dat_f)
	fit5 = lm(T ~ G + L1,data=dat_f)

	pvec[1] = anova(fit0,fit1)$"Pr(>F)"[2]
	pvec[2] = anova(fit2,fit4)$"Pr(>F)"[2]
	pvec[3] = anova(fit1,fit5)$"Pr(>F)"[2]

	f_ = anova(fit3,fit5)$F[2]

	fit1G = lm(G ~ L1,data=dat_f)

	alt = summary(fit1)$coefficients["(Intercept)",1]
	blt1 = summary(fit1)$coefficients["L1",1]

	alg = summary(fit1G)$coefficients["(Intercept)",1]
	blg1 = summary(fit1G)$coefficients["L1",1]

	dat_f$eG = resid(fit1G)
	dat_f$eT = resid(fit1)       

	ss = dim(dat_f)[1]
	fvecr = rep(NA,no.bootstrap)
	fvecr_r = rep(NA,no.bootstrap)

	for(i in 1:no.bootstrap){
		nni <- trunc(1 + ss*runif(ss, 0, 1)) ;
		dat_f$G_ = alg + blg1*dat_f$L1 + dat_f$eG[nni]

		fit_0 = lm(T ~ G_,data=dat_f)
		fit_1 = lm(T ~ G_ + L1,data=dat_f)
		fvecr[i] = anova(fit_0,fit_1)$F[2]

		dat_f$T_ = alt + blt1*dat_f$L1 + dat_f$eT[nni]

		fit_0 = lm(G ~ T_,data=dat_f)
		fit_1 = lm(G ~ T_ + L1,data=dat_f)
		fvecr_r[i] = anova(fit_0,fit_1)$F[2]
	}
   } #End dfL == 1

   #####F Method
   fvecr = fvecr[!is.na(fvecr)]
   df1 = anova(fit3,fit5)$Df[2]
   df2 = anova(fit3,fit5)$Res.Df[2]
   fncp = mean(fvecr,na.rm=TRUE)*(df1/df2)*(df2-df1)-df1
   if(fncp < 0) fncp = 0

   ######### Transform F to normal
   npvals = pf(fvecr,df1,df2,ncp=fncp,lower.tail=TRUE)
   nfvecr = qnorm(npvals)

   npf = pf(f_,df1,df2,ncp=fncp,lower.tail=TRUE) #Transform observed F
   zf = qnorm(npf)
   pvec[4] = pnorm(zf,mean=0,sd=sd(nfvecr))

   pvalc = max(pvec)  ###Causal p-value

   #### Reactive p-value
   fit0G = lm(G ~ 1,data=dat_f)
   pvec1 = rep(NA,4)

	pvec1[1] = anova(fit0G,fit1G)$"Pr(>F)"[2]
	pvec1[2] = anova(fit3,fit5)$"Pr(>F)"[2]
	pvec1[3] = anova(fit1G,fit4)$"Pr(>F)"[2]
	f_ = anova(fit2,fit4)$F[2]

   #####F Method
   fvecr_r = fvecr_r[!is.na(fvecr_r)]
   df1 = anova(fit3,fit5)$Df[2]
   df2 = anova(fit3,fit5)$Res.Df[2]
   fncp = mean(fvecr_r,na.rm=TRUE)*(df1/df2)*(df2-df1)-df1
   if(fncp < 0) fncp = 0

   ######### Transform F to normal
   npvals = pf(fvecr_r,df1,df2,ncp=fncp,lower.tail=TRUE)
   nfvecr = qnorm(npvals)

   npf = pf(f_,df1,df2,ncp=fncp,lower.tail=TRUE) #Transform observed F
   zf = qnorm(npf)
   pvec1[4] = pnorm(zf,mean=0,sd=sd(nfvecr))


   pvalr = max(pvec1)  ###Reactive p-value

   ccall = NA
   ccall = ifelse((pvalc < 0.05/93) & (pvalr > 0.05/93),1,ccall)
   ccall = ifelse((pvalc > 0.05/93) & (pvalr < 0.05/93),2,ccall)
   ccall = ifelse((pvalc > 0.05/93) & (pvalr > 0.05/93),3,ccall)
   ccall = ifelse((pvalc < 0.05/93) & (pvalr < 0.05/93),0,ccall)

   return(c(pvalc,pvalr,ccall))
}




#################################
#################################
#################################

rownames(expression_all) <- expression_all$TargetID
rownames(genotype_cis) <- genotype_cis$clone

expression_all <- expression_all[,5:ncol(expression_all)]
genotype_cis <- genotype_cis[,5:ncol(genotype_cis)]

indivs <- colnames(expression_all)
indivs <- indivs[which(indivs %in% colnames(genotype_cis))]
length(indivs)

expression_all <- expression_all[, indivs]
genotype_cis <- genotype_cis[, indivs]

#test
trans_effects_of_cis_eqtls <- trans_effects_of_cis_eqtls[order(trans_effects_of_cis_eqtls$SNP, decreasing=T),]

Nlines <- nrow(trans_effects_of_cis_eqtls)
cm <- NULL
cm2 <- NULL
for (i in 1:Nlines) {
	# elements
	probe_cis <- trans_effects_of_cis_eqtls[i,4]
	probe_trans <- trans_effects_of_cis_eqtls[i,11]
	SNP_genotype <- trans_effects_of_cis_eqtls[i,1]

	# look for cis gene expression
	cis_exp_tmp <- expression_all[which(rownames(expression_all) %in% probe_cis),]
	# look for snp line
	trans_exp_tmp <- expression_all[which(rownames(expression_all) %in% probe_trans),]
	# look for meth line
	snp_genotype_tmp <- genotype_cis[which(rownames(genotype_cis) %in% SNP_genotype),]

	# linear model with interaction
	TT <- as.numeric(cis_exp_tmp)
	
	GG <- as.numeric(trans_exp_tmp)

	#LL <- as.numeric(snp_genotype_tmp)
	LL <- as.numeric(round(snp_genotype_tmp))

	NApos <- c(which(is.na(TT)), which(is.na(GG)), which(is.na(LL)))
	
	if(length(NApos) > 0){
		GG <- GG[-NApos]
		TT <- TT[-NApos]
		LL <- LL[-NApos]
		print("NA datapoints removed!")
	}
			
	dat <- as.data.frame(cbind(LL, GG, TT))	
	compMod <- compareModels(dat)
	compMod2 <- CausalityTestJM(LL,GG,TT)
	print(compMod2)
	compMod2 <- as.data.frame(matrix(compMod2,ncol=3))
	colnames(compMod2) <- c("STC_ICT_pvalue", "SCT_ICT_pvalue", "model")

	cm <- rbind(cm,as.data.frame(c(trans_effects_of_cis_eqtls[i,], compMod[c("linearGE", "linearM", "indep")])))
	cm2 <- rbind(cm2,as.data.frame(cbind(trans_effects_of_cis_eqtls[i,], compMod2)))
	cat(i)
	cat("\n")
}
	
colnames(cm)[c(18,19,20)] <- c("STC", "SCT", "INDEP")


write.table(cm, file = "/data/jbryois/projects/Alspac/Revisions/causal_model/trans_effect_of_cis_eQTLs_with_models.txt_normal_exp_filtered_probe_with_snp" , quote = FALSE, row.names = FALSE, sep = "\t")
write.table(cm2, file = "/data/jbryois/projects/Alspac/Revisions/causal_model/trans_effect_of_cis_eQTLs_with_models_ICT.txt_normal_exp_filtered_probe_with_snp" , quote = FALSE, row.names = FALSE, sep = "\t")

cm <- read.table("/data/jbryois/projects/Alspac/Revisions/causal_model/trans_effect_of_cis_eQTLs_with_models.txt_normal_exp", header=T)
cm2 <- read.table("/data/jbryois/projects/Alspac/Revisions/causal_model/trans_effect_of_cis_eQTLs_with_models_ICT.txt_normal_exp", header=T)

INDEP <- length(which(cm$INDEP == 1))
STC <- length(which(cm$STC == 1))
SCT <- length(which(cm$SCT == 1))

sum <- c(INDEP, STC, SCT)
names(sum) <- c("INDEP", "STC", "SCT")
print(sum)
percentage_sum_best_model <- sum*100/156
 
sum_indep <- sum(cm$INDEP) 
sum_SCT <- sum(cm$SCT) 
sum_STC <- sum(cm$STC)
percentage <- c(sum_STC*100/(sum_indep+sum_SCT+sum_STC), sum_SCT*100/(sum_indep+sum_SCT+sum_STC),sum_indep*100/(sum_indep+sum_SCT+sum_STC))
 

STC_major_model <- cm[cm$STC==1,]
SCT_major_model <- cm[cm$SCT==1,]
INDEP_major_model <- cm[cm$INDEP==1,]

BATF3 <- cm[cm$SNP_Label=="rs1156058",]
BATF3_indep_length <- length(which(BATF3$INDEP==1)) #47
BATF3_STC_length <- length(which(BATF3$STC==1))#0
BATF3_SCT_length <- length(which(BATF3$SCT==1)) #11
BATF3_SCT <- BATF3[BATF3$SCT==1,]
BATF3_indep <- BATF3[BATF3$INDEP==1,]

cm_STC <- cm[cm$STC==1,]
cm_SCT <- cm[cm$SCT==1,]

HMX2 <- cm[cm$SNP_Label=="rs705170",]
HMX2_STC <- HMX2[HMX2$STC==1,]
HMX2_SCT <- HMX2[HMX2$SCT==1,]
HMX2_INDEP <- HMX2[HMX2$INDEP==1,]
#23 HMX2, 23 SCT

PSMG1 <- cm[cm$SNP_Label=="rs2836950",]
PSMG1_indep_length <- length(which(PSMG1$INDEP==1)) #47
PSMG1_STC_length <- length(which(PSMG1$STC==1)) #47
PSMG1_SCT_length <- length(which(PSMG1$SCT==1)) #47


for_hist <- c(rep(1,STC),rep(2,SCT),rep(3,INDEP))

pdf("/data/jbryois/projects/Alspac/Fine_analysis_trans_eQTL_without_outliers/Independant_models/hist_number_per_model_normal_exp.pdf")
boxplot(STC_major_model[c(21,22,23)], col="skyblue2", main="Likelyhood of other models for STC", ylab="Relative Likelihood")
boxplot(SCT_major_model[c(21,22,23)], col="skyblue2", main="Likelyhood of other models for SCT", ylab="Relative Likelihood")
boxplot(INDEP_major_model[c(21,22,23)], col="skyblue2",main="Likelyhood of other models for INDEP", ylab="Relative Likelihood")
plot(sum, type='h')
pie(percentage, labels=c("STC", "SCT", "INDEP"))
pie(percentage_sum_best_model[c(2,3,1)], labels=c("STC", "SCT", "INDEP"))

hist(for_hist, col="skyblue2", xaxt='n', main="Histogram of Best model")
axis(1,at=c(1.1,1.9,2.9), labels=c("STC","SCT","INDEP"))
dev.off()
