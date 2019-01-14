library(pedigreemm)

#Get start and end line
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2) {
    print("Error, this script requires 1 arguments:")
    print("Number_of_covariate")
    stop("Exiting script")
}

start <- as.numeric(arguments[1])
end <- as.numeric(arguments[2])

#Function that run the mixed model and extracts the variance
#argument so that lmer does not fail due to internal control that is not relevent to pedigreem at the end (n obs=n random effects but here we ahve pedigree)

f <- function(x,gc,age,sex, smoke, date, lane, insert,id,family_id,ped,pheno){
model_full <- tryCatch(pedigreemm(as.numeric(x) ~  (age-1 | id) + gc  + age + sex + smoke + (1 |date) + (1|lane) + (1|insert) + (1 | famid) + pheno, REML=F,pedigree = list(id = ped),control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")),warning=function(err) NA)
model_restricted <- tryCatch(pedigreemm(as.numeric(x) ~ gc + age + sex + smoke + (1 |date) + (1|lane) + (1|insert) + (1 | famid) + pheno, REML=F,control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")),warning=function(err) NA)
if(is.na(model_full) | is.na(model_restricted)){return(NA)}
else{pvalue <- anova(model_full,model_restricted)[[8]][2]
return(pvalue)
}
}

#Load covariates files
covariates <- read.table("/data/project/processed_data/GAIT2/Info/Gait_Allcovariates.tab", header=T,sep="\t")
phenotypes <- read.table("/data/2backup/common/data/funpopgen_projects/GAIT2/Phenotypes/bloodCounts.txt", stringsAsFactors=F, header=T, sep="\t")
#phenotypes$LIMF+phenotypes$EOS+phenotypes$BAS+phenotypes$NEU+phenotypes$MON
phenotypes$SampleID <- paste("GAIT2_",phenotypes$ID,sep="")
phenotypes <- phenotypes[phenotypes$SampleID%in%covariates$sampleID,]
thrombosis <- read.csv("/data/2backup/common/data/funpopgen_projects/GAIT2/Phenotypes/thrombosis.txt",header=T)
thrombosis_matched <- thrombosis[match(phenotypes$ID,thrombosis$ID),]
thrombosis_matched$ID==phenotypes$ID
thrombosis_matched$ID==covariates$SubjectIDs
throm <- thrombosis_matched$Throm
phenotypes <- cbind(thrombosis_matched$Throm,phenotypes[c("BA","EO","HE","LI","MOt","NE")])

#Load expression scaled filtered and standard normalized
expression <- read.table("/data/project/processed_data/GAIT2/Quantification/Gait2_31082015_genes.raw.scaled.filter.rank",header=T)

#Check that covariates are in the same order
length(which(colnames(expression)[-c(1:4)]==covariates$sampleID))==length(colnames(expression)[-c(1:4)])
#TRUE

length(which(colnames(expression)[-c(1:4)]==phenotypes$SampleID))==length(colnames(expression)[-c(1:4)])
#TRUE

#Get Pedigree object
pedigree_info <- covariates[c("FA","MO","SubjectIDs","sampleID","FAMID","date","lane","GC_mean","INSERT_SIZE_MODE","AGE","SMOKIN","SEX")]

#individual missing -> get ind that don't have parents (required by pedigreemm) and puts parent to 0
ind_to_add_f <- unique(na.omit(as.matrix(pedigree_info[1])[!as.matrix(pedigree_info[1])%in%as.matrix(pedigree_info[3])]))
pedigree_info$FA[pedigree_info$FA%in%ind_to_add_f] <-0
pedigree_info$MO[pedigree_info$FA%in%ind_to_add_f] <-0

ind_to_add_m <- unique(na.omit(as.matrix(pedigree_info[2])[!as.matrix(pedigree_info[2])%in%as.matrix(pedigree_info[3])]))
pedigree_info$MO[pedigree_info$MO%in%ind_to_add_m] <-0
pedigree_info$FA[pedigree_info$MO%in%ind_to_add_m] <-0

ped <- pedigree(sire = as.integer(as.matrix(pedigree_info["FA"])),dam  = as.integer(as.matrix(pedigree_info["MO"])), label= as.character(as.matrix(pedigree_info["SubjectIDs"])))

#GET variables in right format for individuals and family
date <- as.factor(as.matrix(pedigree_info["date"]))
lane <- as.factor(as.matrix(pedigree_info["lane"]))
insert <- as.factor(as.matrix(pedigree_info["INSERT_SIZE_MODE"]))
smoke <- as.factor(as.matrix(pedigree_info["SMOKIN"]))
gc <- scale(as.numeric(as.matrix(pedigree_info["GC_mean"])),center=T,scale=T)
age <- scale(as.numeric(as.matrix(pedigree_info["AGE"])),center=T,scale=T)
#age <- as.numeric(as.matrix(pedigree_info["AGE"]))
sex <- as.numeric(as.factor(as.matrix(pedigree_info["SEX"])))-1
id <- as.character(as.matrix(pedigree_info["SubjectIDs"]))
famid <- as.factor(as.matrix(pedigree_info$FAMID))

######### FAMILY INTERACTION WITH CELL TYPE (probably not very powerfull)

#model_full <- pedigreemm(as.numeric(expression[-c(1:4)][start:end,]) ~  (sex-1 | id) + gc  + age + sex + smoke + (1 |date) + (1|lane) + (1|insert) + (1 | famid) + scale(phenotypes,center=T,scale=T), REML=F,pedigree = list(id = ped),control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore"))
#model_restricted <- tryCatch(pedigreemm(as.numeric(expression[-c(1:4)][start:end,]) ~ gc + age + sex + smoke + (1 |date) + (1|lane) + (1|insert) + (1 | famid) + scale(phenotypes,center=T,scale=T), REML=F,control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore")),warning=function(err) NA)
#anova(model_full,model_restricted)

#f <- function(x,gc,age,sex, smoke, date, lane, insert,id,family_id,ped,pheno){
family_interaction_2 <- apply(expression[-c(1:4)][start:end,],1,f,gc,age,sex,smoke,date,lane,insert,id,famid,ped,scale(phenotypes,center=T,scale=T))

final_variance_component <- cbind(expression[c(1:4)][start:end,],family_interaction_2)

write.table(final_variance_component, paste("/data/jbryois/projects/Ana_help/GAIT2_polygenic_interaction/interaction_individual_age_absolut",start,end,sep="_"),col.names=F,quote=F,sep="\t",row.names=F)