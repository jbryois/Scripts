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

f <- function(x,date,lane,gc,insert,age,sex,ind_id,family_id,ped){
model <- pedigreemm(as.numeric(x) ~ (1 | ind_id) + (1 | date) + (1 | lane) + gc + (1 | insert) + age + sex + (1 | famid), pedigree = list(ind_id = ped),control=lmerControl(check.nobs.vs.nlev="ignore",check.nobs.vs.nRE="ignore"))
list <- list("A" = attr(VarCorr(model)$ind_id,"stddev")^2, "C" = attr(VarCorr(model)$famid,"stddev")^2,"E" = attr(VarCorr(model),"sc")^2, "date" = attr(VarCorr(model)$date,"stddev")^2, "lane"=attr(VarCorr(model)$lane,"stddev")^2, "insert"= attr(VarCorr(model)$insert,"stddev")^2)
return(list)
}

#Load covariates files
covariates <- read.table("/data/project/processed_data/GAIT2/Info/Gait_Allcovariates.tab", header=T)

#Load expression scaled filtered and standard normalized
#expression <- read.table(pipe("head -n 10 /data/project/processed_data/GAIT2/Quantification/Gait2_31082015_genes.raw.scaled.filter.rank"),header=T)
expression <- read.table("/data/project/processed_data/GAIT2/Quantification/Gait2_31082015_genes.raw.scaled.filter.rank",header=T)

#Check that covariates are in the same order
length(which(colnames(expression)[-c(1:4)]==covariates$sampleID))==length(colnames(expression)[-c(1:4)])

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
gc <- scale(as.numeric(as.matrix(pedigree_info["GC_mean"])),center=T,scale=T)
insert <- as.factor(as.matrix(pedigree_info["INSERT_SIZE_MODE"]))
age <- scale(as.numeric(as.matrix(pedigree_info["AGE"])),center=T,scale=T)
smoke <- as.factor(as.matrix(pedigree_info["SMOKIN"]))
sex <- as.factor(as.matrix(pedigree_info["SEX"]))
id <- as.character(as.matrix(pedigree_info["SubjectIDs"]))
famid <- as.factor(as.matrix(pedigree_info$FAMID))

variance_component <- apply(expression[-c(1:4)][start:end,],1,f,date,lane,gc,insert,age,sex,id,famid,ped)

A_list <- lapply(variance_component,function(x) x$A)
A_df <- data.frame(matrix(unlist(A_list)))
colnames(A_df) <- "A"

C_list <- lapply(variance_component,function(x) x$C)
C_df <- data.frame(matrix(unlist(C_list)))
colnames(C_df) <- "C"

E_list <- lapply(variance_component,function(x) x$E)
E_df <- data.frame(matrix(unlist(E_list)))
colnames(E_df) <- "E"

date_list <- lapply(variance_component,function(x) x$date)
date_df <- data.frame(matrix(unlist(date_list)))
colnames(date_df) <- "date"

lane_list <- lapply(variance_component,function(x) x$lane)
lane_df <- data.frame(matrix(unlist(lane_list)))
colnames(lane_df) <- "lane"

insert_list <- lapply(variance_component,function(x) x$insert)
insert_df <- data.frame(matrix(unlist(insert_list)))
colnames(insert_df) <- "insert"

h2 <- A_df/(A_df+C_df+E_df)
colnames(h2) <- "h2"
final_variance_component <- cbind(expression[c(1:4)][start:end,], A_df, C_df,E_df,h2,date_df,lane_df,insert_df)

write.table(final_variance_component, paste("/data/jbryois/projects/Ana_help/GAIT2_pedigreeMM/variance_component_gene_level",start,end,sep="_"),col.names=F,quote=F,sep="\t",row.names=F)