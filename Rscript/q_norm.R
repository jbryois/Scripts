library(preprocessCore)
#install.packages("preprocessCore")

#Example
#d <- data.frame(col_names=paste0("TAD_",1:1000),col1=rnorm(1000,0,1),col2=rnorm(1000,0,10),col3=rnorm(1000,10,15),col4=rnorm(1000,20,5),col5=rnorm(1000,100,5))
#write.table(d,"~/Desktop/test.txt",col.names=T,row.names=F,sep="\t",quote=F)

#For your data put your file name here
file = "~/Desktop/test.txt"

d <- read.table(file,header = T)

#get numerical data only 
d_num <- d[sapply(d,is.numeric)]

#get non_numeric_columns 
d_annot <- d[!sapply(d,is.numeric)]

#See distribution before quantile normalization
boxplot(d_num)

#Quantile normalise
qnorm = normalize.quantiles(as.matrix(d_num))

#See distribution after quantile normalization
boxplot(qnorm)

#add annotation back to quantile normalised data

d_qnorm <- cbind(d_annot,qnorm)
colnames(d_qnorm) <- colnames(d)

#Save file
write.table(d,paste0(file,".qnorm.txt"),col.names = T,row.names = F,sep = "\t",quote = F)
