rm(list=ls())
library(Matrix)
args=(commandArgs(TRUE))
if(length(args)<2){
    print("Usage: R CMD BATCH '--args -g=genotypePath -o=outputPath' generateLD.R") 
    quit()
}
topNum = 50
for(i in c(1:length(args))){
        param = strsplit(args[[i]], "=")
        option = param[[1]][1];
        value = param[[1]][2];
        if(option=="-g") {
                genotypePath = value
        }else if(option=="-o"){
                outputPath = value
        }else{
                print("Usage: R CMD BATCH '--args -g=genotypePath -o=outputPath' generateLD.R")
                quit()
        }
}
x=as.matrix(read.table(genotypePath))
n = dim(x)[2];# sample num
x = x[,2:n];
d <-dim(x)
x<-as.numeric(x)
dim(x)<-d
ld = cor(t(x));# m by m, correlation between SNPS in terms of samples
ldPD = data.matrix(nearPD(ld)$mat)
write.table(ldPD,outputPath,row.names=FALSE, col.names=FALSE, quote=FALSE)
