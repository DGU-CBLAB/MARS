rm(list=ls())
library(mvtnorm)
args=(commandArgs(TRUE))
if(length(args)<4){
    print("Usage: R CMD BATCH '--args -g=genotypePath -s=statPath -o=output_genotype -u=output_stat [-t=topNum(default:50)]' MARS_alt.R")
    quit();
}
topNum = 50
for(i in c(1:length(args))){
        param = strsplit(args[[i]], "=")
        option = param[[1]][1];
        value = param[[1]][2];
	if(option=="-g") {
                genotypePath = value
        }else if(option=="-s"){
                statPath = value
        }else if(option=="-o"){
                outputPath_geno = value
        }else if(option=="-u"){
                outputPath_stat = value
	}else if(option=="-t"){
                topNum = value
        }else{
                print("Usage: R CMD BATCH '--args -g=genotypePath -s=statPath -o=output_genotype -u=output_stat [topNum(default:50)]' MARS_alt.R")
                quit()
        }
}

x = as.matrix(read.table(genotypePath));
S = as.matrix(read.table(statPath));

I=(order(abs(S), decreasing=TRUE))[1:topNum]
write.table(t(S[I]),outputPath_stat, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(x[I,],outputPath_geno, row.names=FALSE, col.names=FALSE, quote=FALSE)


