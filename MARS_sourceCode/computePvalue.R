rm(list=ls())
args=(commandArgs(TRUE))
if(length(args)<3){
	print("Usage: R CMD BATCH '--args -a=LRT_data -n=LRT_null -o=output [-t=threshold(default:0.05)] computePvalue.R")
	exit()
}
pThreshold = 0.05
for(i in c(1:length(args))){
        param = strsplit(args[[i]], "=")
        option = param[[1]][1]
        value = param[[1]][2]
        if(option=="-a") {
        	altPath = value
	}else if(option=="-n"){
		nullPath = value
	}else if(option=="-o"){
		outPath = value
        }else if(option=="-t"){
		pThreshold = as.double(value)
	}else{
		print("Usage: R CMD BATCH '--args -a=LRT_data -n=LRT_null -o=output [-t=threshold(default:0.05)] computePvalue.R")
	        exit()
	}
}
sink(paste(outPath))
cat("command: ",args,"\n\n")
cat("pvalue_LRT\tpvalue_UNI\tpThreshold\tsignificance_LRT\tsignificance_UNI\n")
null = as.matrix(read.table(nullPath))
alt = as.matrix(read.table(altPath))
simNum = dim(null)[[1]]
nullLRT = null[,2]
nullUNI = null[,3]
altLRT = alt[[1]]
altUNI = alt[[2]]
LRT_pvalue = length(which(abs(nullLRT)>abs(altLRT)))/simNum
UNI_pvalue = length(which(nullUNI<altUNI))/simNum
cat(LRT_pvalue,"\t",UNI_pvalue, "\t", pThreshold,"\t",(LRT_pvalue<=pThreshold),"\t",(UNI_pvalue<=pThreshold),"\n")
sink();
