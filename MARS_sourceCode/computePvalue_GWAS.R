rm(list=ls())
args=(commandArgs(TRUE))
if(length(args)<3){
	print("Usage: R CMD BATCH '--args -a=LRT_data -n=LRT_null -o=output [-f=MARS/fastMARS(0/1, default:0)] [-u=univariate threshold(default:5e-08)]' computePvalue_GWAS.R")
	exit()
}
UNIthreshold = 5e-08
fast = 0
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
	}else if(option=="-f"){
                fast = as.integer(value)
        }else if(option=="-u"){
		UNIthreshold = as.double(value)
	}else{
		print("Usage: R CMD BATCH '--args -a=LRT_data -n=LRT_null -o=output [-f=MARS/fastMARS(0/1, default:0)] [-u=univariate threshold(default:5e-08)]' computePvalue_GWAS.R")
	        exit()
	}
}

sink(paste(outPath))
cat("command : ", args, "\n\n")
null = as.matrix(read.table(nullPath))
alt = as.matrix(read.table(altPath))
simNum = dim(null)[[1]]
nullLRT = null[,2]
nullUNI = null[,3]
altLRT = alt[[1]]
altUNI = alt[[2]]
if(fast==0){
	LRT_pvalue = length(which(abs(nullLRT)>abs(altLRT)))/simNum
	UNI_pvalue = length(which(nullUNI<altUNI))/simNum
	quantile = length(which(nullUNI<UNIthreshold))
	LRTthreshold = sort(nullLRT)[simNum-quantile+1]
	cat("LRT_score\tunivariate_pvalue\tthreshold_LRT\tthreshold_UNI\tsignificance_LRT\tsignificance_UNI\n")
	cat(altLRT,"\t",altUNI, "\t", LRTthreshold,"\t", UNIthreshold,"\t",(altLRT>LRTthreshold),"\t",(altUNI<UNIthreshold),"\n")
}else{
	w = null[,1]
### pthreshold ###
        sum = 0;
        wsum =  0;
        for(i in c(1:simNum)){
                if(nullUNI[i]<UNIthreshold) wsum=wsum+w[i]
                sum=sum+w[i]
        }
        pvalue_threshold=wsum/sum;
### UNI p pvalue ###
        sum = 0;
        wsum =  0;
        for(i in c(1:simNum)){
                if(nullUNI[i]<altUNI) wsum=wsum+w[i]
                sum=sum+w[i]
        }
        UNI_pvalue = wsum/sum
### UNI LRT value ###	
	sum = 0;
        wsum =  0;
        for(i in c(1:simNum)){
                if(nullLRT[i]>altLRT) wsum=wsum+w[i]
                sum=sum+w[i]
        }
        LRT_pvalue = wsum/sum
	cat("LRT_pvalue\tunivariate_pvalue\tthreshold_pvalue\tthreshold_UNI\tsignificance_LRT\tsignificance_UNI\n")
	cat(LRT_pvalue,"\t",UNI_pvalue, "\t", pvalue_threshold,"\t", UNIthreshold,"\t",(pvalue_threshold > LRT_pvalue),"\t",(pvalue_threshold > UNI_pvalue),"\n")
}
sink();
