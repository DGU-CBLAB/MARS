rm(list=ls())
library(mvtnorm)
args=(commandArgs(TRUE))
if(length(args)<3){
	print("Usage: R CMD BATCH '--args [-n=number_of_simulations] [-g=genotype] [-o=output] [-f=MARS/fastMARS(0/1, default:0)] [-t=topNum]' MARS_NULL.R")
	quit()
}
fast=0
topNum = 50
for(i in c(1:length(args))){
        param = strsplit(args[[i]], "=")
        option = param[[1]][1];
        value = param[[1]][2];
        if(option=="-n") {
                simNum=as.integer(value)
        }else if(option=="-g"){
		genotypePath = value
	}else if(option=="-o"){
                outputPath = value
        }else if(option=="-f"){
                fast = value
        }else if(option=="-t"){
                topNum = as.integer(value)
	}else{
		print("Usage: R CMD BATCH '--args [-n=number_of_simulations] [-g=genotype] [-o=output] [-f=MARS/fastMARS(0/1, default:0)] [-t=topNum]' MARS_NULL.R")
	        quit()
	}
}
top50<-function(fS,fR, topNum){
         I50=(order(abs(fS), decreasing=TRUE))[1:topNum]
         S50=fS[I50]
         fR=c(rbind(S50,(I50-1)))
         return(fR)
}

start_time <- Sys.time()
R = matrix(0,nrow=simNum, ncol=(topNum*2))
x = as.matrix(read.table(genotypePath));
m = dim(x)[1]
n = dim(x)[2];
x = x[,2:n];
n = n-1;
d <-dim(x)
x<-as.numeric(x)
dim(x)<-d
x=t(x)
xs<-scale(x,center=TRUE, scale=TRUE)
I = diag(n)
mu = mat.or.vec(n,1);
if(fast==0){
	print("Genering null for MARS")
	Sall=rmvnorm(simNum, mu, I)#simNum, mu, I)
	Sall_new = (t(xs)/sqrt(n))%*%t(Sall)
	R = t(apply(Sall_new,2,top50, fR=R, topNum=topNum))
	write.table(cbind(matrix(0,nrow=simNum),R), outputPath, row.names=FALSE, col.names=FALSE, quote=FALSE)
	rm(R, Sall_new, Sall, mu, x, xs, I, d, n, m)
}else{
	print("Genering null for fastMARS")
	I2 = sqrt(2)*I
	Sall=rmvnorm(simNum, mu, I2)#simNum, mu, I)
	W = dmvnorm(Sall, mu, I,log=TRUE)/dmvnorm(Sall, mu, I2, log=TRUE)
	Sall_new = (t(xs)/sqrt(n))%*%t(Sall)
	R = t(apply(Sall_new,2,top50, fR=R, topNum=topNum))
	write.table(cbind(W,R), outputPath, row.names=FALSE, col.names=FALSE, quote=FALSE)
	rm(R, Sall_new, Sall, mu, x, xs, I, I2, W, d, n, m)	
}
end_time <- Sys.time()
cat("Running time", (end_time - start_time), "\n")

