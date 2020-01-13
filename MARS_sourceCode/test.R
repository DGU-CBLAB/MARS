rm(list=ls())
args=(commandArgs(TRUE))
for(i in c(1:length(args))){
	param = strsplit(args[[i]], "=")
	print(param)
        option = param[[1]][1]
	value = param[[1]][2] 
        print(option)
	print(value)
	if(option=="-z") {print("z selected"); print(value);
	}else if(option=="-l"){print("l selected"); print(value)}
}
