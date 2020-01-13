#include <vector>
#include <algorithm>
#include <set>
#include <iostream>
#include <armadillo>

#include "Util.h"
#include "PostCal.h"

using namespace arma;


void printGSLPrint(mat &A, int row, int col) {
	for(int i = 0; i < row; i++) {
		for(int j = 0; j < col; j++)
			printf("%g ", A(i, j));
		printf("\n");
	}	
}

string PostCal::convertConfig2String(int * config, int size) {
	string result = "0";
	for(int i = 0; i < size; i++)
		if(config[i]==1)
			result+= "_" + convertInt(i);
	return result;
}

// We compute dmvnorm(Zcc, mean=rep(0,nrow(Rcc)), Rcc + Rcc %*% Rcc) / dmvnorm(Zcc, rep(0, nrow(Rcc)), Rcc))
// // togheter to avoid numerical over flow
double PostCal::fracdmvnorm(mat Z, mat mean, mat R, mat diagC, double NCP) {
        mat newR = R + R * diagC  * R;
        mat ZcenterMean = Z - mean;
        mat res1 = trans(ZcenterMean) * inv(R) * (ZcenterMean);
        mat res2 = trans(ZcenterMean) * inv(newR) *  (ZcenterMean);

        double v1 = res1(0,0)/2-res2(0,0)/2-baseValue/2;
        return(exp(v1)/sqrt(det(newR))* sqrt(det(R)));
}

double PostCal::dmvnorm(mat Z, mat mean, mat R) {
        mat ZcenterMean = Z - mean;
        mat res = trans(ZcenterMean) * inv(R) * (ZcenterMean);
        double v1 = res(0,0);
        double v2 = log(sqrt(det(R)));
        return (exp(-v1/2-v2));
}

// cc=causal SNPs
// Rcc = LD of causal SNPs
// Zcc = Z-score of causal SNPs
// dmvnorm(Zcc, mean=rep(0,nrow(Rcc)), Rcc + Rcc %*% Rcc) / dmvnorm(Zcc, rep(0, nrow(Rcc)), Rcc))
//
double PostCal::fastLikelihood(int * configure, double * stat, double NCP) {
	int causalCount = 0;
	vector <int> causalIndex;

	for(int i = 0; i < snpCount; i++) {
		causalCount += configure[i];
		if(configure[i] == 1)
			causalIndex.push_back(i);
	}
	
	if (causalCount == 0) {
		int maxVal = 0;
		for(int i = 0; i < snpCount; i++) {
			if (maxVal < abs(stat[i]))
				maxVal = stat[i];
		}
		baseValue = maxVal * maxVal;
	}

	mat Rcc(causalCount, causalCount, fill::zeros);
	mat Zcc(causalCount, 1, fill::zeros);
	mat mean(causalCount, 1, fill::zeros);
	mat diagC(causalCount, causalCount, fill::zeros);
	for (int i = 0; i < causalCount; i++){
		for(int j = 0; j < causalCount; j++) {
			Rcc(i,j) = sigmaMatrix(causalIndex[i], causalIndex[j]);
		}
		Zcc(i,0) = stat[causalIndex[i]];
		diagC(i,i) = NCP;
	}
	
	while (det(Rcc) <= 0.01) {	
		mat toAdd(causalCount, causalCount);
		toAdd.eye();
		Rcc = Rcc + 0.1 * toAdd;
	}
	return fracdmvnorm(Zcc, mean, Rcc, diagC, NCP);
}

/*double PostCal::likelihood(int * configure, double * stat, double NCP) {//JOO remove start
	int causalCount = 0;
	int index_C = 0;
        double matDet = 0;
	double res    = 0;

	for(int i = 0; i < snpCount; i++) 
		causalCount += configure[i];
	if(causalCount == 0){
		mat tmpResultMatrix1N = statMatrixtTran * invSigmaMatrix;
		mat tmpResultMatrix11 = tmpResultMatrix1N * statMatrix;
		res = tmpResultMatrix11(0,0);	
		baseValue = res;
		matDet = sigmaDet;
		res = res - baseValue;
		return( exp(-res/2)/sqrt(abs(matDet)) );
	}
	mat U(snpCount, causalCount, fill::zeros);
	mat V(causalCount, snpCount, fill::zeros);
	mat VU(causalCount, causalCount, fill::zeros);
		
	for(int i = 0; i < snpCount; i++) {
                if (configure[i] == 0)	continue;
                else {
                        for(int j = 0; j < snpCount; j++) 
                                U(j, index_C) = sigmaMatrix(j,i);
			V(index_C, i) = NCP;
                        index_C++;
                }
        }
	VU = V * U;
	mat I_AA   = mat(snpCount, snpCount, fill::eye);
	mat tmp_CC = mat(causalCount, causalCount, fill::eye)+ VU;
	matDet = det(tmp_CC) * sigmaDet;
	mat tmp_AA = invSigmaMatrix - (invSigmaMatrix * U) * pinv(tmp_CC) * V ;
	//tmp_AA     = invSigmaMatrix * tmp_AA;
	mat tmpResultMatrix1N = statMatrixtTran * tmp_AA;
        mat tmpResultMatrix11 = tmpResultMatrix1N * statMatrix;
        res = tmpResultMatrix11(0,0);  

	res = res - baseValue;
	if(matDet==0) {
		cout << "Error the matrix is singular and we fail to fix it." << endl;
		exit(0);
	}
	
	//	We compute the log of -res/2-log(det) to see if it is too big or not. 
	//	In the case it is too big we just make it a MAX value.
	
	double tmplogDet = log(sqrt(abs(matDet)));
	double tmpFinalRes = -res/2 - tmplogDet;
	if(tmpFinalRes > 700) 
		return(exp(700));
	return( exp(-res/2)/sqrt(abs(matDet)) );	
}*///JOO remove end

int PostCal::nextBinary(int * data, int size) {
	int i = 0;
	int total_one = 0;	
	int index = size-1;
        int one_countinus_in_end = 0;

        while(index >= 0 && data[index] == 1) {
                index = index - 1;
                one_countinus_in_end = one_countinus_in_end + 1;
	}
	if(index >= 0) {
        	while(index >= 0 && data[index] == 0) {
               	 index = index - 1;	
		}
	}
        if(index == -1) {
                while(i <  one_countinus_in_end+1 && i < size) {
                        data[i] = 1;
                        i=i+1;
		}
                i = 0;
                while(i < size-one_countinus_in_end-1) {
                        data[i+one_countinus_in_end+1] = 0;
                        i=i+1;
		}
	}
        else if(one_countinus_in_end == 0) {
                data[index] = 0;
                data[index+1] = 1;
	} else {
                data[index] = 0;
                while(i < one_countinus_in_end + 1) {
                        data[i+index+1] = 1;
			if(i+index+1 >= size)
				printf("ERROR3 %d\n", i+index+1);
                        i=i+1;
		}
                i = 0;
                while(i < size - index - one_countinus_in_end - 2) {
                        data[i+index+one_countinus_in_end+2] = 0;
			if(i+index+one_countinus_in_end+2 >= size) {
				printf("ERROR4 %d\n", i+index+one_countinus_in_end+2);
			}
                        i=i+1;
		}
	}
	i = 0;
	total_one = 0;
	for(i = 0; i < size; i++)
		if(data[i] == 1)
			total_one = total_one + 1;
	
	return(total_one);		
}

double PostCal::computeTotalLikelihood(double * stat, double NCP) {	
	int num = 0;
	double sumLikelihood = 0;
	double tmp_likelihood = 0;
	long int total_iteration = 0 ;
	int * configure = (int *) malloc (snpCount * sizeof(int *)); // original data	

	for(long int i = 0; i <= maxCausalSNP; i++)
		total_iteration = total_iteration + nCr(snpCount, i);
	cout << snpCount << endl;
	cout << "Max Causal=" << maxCausalSNP << endl;
	cout << "Total="      << total_iteration << endl;
	for(long int i = 0; i < snpCount; i++) 
		configure[i] = 0;
	for(long int i = 0; i < total_iteration; i++) {
                tmp_likelihood = fastLikelihood(configure, stat, NCP) * (pow(gamma, num))*(pow(1-gamma, snpCount-num));
                sumLikelihood += tmp_likelihood;
		for(int j = 0; j < snpCount; j++) {
                        postValues[j] = postValues[j] + tmp_likelihood * configure[j];
		}
		histValues[num] = histValues[num] + tmp_likelihood;
                num = nextBinary(configure, snpCount);
       		if(i % 100000 == 0)
			cout << i << " "  << sumLikelihood << endl;
	}
	for(int i = 0; i <= maxCausalSNP; i++)
		histValues[i] = histValues[i]/sumLikelihood;
        free(configure);
        return(sumLikelihood);
}

bool PostCal::validConfigutation(int * configure, char * pcausalSet) {
	for(int i = 0; i < snpCount; i++){
		if(configure[i] == 1 && pcausalSet[i] == '0')
			return false;
	}
	return true;	
}

/*
 * This is a auxilary function used to generate all possible causal set that 
 * are selected in the p-causal set
*/
void PostCal::computeALLCausalSetConfiguration(double * stat, double NCP, char * pcausalSet, string outputFileName) {
	int num = 0;
        double sumLikelihood = 0;
        double tmp_likelihood = 0;
        long int total_iteration = 0 ;
        int * configure = (int *) malloc (snpCount * sizeof(int *)); // original data   

        for(long int i = 0; i <= maxCausalSNP; i++)
                total_iteration = total_iteration + nCr(snpCount, i);
        for(long int i = 0; i < snpCount; i++)
                configure[i] = 0;
        for(long int i = 0; i < total_iteration; i++) {
		if (validConfigutation(configure, pcausalSet)) {
                	tmp_likelihood = fastLikelihood(configure, stat, NCP) * (pow(gamma, num))*(pow(1-gamma, snpCount-num));
			//printVector(configure, snpCount);
			//cout << " " << tmp_likelihood << endl;
			exportVector2File(outputFileName, configure, snpCount);
			export2File(outputFileName, tmp_likelihood);
		}
		num = nextBinary(configure, snpCount);
	}
}

/*
	stat is the z-scpres
	sigma is the correaltion matrix
	G is the map between snp and the gene (snp, gene)
*/
double PostCal::findOptimalSetGreedy(double * stat, double NCP, char * pcausalSet, int *rank,  double inputRho, string outputFileName) {
	int index = 0;
        double rho = 0;
        double total_post = 0;

        totalLikeLihood = computeTotalLikelihood(stat, NCP);
	
	export2File(outputFileName+".log", totalLikeLihood); //Output the total likelihood to the log File
	for(int i = 0; i < snpCount; i++)
		total_post += postValues[i];
	printf("Total Likelihood= %e SNP=%d \n", total_post, snpCount);
	
        std::vector<data> items;
        std::set<int>::iterator it;
	//output the poster to files
        for(int i = 0; i < snpCount; i++) {
             //printf("%d==>%e ",i, postValues[i]/total_likelihood);
             items.push_back(data(postValues[i]/total_post, i, 0));
        }
        printf("\n");
        std::sort(items.begin(), items.end(), by_number());
        for(int i = 0; i < snpCount; i++)
                rank[i] = items[i].index1;

        for(int i = 0; i < snpCount; i++)
                pcausalSet[i] = '0';
        do{
                rho += postValues[rank[index]]/total_post;
                pcausalSet[rank[index]] = '1';
                printf("%d %e\n", rank[index], rho);
                index++;
        } while( rho < inputRho);

        printf("\n");
	return(0);
}
double PostCal::computePvalue(double * stat, int snpCount){
//      cout<<"***INSIDE POST***"<<endl;
                        /*cout<<"snpCount="<<snpCount<<endl;
                cout<<"***INSIDE POST ,Top 50 Stat=***"<<endl;
                        for(int kk= 0;kk<snpCount;kk++) cout<<stat[kk]<<" ";
                        cout<<" " <<endl;*/
        double max;
        max = findAbsMax(stat,snpCount);
        pvalue = 2*normalCDF(-max);
//      cout<<"max= "<<max<<" pvalue= "<<pvalue<<endl;
        return pvalue;
}
//JOO add start
double PostCal::findAbsMax(double* v, int m){
        double max = 0;
        double val=0;
        for (int i=0;i<m; i++) {
                val = fabs(v[i]);
                if (max < val) max = val;
        }
        return max;
}
double PostCal::normalCDF(double value)
{
   return 0.5 * erfc(-value * M_SQRT1_2);
}
double PostCal::computeLRT(double * stat, double NCP) {
        int num = 0; //int total = 0;
        double sumLikelihood = 0;
        double allZero_likelihood = 0;//define allZero_likelihood
        double tmp_likelihood = 0;
        long int total_iteration = 0 ;
        int * configure = (int *) malloc (snpCount * sizeof(int *)); // original data

        for(long int i = 0; i <= maxCausalSNP; i++)
                total_iteration = total_iteration + nCr(snpCount, i);
        for(long int i = 0; i < snpCount; i++)
                configure[i] = 0;
        for(long int i = 0; i < total_iteration; i++) {
                tmp_likelihood = fastLikelihood(configure, stat, NCP) * (pow(gamma, num))*(pow(1-gamma, snpCount-num));
                if(i == 0){ allZero_likelihood =tmp_likelihood; }//save allzero likelihood
                sumLikelihood += tmp_likelihood;
                for(int j = 0; j < snpCount; j++) {
                        postValues[j] = postValues[j] + tmp_likelihood * configure[j];
                }
                histValues[num] = histValues[num] + tmp_likelihood;
                num = nextBinary(configure, snpCount);
        }//cout<<"here="<<sigmaMatrix(1,0)<<endl;
        for(int i = 0; i <= maxCausalSNP; i++)
                histValues[i] = histValues[i]/sumLikelihood;
        LRTscore = (sumLikelihood-allZero_likelihood)/allZero_likelihood;//compute LRTscore
        free(configure);
	//for(int i=0; i<snpCount; i++) total = total+stat[i];
	//if(total==0) LRTscore= 10000000; 
        return(LRTscore);
}
double PostCal::computeLRT(double * stat, double NCP, int subSize) {
        int num = 0; //int total = 0;
        double sumLikelihood = 0;
        double allZero_likelihood = 0;//define allZero_likelihood
        double tmp_likelihood = 0;
        long int total_iteration = 0 ;
        int * configure = (int *) malloc (subSize * sizeof(int *)); // original data

        for(long int i = 0; i <= maxCausalSNP; i++)
                total_iteration = total_iteration + nCr(subSize, i);
        for(long int i = 0; i < subSize; i++)
                configure[i] = 0;
        for(long int i = 0; i < total_iteration; i++) {
                tmp_likelihood = fastLikelihood(configure, stat, NCP) * (pow(gamma, num))*(pow(1-gamma, subSize-num));
                if(i == 0){ allZero_likelihood =tmp_likelihood; }//save allzero likelihood
                sumLikelihood += tmp_likelihood;
                for(int j = 0; j < subSize; j++) {
                        postValues[j] = postValues[j] + tmp_likelihood * configure[j];
                }
                histValues[num] = histValues[num] + tmp_likelihood;
                num = nextBinary(configure, subSize);
        }//cout<<"here="<<sigmaMatrix(1,0)<<endl;
        for(int i = 0; i <= maxCausalSNP; i++)
                histValues[i] = histValues[i]/sumLikelihood;
        LRTscore = (sumLikelihood-allZero_likelihood)/allZero_likelihood;//compute LRTscore
        free(configure);
        return(LRTscore);
}
//JOO add end
