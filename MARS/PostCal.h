#ifndef POSTCAL_H
#define POSTCAL_H

#include <iostream>
#include <fstream>
#include <map>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <armadillo>

using namespace std;
using namespace arma;

void printGSLPrint(mat A, int row, int col);
 
class PostCal{

private:
	double gamma;		// the probability of SNP being causal
	double * postValues;	//the posterior value for each SNP being causal
//	double * sigma;		//the LD matrix //JOO remove
	double * histValues;	//the probability of the number of causal SNPs, we make the histogram of the causal SNPs
	int snpCount;		//total number of variants (SNP) in a locus
	int maxCausalSNP;	//maximum number of causal variants to consider in a locus
	double sigmaDet;	//determinie of matrix
	double totalLikeLihood; //Compute the total likelihood of all causal status (by likelihood we use prior)
	double baseValue;
	double LRTscore; 	//JOO add
	double pvalue;
	mat sigmaMatrix;	
	mat invSigmaMatrix;
	mat statMatrix;
        mat statMatrixtTran;
public:
	double LRTscore_;        //JOO add
	double pvalue_;
	 PostCal(mat sigmaMatrix, double * stat, int snpCount, int maxCausalSNP, double gamma) {//JOO change
		baseValue = 0;
		LRTscore = 0;	//JOO add
		LRTscore_=0;	//JOO add
		pvalue= 0;
		pvalue_ = 0;
		this->gamma = gamma;
		this-> snpCount = snpCount;
		this-> maxCausalSNP = maxCausalSNP;
//                this->sigma = new double[snpCount * snpCount];//JOO remove
		this-> postValues = new double [snpCount];
		this-> histValues = new double [maxCausalSNP+1];              
//		this->invSigmaMatrix = invSigmaMatrix;	//JOO add
//		this->sigmaDet = sigmaDet;		//JOO add
		this->sigmaMatrix = sigmaMatrix;	//JOO add
		statMatrix                 = mat (snpCount, 1);
		statMatrixtTran            = mat (1, snpCount);
		//sigmaMatrix         	   = mat (snpCount, snpCount);// JOO remove
	
//		for(int i = 0; i < snpCount*snpCount; i++)//JOO remove
//			this->sigma[i] = sigma[i];	//JOO remove
		for(int i = 0; i < snpCount; i++)
                        this->postValues[i] = 0;
		for(int i= 0; i <= maxCausalSNP;i++)
			this->histValues[i] = 0;
		for(int i = 0; i < snpCount; i++) {
                	statMatrix(i,0) = stat[i];
        	        statMatrixtTran(0,i) = stat[i];
	        }
		
/*		for(int i = 0; i < snpCount; i++) {	//JOO remove start	
                	for (int j = 0; j < snpCount; j++)
                       		sigmaMatrix(i,j) = sigma[i*snpCount+j];
       		}
		invSigmaMatrix = inv(sigmaMatrix);
		sigmaDet       = det(sigmaMatrix);	//JOO remove end
		*///cout << sigmaMatrix(0,2)<<"\t"<<invSigmaMatrix(0,0)<<"\t"<<invSigmaMatrix(3,4)<<"\t"<<sigmaDet<<endl;
	
	}
        ~PostCal() {
		delete [] histValues;
		delete [] postValues;
//                delete [] sigma;//JOO remove
	}

	bool validConfigutation(int * configure, char * pcausalSet);
	void computeALLCausalSetConfiguration(double * stat, double NCP, char * pcausalSet, string outputFileName);


	double dmvnorm(mat Z, mat mean, mat R);
        double fracdmvnorm(mat Z, mat mean, mat R, mat diagC, double NCP);

        double fastLikelihood(int * configure, double * stat, double NCP);
//	double likelihood(int * configure, double * stat, double NCP) ;//JOO remove
	int nextBinary(int * data, int size) ;
	double computeTotalLikelihood(double * stat, double NCP) ;	
	double findOptimalSetGreedy(double * stat, double NCP, char * pcausalSet, int *rank,  double inputRho, string outputFileName);
        double findAbsMax(double* v, int m);			//JOO add
        double normalCDF(double value);				//JOO add
        double computeLRT(double * stat, double NCP) ;		//JOO add
	double computeLRT(double * stat, double NCP, int subSize);
	double computePvalue(double * stat, int snpCount);
	string convertConfig2String(int * config, int size);
	void printHist2File(string fileName) {
		exportVector2File(fileName, histValues, maxCausalSNP+1);
	}

	void saveLRT(){
                LRTscore_ = LRTscore;
		pvalue_ = pvalue;
	}
	void printLRT2File(string fileName) {	//JOO add start
                ofstream outfile(fileName.c_str(), ios::out | ios::app);
                outfile << LRTscore << "\t"<<pvalue<<endl;
	}					//JOO add end
};
 
#endif
