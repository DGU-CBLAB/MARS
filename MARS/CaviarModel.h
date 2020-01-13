#ifndef CAVIARMODEL_H
#define CAVIARMODEL_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>	
#include <armadillo>

#include "PostCal.h"

using namespace std;
using namespace arma;

 
class CaviarModel{
public:
	bool ldoption;
	bool alt;
	double rho;
	double NCP;
	double gamma;
	int snpCount;
	int totalCausalSNP;
	int seed;		
	int quantile;		
	string snpFile;
	string line;
	string data;
	double * sigma;
	double * sigma50;;
        double * stat;
	double * snp;		
	double * snp50;
        char * pcausalSet;
        int * rank;
	int subSize;
	int simNum;
	PostCal * post;
	string ldFile;
        string zFile;
        string outputFileName;
        string geneMapFile;	
	mat sigmaMatrix;	
	int * topIndex; //for null
        double * LRTlist;       
        double * Plist;         
        double * wlist;
	int sampleSize;

	CaviarModel(bool ldoption, string snpFile, string ldFile, string zFile, string outputFileName, int sampleSize, int simNum, int totalCausalSNP, bool alt, double NCP, double rho, double gamma=0.01) {
		int tmpSize = 0;
		subSize = 50;
		this->sampleSize = sampleSize;
		this->simNum = simNum;
		this->alt = alt;
		this->ldoption = ldoption;
		this->NCP = NCP;
		this->rho = rho;
		this->gamma = gamma;
		this->ldFile = ldFile;
		this->snpFile = snpFile;
		this->zFile  = zFile;
		this->outputFileName = outputFileName;
		this->totalCausalSNP = totalCausalSNP;
		if(ldoption){
			fileSize(ldFile, tmpSize);
			snpCount   = (int)sqrt(tmpSize);
		}
		else{
			fileSize_snp(snpFile, tmpSize);
			snpCount   = (int)(tmpSize/(sampleSize+1));
		}
		cout<<"snpCount="<<snpCount<<endl;
         	sigma      = new double[snpCount * snpCount];
		sigma50   = new double[subSize * subSize];
		stat       = new double[snpCount];				
		snp      = new double[snpCount * sampleSize];	
		snp50     = new double[subSize*sampleSize];
		pcausalSet = new char[snpCount];
		rank       = new int[snpCount];
		sigmaMatrix = mat (snpCount, snpCount);
		if(alt){
			if(ldoption){
				importData(ldFile, sigma);
			}
			else{
				loadSNP(snpFile, snp, snpCount, sampleSize);	
				makeSubSigma(snp, sigma, subSize, sampleSize);
			}
//		makeSigmaPositiveSemiDefinite(sigma, snpCount);
		//importStat(zFile, stat);
			importStat2(zFile, stat, subSize);
			sigmaMatrix = mat (snpCount, snpCount);
			for(int i = 0; i < snpCount; i++) {	
				for (int j = 0; j < snpCount; j++){
                                	sigmaMatrix(i,j) = sigma[i*snpCount+j];
			//	if(j==0) cout<< sigmaMatrix(i,j)<<" ";
				}
                	}
			post = new PostCal(sigmaMatrix, stat, snpCount, totalCausalSNP, gamma);       					
			/*for(int i=0; i<5 ; i++){//test!!! 
				for(int j=0; j<5 ; j++){
					cout<< sigmaMatrix(i,j)<<" ";
				}cout<<endl;
			}
			cout<<endl;
			for(int i=0; i<5 ; i++){//test!!!
                                        cout<< stat[i]<<" ";
                        }*/
		}else{
			cout<<"Number of simulations="<<simNum<<endl;
			loadSNP(snpFile, snp, snpCount, sampleSize);
			LRTlist    = new double[simNum];  
        	        Plist      = new double[simNum]; 
	                wlist      = new double[simNum];	
			topIndex = new int[subSize];
		}
	}
	void run() {
		if(alt){
        		post->computePvalue(stat, NCP);
			post->computeLRT(stat, NCP);
			post->saveLRT();
		}else{
			ifstream fin(zFile.c_str(), std::ifstream::in);
			for(int i=0; i<simNum; i++){
				stat       = new double[snpCount];
				getline(fin, line);
				istringstream iss(line);
				iss>>data;
				wlist[i] = atof(data.c_str());
				for(int j = 0; j<subSize; j++){
					iss >> data;
					stat[j] = atof(data.c_str());
					iss >> data;
					topIndex[j] = atoi(data.c_str());	
				}
				for(int j=0;j<subSize;j++){
                                	for(int k = 0;k<sampleSize; k++){
                                        	snp50[j*sampleSize+k]=snp[topIndex[j]*sampleSize+k];
                                	}
                        	}
				makeSubSigma(snp50, sigma50, subSize, sampleSize);
				makeSigmaPositiveSemiDefinite(sigma50, subSize);
	                        for(int j = 0; j < subSize; j++) {
        	                        for (int k = 0; k < subSize; k++)
                	                        sigmaMatrix(j,k) = sigma50[j*subSize+k];
                        	}
                        	post = new PostCal(sigmaMatrix, stat, subSize, totalCausalSNP, gamma);
                        	post->computeLRT(stat, NCP, subSize);//JOO add
                        	post->computePvalue(stat, subSize);
                        	post->saveLRT();//string(outputFileName)+"_LRT");
                        	LRTlist[i] = post->LRTscore_;
                        	Plist[i] = post->pvalue_;
                                delete[] stat;
                                delete post;			
			}
			fin.close();		
		}
	}
	void finishUp(string fileName){
		ofstream outfile(fileName.c_str(), ios::out | ios::app);
		if(alt){
			cout<< post->LRTscore_ << "\t"<<post->pvalue_<<endl;
			outfile << post->LRTscore_ << "\t"<<post->pvalue_<<endl;
		}else{
			for(int i=0; i<simNum; i++){
                        	outfile << wlist[i]<<"\t"<<LRTlist[i] <<"\t"<<Plist[i]<< endl;
                        }
		}
		outfile.close();
	}

	~CaviarModel() { 
		delete[] sigma; 
		delete[] sigma50; 
		delete[] snp; 
		delete[] snp50;
		if(!alt){
			delete[] LRTlist;
        	        delete[] wlist;
                	delete[] Plist;
			delete[] topIndex;
		}
		else delete post;
	}

};
 
#endif
