#ifndef UTIL_H
#define UTIL_H

#include <cmath>
#include <string>
#include <vector>

using namespace std;

	struct data {
                data(double num, int ind1, int ind2) {
                        number = num;
                        index1 = ind1;
                        index2 = ind2;
                }
                double number;
                int index1;
                int index2;
        };

        struct by_number {
            bool operator()(data const &left, data const &right) {
                return abs(left.number) > abs(right.number);
            }
        };



	string convertInt(int number);
	long int fact(int n) ;
	void copyConfigure(double *dest, double *src, int size) ;
	double min(double a, double b) ;
	long int nCr(int n, int r) ;
	void printVector(char * data, int size) ;
	void printVector(int * data, int size) ;
	void printVector(double * data, int size) ;
	void diffVector(double * data1, double * data2, int size, double * result) ;
	void sumVector(double * data1, double * data2, int size, double * result) ;
	double multVector(double * data1, double * data2, int size) ;
	void dotVector(double * data1, double * data2, int size, double * result) ;
	void multVectorMatrix(double *vector, double * matrix, int size, double * result) ;
	void fileSize(string fileName, int & size);
	void fileSize_snp(string fileName, int & size);
	void importData(string fileName, double * vector);
	void importData(string fileName, int * vector);
	void importStat(string fileName, double * vector);
	void importStat2(string fileName, double * vector, int subSize);
	void importDataSecondColumn(string fileName, double * vector);
	void importDataNthColumn(string fileName, double * vector, int colNum, int ignore=0);
	void importDataFirstColumn(string fileName, string * list, int ignore=0);
	void loadSNP(string fileName, double* snp, int snpCount, int sampleSize);
	double correlationCoefficient(double* X, double* Y, int n);
	void makeSubSigma(double* snp_50, double* sigma_50, int subSize, int sampleSize);
	void rmvnorm(double * mean, double * sigma, int size, double * results);
	void resetVector(char *data, int size);
	void resetVector(int *data, int size);
	void resetVector(double *data, int size);
	void generateMean(int * causalSNP, double * sigma, int snpCount, double * result);
	void exportVector2File(string fileName, char * data, int size);
	void exportVector2File(string fileName, int * data, int size);
	void exportVector2File(string fileName, double * data, int size);
	void export2File(string fileName, int data);
	void export2File(string fileName, double data);
	void matrixMul(int * aData, int * bData, int * cData, int row1, int col1, int row2, int col2);
	int snp2Gene(int * G, int snpId, int snpCount, int geneCount);
	void setIdentitymatrix(int * G, int snpCount, int geneCount);
	void makeSigmaPositiveSemiDefinite(double * sigma, int size);
	void rmvnormC(double * mean, double * sigma, int seed, int size, double* results);
	double* generateMu(int snpCount, double NCP);
	double* generateStat(double * sigma, int m, int seed, int* quantile, double threshold, double mu[]);
	double normalCFD(double value);
        double findAbsMax(double* v, int m);
        int i4_max ( int i1, int i2 );
        int i4_min ( int i1, int i2 );
        double *multinormal_sample ( int m, int n, double a[], double mu[], int *seed );
        double r8_uniform_01 ( int *seed );
        double *r8po_fa ( int n, double a[] );
        double *r8vec_normal_01_new ( int n, int *seed );
        double *r8vec_uniform_01_new ( int n, int *seed );

#endif

