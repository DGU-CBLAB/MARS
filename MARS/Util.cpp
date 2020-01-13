#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream> 
#include <sstream>
#include <string>
#include <vector>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Util.h"	
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
using namespace std;

long int fact(int n) {
        if(n==0)
                return 1;
        return n* fact(n-1);
}

void copyConfigure(double *dest, double *src, int size) {
	for(int i = 0; i < size; i++) 
		dest[i] = src[i];
}

double min(double a, double b) {
	if(a>b)
		return b;
	else
		return a;
}

long int nCr(int n, int r) {
        long int result = 1;
        for(int i = n; i > n-r; i--)
                result *= i;
        return result/fact(r);
}

void printVector(char * data, int size) {
        for(int i = 0; i < size; i++)
                printf("%c, ", data[i]);
}

void printVector(int * data, int size) {
        for(int i = 0; i < size; i++)
                printf("%d, ", (int)data[i]);
}

void printVector(double * data, int size) {
        for(int i = 0; i < size; i++)
                printf("%lf, ", data[i]);
}

void diffVector(double * data1, double * data2, int size, double * result) {
	for(int i = 0; i < size; i++ ) 
		result[i] = data1[i] - data2[i];
}

void sumVector(double * data1, double * data2, int size, double * result) {
        for(int i = 0; i < size; i++ ) 
                result[i] = data1[i] + data2[i];
}

double multVector(double * data1, double * data2, int size) {
	double res = 0;
	for(int i = 0; i < size; i++ ) 
                res += data1[i] * data2[i];
	return res;
}

void dotVector(double * data1, double * data2, int size, double * result) {
	for(int i = 0; i < size; i++ ) 
                result[i] = data1[i] * data2[i];
}

void multVectorMatrix(double *vector, double * matrix, int size, double * result) {
	double total_row = 0;
	for(int i = 0; i < size; i++) {
		total_row = 0;
		for(int j = 0; j < size; j++) {
			total_row += vector[j] * matrix[i + j * size];
		}
		result[i]= total_row;
	}
}

void importData(string fileName, double * vector) {
	int index = 0;
	double data = 0;
        ifstream fin(fileName.c_str(), std::ifstream::in);
        fin >> data;
        while (fin.good()) {
                vector[index] = data;
                index++;
                fin >> data;
        }
        fin.close();
}

void importData(string fileName, int * vector) {
        int index = 0;
        double data = 0;
        ifstream fin(fileName.c_str(), std::ifstream::in);
        while( fin.good()  ){
                fin >> data;
                vector[index] = (int)data;
                index++;
        }
        fin.close();
}

/*
	The column index starts by 1 in this implemenation
*/
void importStat(string fileName, double * vector) {
        int index = 0;
        string line = "";
        double data = 0.0;
        ifstream fin(fileName.c_str(), std::ifstream::in);
        while( getline(fin, line) ){
                istringstream iss(line);
                iss >> data;
                vector[index] = (double)data;
                index++;
        }
        cout << "reach=" << index << endl;
        fin.close();
}
void importStat2(string fileName, double * vector, int subSize) {
	cout<<"subSize="<<subSize<<endl;
        int index = 0;
        string line = "";
        double data = 0.0;
        ifstream fin(fileName.c_str(), std::ifstream::in);
        getline(fin, line);
	istringstream iss(line);
	for(int i=0;i<subSize; i++){
                iss >> data;
                vector[i] = (double)data;
		index++;
        }
        cout << "reach=" << index << endl;
        fin.close();
}
void importDataSecondColumn(string fileName, double * vector) {
	int index = 0;
	string line = "";
	string dataS = "";
	double data = 0.0;	
	ifstream fin(fileName.c_str(), std::ifstream::in);
	while( getline(fin, line) ){
		istringstream iss(line);
		iss >> dataS;
		iss >> data;
	        vector[index] = (double)data;
                index++;
        }
	cout << "reach=" << index << endl;
        fin.close();
}

/*
	The column index starts by 1 in this implemenation
*/
void importDataNthColumn(string fileName, double * vector, int colNum, int ignore) {
        int index = 0;
        string line = "";
        string dataS = "";
        double data = 0.0;
	ifstream fin(fileName.c_str(), std::ifstream::in);
	for(int i = 0; i < ignore; i++)
        	getline(fin, line);

        while( getline(fin, line) ){
                istringstream iss(line);
                iss >> dataS;
                for(int i = 0; i < colNum-1;i++)
			iss >> data;
                vector[index] = (double)data;
                index++;
        }
        cout << "reach=" << index << endl;
        fin.close();
}

void importDataFirstColumn(string fileName, string * list, int ignore) {
 	int index = 0;
        string data = "";
        string line = "";
	ifstream fin(fileName.c_str(), std::ifstream::in);
	for(int i = 0; i < ignore; i++)
                getline(fin, line);
	
        while( getline(fin, line) ){
		istringstream iss(line);
                iss >> data;
                list[index] = data;
		index++;
        }
	//cout << "FINISH" << endl;
        fin.close();
}
void loadSNP(string fileName, double* snp, int snpCount, int sampleSize){
        string data_s;
        double data_d;
        int size = 0;
        ifstream fin(fileName.c_str(), std::ifstream::in);
        for(int i=0; i<snpCount; i++){
                fin >>data_s;//snpID
                for(int j=0; j<sampleSize; j++){
                        fin>>data_d;
                        snp[i*sampleSize+j] = data_d;
                }
        }
        fin.close();
}
double correlationCoefficient(double* X, double* Y, int n){

    double sum_X = 0, sum_Y = 0, sum_XY = 0;
    double squareSum_X = 0, squareSum_Y = 0;
    for (int i = 0; i < n; i++){
        sum_X = sum_X + X[i];
        sum_Y = sum_Y + Y[i];
        sum_XY = sum_XY + X[i] * Y[i];
        squareSum_X = squareSum_X + X[i] * X[i];
        squareSum_Y = squareSum_Y + Y[i] * Y[i];
    }
    double corr = (double)(n * sum_XY - sum_X * sum_Y) / sqrt((n * squareSum_X - sum_X * sum_X) * (n * squareSum_Y - sum_Y * sum_Y));
    return corr;
}
void makeSubSigma(double* snp, double* sigma, int subSize, int sampleSize){
    double* X = new double[sampleSize];
    double* Y= new double[sampleSize];
    for(int i=0;i<subSize; i++){
        for(int j=0; j<subSize; j++){
                for(int k=0; k<sampleSize; k++){
                        X[k]=snp[i*sampleSize+k];
                        Y[k]=snp[j*sampleSize+k];
                }
                sigma[i*subSize+j]=correlationCoefficient(X,Y,sampleSize);
        }
    }
}
void fileSize(string fileName, int & size) {
	size = 0;
	double data = 0;	
	ifstream fin(fileName.c_str(), std::ifstream::in);
	while( fin.good()  ){
		fin >> data;
		size++;
	}
	fin.close();
}

void fileSize_snp(string fileName, int & size) {
        size = 0;
        string data;
//      cout<<"fileSize check="<<fileName<<endl;
        ifstream fin(fileName.c_str(), std::ifstream::in);
        while( fin.good()  ){
                fin >> data;//cout<<"data="<<data<<" ";
                size++;
        }
        fin.close();
}

string convertInt(int number) {
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

void resetVector(char *data, int size){
	for(int i = 0; i < size; i++)
		data[i] = '0';
}

void resetVector(int * data, int size) {
	for(int i = 0; i < size; i++)
		data[i] = 0;
}

void resetVector(double * data, int size) {
        for(int i = 0; i < size; i++)
                data[i] = 0;
}

void exportVector2File(string fileName, char * data, int size) {
	ofstream outfile(fileName.c_str(), ios::out | ios::app);
	for (int i = 0; i < size; i++)
		outfile << data[i] << " ";
	//outfile << endl;
	outfile.close();
}

void exportVector2File(string fileName, double * data, int size) {
        ofstream outfile(fileName.c_str(), ios::out | ios::app);
        for (int i = 0; i < size; i++)
                outfile << data[i] << " ";
        //outfile << endl;
        outfile.close();
}

void exportVector2File(string fileName, int * data, int size) {
        ofstream outfile(fileName.c_str(), ios::out | ios::app);
        for (int i = 0; i < size; i++)
                outfile << data[i] << " ";
        //outfile << endl;
        outfile.close();
}

void export2File(string fileName, double data) {
        ofstream outfile(fileName.c_str(), ios::out | ios::app);
        outfile << data << endl;
        outfile.close();
}

void export2File(string fileName, int data) {
        ofstream outfile(fileName.c_str(), ios::out | ios::app);
        outfile << data << endl;
        outfile.close();
}

int snp2Gene(int * G, int snpId, int snpCount, int geneCount) {
	for(int i = 0; i < geneCount; i++) {
		if(G[snpId*geneCount + i] == 1)
			return i;
	}
	return -1;
}

void setIdentitymatrix(int * G, int snpCount, int geneCount) {
	for(int i = 0; i < snpCount; i++) {
		for(int j = 0; j < geneCount; j++) {
			G[i*geneCount + j] = 0;
		}
		G[i*geneCount + (i/(snpCount/geneCount))] = 1;
	}
	  for(int i = 0; i < snpCount; i++) {
                for(int j = 0; j < geneCount; j++) {
                        printf("%d ", G[i*geneCount+j]);
                }
		printf("\n");
        }
}

void makeSigmaPositiveSemiDefinite(double * sigma, int size) {
	int gsl_tmp = 0;
	double matDet  = 0;
        double addDiag = 0;
        bool positive = false;
	
	//gsl_set_error_handler_off();	
	gsl_matrix * tmpResultMatrix = gsl_matrix_calloc (size, size);	
	gsl_permutation *p = gsl_permutation_alloc(size);
        do{
		for(int i = 0; i < size; i++) {
                	for (int j = 0; j < size; j++) {
                      		if(i==j)
					gsl_matrix_set(tmpResultMatrix,i,j,sigma[i*size+j]+addDiag);
				else
					gsl_matrix_set(tmpResultMatrix,i,j,sigma[i*size+j]);
			}
        	}
	
		gsl_linalg_LU_decomp(tmpResultMatrix, p, &gsl_tmp);
       		matDet = gsl_linalg_LU_det(tmpResultMatrix,gsl_tmp);	
		//cout << matDet << "\t" << addDiag << endl;
		if(matDet > 0 ) 
			positive = true;
		else {
//			cout << "add ";
			addDiag+=0.1;		
		}
	} while(!positive);
	for(int i = 0; i < size*size; i++){
                if(i%(size+1) == 0)
                        sigma[i] = sigma[i] + addDiag;
        }
}


void rmvnormC(double * mean, double * sigma, int seed, int size, double* results) {
        gsl_matrix * work = gsl_matrix_alloc(size,size);
        gsl_vector * resultVector = gsl_vector_alloc(size);
        gsl_vector * meanVector  = gsl_vector_alloc(size);
        gsl_matrix * sigmaMatrix = gsl_matrix_alloc(size,size);
        const gsl_rng_type * T;
        gsl_rng * r;
        gsl_rng_env_setup();
        gsl_rng_default_seed = seed;
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);

        for(int i = 0; i < size; i++) {
                for(int j = 0; j < size; j++) {
                        gsl_matrix_set(sigmaMatrix,i,j,sigma[i*size + j]);
                }
                gsl_vector_set(meanVector,i ,mean[i]);
        }

        gsl_matrix_memcpy(work, sigmaMatrix);
        gsl_linalg_cholesky_decomp(work);
        for(int i=0; i<size; i++)
                gsl_vector_set( resultVector, i, gsl_ran_ugaussian(r) );

        gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, resultVector);
        gsl_vector_add(resultVector,meanVector);
	//cout<<"results="; for (int ii = 0; ii < size; ii++ ) cout<<results[ii]<<" ";cout<<"\n"<<endl;
        for(int i = 0; i < size; i++)
                results[i] = gsl_vector_get(resultVector, i);
	//cout<<"results="; for (int ii = 0; ii < size; ii++ ) cout<<results[ii]<<" ";cout<<"\n"<<endl;
        gsl_matrix_free(work);
        gsl_vector_free(resultVector);
        gsl_vector_free(meanVector);
        gsl_matrix_free(sigmaMatrix);
        gsl_rng_free(r);
}
double* generateMu(int snpCount, double NCP){//implant 2 causal SNPs
        int Csnp = 2;
        double* x = new double[snpCount];
        int* ranc = new int[Csnp];
        for (int i = 0; i < snpCount; i++ ) x[i] = 0;
        ranc[0] = rand() % snpCount;
        ranc[1] = rand() % snpCount;
        while (ranc[0]==ranc[1]){
                ranc[1] =rand() % snpCount;
        }
        for(int i=0; i<Csnp; i++) x[ranc[i]]=NCP;
        return x;
}
double* generateStat(double * sigma, int m, int seed, int* quantile, double threshold, double mu[]){
	double *x;
	x = multinormal_sample ( m, 1, sigma, mu, &seed );
	return x;
}
double *multinormal_sample ( int m, int n, double a[], double mu[], int *seed ){
	int i;
  	int j;
  	int k;
  	double *r;
 	double *x;
  	double *y;
  	r = r8po_fa ( m, a );

  	if ( !r ){
    		cout << "\n";
    		cout << "MULTINORMAL_SAMPLE - Fatal error!\n";
    		cout << "  The variance-covariance matrix is not positive definite symmetric.\n";
    		exit ( 1 );
  	}
  	y = r8vec_normal_01_new ( m*n, seed );
  	x = new double[m*n];

  	for( j = 0; j < n; j++ ){
    		for ( i = 0; i < m; i++ ){
      			x[i+j*m] = mu[i];
      			for ( k = 0; k < m; k++ ){
        			x[i+j*m] = x[i+j*m] + r[k+i*m] * y[k+j*m];
      			}
    		}
  	}

  	delete [] r;
  	delete [] y;
  	return x;
}
int i4_max ( int i1, int i2 )
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
int i4_min ( int i1, int i2 )
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
double r8_uniform_01 ( int *seed )
{
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
double *r8po_fa ( int n, double a[] )
{
  double *b;
  int i;
  int j;
  int k;
  double s;

  b = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i+j*n] = a[i+j*n];
    }
  }

  for ( j = 0; j < n; j++ )
  {
    for ( k = 0; k <= j-1; k++ )
    {
      for ( i = 0; i <= k-1; i++ )
      {
        b[k+j*n] = b[k+j*n] - b[i+k*n] * b[i+j*n];
      }
      b[k+j*n] = b[k+j*n] / b[k+k*n];
    }

    s = b[j+j*n];
    for ( i = 0; i <= j-1; i++ )
    {
      s = s - b[i+j*n] * b[i+j*n];
    }

    if ( s <= 0.0 )
    {
      delete [] b;
      return NULL;
    }

    b[j+j*n] = sqrt ( s );
  } 
	for ( i = 0; i < n; i++ ){
    		for ( j = 0; j < i; j++ ){
      			b[i+j*n] = 0.0;
    		}
  	}

  	return b;
}
double *r8vec_normal_01_new ( int n, int *seed )
{
# define PI 3.141592653589793

  int i;
  int m;
  static int made = 0;
  double *r;
  static int saved = 0;
  double *x;
  int x_hi;
  int x_lo;
  static double y = 0.0;

  x = new double[n];
  if ( n < 0 )
  {
    made = 0;
    saved = 0;
    y = 0.0;
    return NULL;
  }
  else if ( n == 0 )
  {
    return NULL;
  }
  x_lo = 1;
  x_hi = n;
  if ( saved == 1 )
  {
    x[0] = y;
    saved = 0;
    x_lo = 2;
  }
  if ( x_hi - x_lo + 1 == 0 )
  {
  }
else if ( x_hi - x_lo + 1 == 1 )
  {
    r = r8vec_uniform_01_new ( 2, seed );

    x[x_hi-1] = sqrt ( - 2.0 * log ( r[0] ) ) * cos ( 2.0 * PI * r[1] );
    y =         sqrt ( - 2.0 * log ( r[0] ) ) * sin ( 2.0 * PI * r[1] );
    saved = 1;

    made = made + 2;

    delete [] r;
  }
 else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
  {
    m = ( x_hi - x_lo + 1 ) / 2;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * PI * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * PI * r[i+1] );
    }
    made = made + x_hi - x_lo + 1;

    delete [] r;
  }
  else
  {
    x_hi = x_hi - 1;

    m = ( x_hi - x_lo + 1 ) / 2 + 1;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * PI * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * PI * r[i+1] );
    }
i = 2*m - 2;

    x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * PI * r[i+1] );
    y           = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * PI * r[i+1] );

    saved = 1;

    made = made + x_hi - x_lo + 2;

    delete [] r;
  }

  return x;
# undef PI
}
double *r8vec_uniform_01_new ( int n, int *seed )
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}

