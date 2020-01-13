#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include <unistd.h> 
#include "Util.h"
#include "PostCal.h"
#include "CaviarModel.h"


using namespace std;

int main( int argc, char *argv[]  ){
	int totalCausalSNP = 2;
	double NCP = 5.7;
	double gamma = 0.01;
	double rho = 0.95;
	int oc = 0;	
	string ldFile = "";
	string zFile  = "";
	string outputFileName = "";
	string geneMapFile = "";	
	string snpFile = "";
	int sampleSize = 0;
	int optionNum = 0;
	int simNum = 0;
	bool ldoption = false;
	bool alt = true;

	while ((oc = getopt(argc, argv, "vhl:x:o:z:n:m:g:r:c:a:")) != -1) {
		switch (oc) {
			case 'v':
				cout << "version 1.0:" << endl;
			case 'h':
				cout << "Usage: ./MARS -z statFile [-x snpFile or -l LDFile] -n sampleSize -o outputFile" << endl;
				cout << "Options: " << endl;
  				cout << "-h, --help            		show this help message and exit " << endl;
  				cout << "-z ZFILE                       summary statistics file" << endl;
				cout << "-l LDFILE, --ld_file=LDFILE  	the ld input file" << endl;
  				cout << "-x snpFILE			genotype file" << endl;	
				cout << "-n sampleSize			number of SNPs" <<endl;
				cout << "-m number of simulations for the null case" << endl;
                                cout << "-o OUTFILE, --out=OUTFILE      specify the output file" << endl;
				cout << "-c causal                      set the maximum number of causal SNPs (default 2)" << endl;
				cout << "-a null/alt" << endl;
  				cout << "-r RHO, --rho-prob=RHO		set $pho$ probability (default 0.95)" << endl;
				cout << "-g GAMMA, --gamma		set $gamma$ the prior of a SNP being causal (default 0.01)" << endl;
				exit(0);
			case 'l':
				ldFile = string(optarg);
				optionNum++;
				ldoption = true;
				break;
			case 'x':
				snpFile = string(optarg);
				optionNum++;
				break;
			case 'n':
                                sampleSize = atoi(optarg);
                                optionNum++;
				break;
			case 'm':
				simNum = atoi(optarg);
				break;
			case 'o':
				outputFileName = string(optarg);
				optionNum++;
				break;
			case 'z':
				zFile = string(optarg);
				optionNum++;
				break;
			case 'r':
				rho = atof(optarg);
				break;
			case 'c':
				totalCausalSNP = atoi(optarg);
				break;
			case 'a':
				if(atoi(optarg)==1)
					alt = true;
				else alt = false;
				break;
			case 'g':
				gamma = atof(optarg);
				break;
			case ':':
			case '?':
			default:
				cout << "Strange" << endl;
				break;
		}
	}
	if(optionNum<4){
		cout << "Usage: ./MARS -z statFile [-x snpFile or -l LDFile] -n sampleSize -o outputFile" << endl;
		exit(0);
	}
	
	cout << "@-------------------------------------------------------------@" << endl;
	cout << "| MARS		| 	   v1.0         |  01/Sep/2018         |" << endl;
	cout << "| (C) 2018 Farhad Hormozdiari and Jong Wha Joanne Joo         |" << endl;
	cout << "|              GNU General Public License, v3                 |" << endl;
	cout << "|-------------------------------------------------------------|" << endl;
	cout << "| For documentation, citation & bug-report instructions:      |" << endl;
	cout << "| 		http://genetics.cs.ucla.edu/MARS/              |" << endl;
	cout << "@-------------------------------------------------------------@" << endl;	
	std::clock_t start = std::clock();
	CaviarModel caviar(ldoption, snpFile, ldFile, zFile, outputFileName, sampleSize, simNum, totalCausalSNP, alt, NCP, rho, gamma);
	caviar.run();
	caviar.finishUp(string(outputFileName)+"_LRT");
        double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC; 	
        cout <<"running time = "<<duration<<" secs\n";				
	return 0;
}
