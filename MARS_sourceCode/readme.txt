###########################################################################################################
###             			 readme for MARS software			        	###
###				MARS         |          v1.0         |  01/Sep/2018			###
###				(C) 2018 Farhad Hormozdiari and Jong Wha Joanne Joo			###
###					   GNU General Public License, v3				###
###				For documentation, citation & bug-report instructions:			###
###					  http://genetics.cs.ucla.edu/MARS/				###
###########################################################################################################

###########################################################################################################
###		              http://genetics.cs.ucla.edu/MARS/install.html			        ###
###            follow instructions in readme_testData.txt for running the test dataset !!!              ###
###########################################################################################################

###########################################################################################################
0. Test environment
macOS High Sierra version 10.13.6(17G3025)
iMac (Retina 5K, 29-inch, 2017), 64-bit Processor
Processor 3.8 GHz Intel Core i5
Memory 32 GB 2400 MHz DDR4  
R version 3.5.1
GCC 4.8.5 20150623

###########################################################################################################
1. Installation (Requirement)
R(R library : Matrix, mvtnorm)
g++, gsl

2. run MARS
Goal: Compute LRT score for the given statistics and genotypes
Usage   : ./MARS -z stat [-x genotype or -l LD_matrix] -n number_of_samples -m number_of_simulations -o Output [-a 0_for_null/1_for_alternative(default:1)] [-c number_of_causal_variants(default:2)]
Option:
m: SNP, n: Sample
-z stat file path (1x50, 50 number of summary statistics) or (number_of_simulationsx101, for the null analysis)
-x genotype file path (m x n)
-l ld file path (m x m)
-n number of samples (individuals)
-o Output file path
-a if you want to analyze null statistics, set 1, default:0
-c number of causal variants to consider in the analysis
Output: (LRT_score, univariate_pvalue) or (weight LRT_score, univariate_pvalue for the null analysis)

###########################################################################################################
3. Analysis example with a test data
 
3.1. Compute LRT score for the data

3.1.1. Select top 50 stat and corresponding SNPs
Usage    : R CMD BATCH '--args -g=genotypePath -s=statPath -o=Output_genotype -u=Output_stat [-t=topNum(default:50)]' MARS_alt.R
Input    : stat(m x 1),  geno(m x n)
Output   : stat(1 x 50) geno(50 x n)  
Command  : R CMD BATCH '--args -g=sample_data/test_GENO -s=sample_data/test_STAT -o=sample_data/test_GENO50 -u=sample_data/test_STAT50 -t=50' MARS_alt.R

3.1.2. Run MARS to compute LRT score
Usage    : ./MARS -z stat [-x genotype or -l LD_matrix] -n number_of_samples -o Output [-a 0_for_null/1_for_alternative(default:1)] [-c number_of_causal_variants(default:2)]  
Input    : stat (1x50) geno(m x n) or ld(m x m)
Output   : pvalue_UNI LRT_score (2x1)
3.1.2.1. Use genotypes: #use snp to generate ld, note that the snp should make positive semidefinite ld matrix, if not, use ld matrix option -l
Command  : ./MARS -z sample_data/test_STAT50 -x sample_data/test_GENO50 -n 338 -o sample_data/test_Output -a 1
3.1.2.2. Use ld matrix:
Command  : R CMD BATCH '--args -g=sample_data/test_GENO50 -o=sample_data/test_LD50' generateLD.R
	   ./MARS -z sample_data/test_STAT50 -l sample_data/test_LD50 -n 338 -o sample_data/test_Output2 -a 1

3.2 Compute LRT score for the null

3.2.1. Generate null samples
Usage   : R CMD BATCH '--args number_of_simulations genotype Output [MARS/fastMARS(0/1, default:0)]' MARS_NULL.R
Input   : genotype(mxn) 
Output  : weight [stat1 index1 stat2 index2 ... stat50 index50](number_of_simulations x 101)
3.2.1.1. MARS
R CMD BATCH '--args -n=10000 -g=sample_data/test_GENO -o=sample_data/test_NULL -f=0 -t=50' MARS_NULL.R
3.2.1.2. fastMARS
R CMD BATCH '--args -n=10000 -g=sample_data/test_GENO -o=sample_data/test_NULL2 -f=1 -t=50' MARS_NULL.R

3.2.2. Run MARS on the null samples 
Usage   : ./MARS -z stat [-x genotype or -l LD_matrix] -n number_of_samples -m number_of_simulations -o Output [-a 0_for_null/1_for_alternative(default:1)] [-c number_of_causal_variants(default:2)]
Input   : stat (number_of_simulations x 101) geno(m x n)
Output  : weight pvalue_UNI LRT_score (number_of_simulationsx3)
3.2.2.1. MARS
Command: ./MARS -z sample_data/test_NULL -x sample_data/test_GENO -n 338 -o sample_data/test_NULL_Output -a 0 -m 10000
3.2.2.2. fastMARS
Command: ./MARS -z sample_data/test_NULL2 -x sample_data/test_GENO -n 338 -o sample_data/test_NULL2_Output -a 0 -m 10000

3.3. Compute Pvalue
Usage: R CMD BATCH '--args -a=LRT_data -n=LRT_null -o=Output [-t=threshold(default:0.05)] computePvalue.R
3.3.1. Compute pvalue to identify the multiloci association
Description: order the LRT_scores from the null samples and find the qunatile of the LRT_score from the data to compute the pvalue 
Command: R CMD BATCH '--args -a=sample_data/test_Output_LRT -n=sample_data/test_NULL_Output_LRT -o=sample_data/test_result -t=0.02797203' computePvalue.R

Usage: R CMD BATCH '--args -a=LRT_data -n=LRT_null -o=Output [-f=MARS/fastMARS(0/1, default:0)] [-u=univariate threshold(default:5e-08)]' computePvalue_GWAS.R
3.3.2. MARS for GWAS anlaysis 
Description: find the quantile of a univariate threshold(default:5e-08) from univariate pvalues from null samples. Find the LRT_threshold by finding the LRT_score of the quantile from the LRT_scores of null samples. Check if the LRT_score of the data is greater than the LRT_threshold to find the significance. Check the manuscript for the details.
Command: R CMD BATCH '--args -a=sample_data/test_Output_LRT -n=sample_data/test_NULL_Output_LRT -o=sample_data/test_result_GWAS -f=0 -u=5e-6' computePvalue_GWAS.R
3.3.3. fastMARS for GWAS anlaysis
Description: Use the weights to find the significance. Check the manuscript for the details.
Command: R CMD BATCH '--args -a=sample_data/test_Output_LRT -n=sample_data/test_NULL2_Output_LRT -o=sample_data/test_result_GWAS2 -f=1 -u=5e-6' computePvalue_GWAS.R
###########################################################################################################


###########################################################################################################
### Additional analysis scripts ###
###########################################################################################################

### generate LD file
Usage: R CMD BATCH '--args -g=genotypePath -o=OutputPath' generateLD.R
Input: genotypefile(50xn)
Output: LDmatrix(50x50)
R CMD BATCH '--args -g=sample_data/test_GENO50 -o=sample_data/test_LD50' generateLD.R
