###########################################################################################################
###             			 readme for MARS software			        	###
###				MARS         |          v1.0         |  01/Sep/2018			###
###				(C) 2018 Farhad Hormozdiari and Jong Wha Joanne Joo			###
###					   GNU General Public License, v3				###
###				For documentation, citation & bug-report instructions:			###
###					  http://genetics.cs.ucla.edu/MARS/				###
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

###########################################################################################################
2. Running example (copy and paste each line in the MARS directory, expected run time on a normal desktop computer is described at #3 and the expected outputs are in the MARS/expectedResult/ directory)

<MARS>
R CMD BATCH '--args -g=sample_data/test_GENO -s=sample_data/test_STAT -o=sample_data/test_GENO50 -u=sample_data/test_STAT50 -t=50' MARS_alt.R
./MARS -z sample_data/test_STAT50 -x sample_data/test_GENO50 -n 338 -o sample_data/test_output -a 1
R CMD BATCH '--args -n=10000 -g=sample_data/test_GENO -o=sample_data/test_NULL -f=0 -t=50' MARS_NULL.R
./MARS -z sample_data/test_NULL -x sample_data/test_GENO -n 338 -o sample_data/test_NULL_output -a 0 -m 10000
R CMD BATCH '--args -a=sample_data/test_output_LRT -n=sample_data/test_NULL_output_LRT -o=sample_data/test_result -t=0.02797203' computePvalue.R
R CMD BATCH '--args -a=sample_data/test_output_LRT -n=sample_data/test_NULL_output_LRT -o=sample_data/test_result_GWAS -f=0 -u=5e-6' computePvalue_GWAS.R

<fastMARS>
R CMD BATCH '--args -g=sample_data/test_GENO -s=sample_data/test_STAT -o=sample_data/test_GENO50 -u=sample_data/test_STAT50 -t=50' MARS_alt.R
./MARS -z sample_data/test_STAT50 -x sample_data/test_GENO50 -n 338 -o sample_data/test_output -a 1
R CMD BATCH '--args -n=10000 -g=sample_data/test_GENO -o=sample_data/test_NULL2 -f=1 -t=50' MARS_NULL.R
./MARS -z sample_data/test_NULL2 -x sample_data/test_GENO -n 338 -o sample_data/test_NULL2_output -a 0 -m 10000
R CMD BATCH '--args -a=sample_data/test_output_LRT -n=sample_data/test_NULL2_output_LRT -o=sample_data/test_result_GWAS2 -f=1 -u=5e-6' computePvalue_GWAS.R

###########################################################################################################
3. Details of the run example

3.1. Compute LRT score for the data

3.1.1. select top 50 stat and corresponding SNPs
R CMD BATCH '--args -g=sample_data/test_GENO -s=sample_data/test_STAT -o=sample_data/test_GENO50 -u=sample_data/test_STAT50 -t=50' MARS_alt.R
(expected run time for the test data: < 1 sec)

3.1.2. run MARS to compute LRT score
./MARS -z sample_data/test_STAT50 -x sample_data/test_GENO50 -n 338 -o sample_data/test_output -a 1
(expected run time for the test data: < 1 sec)

3.2 Compute LRT score for the null
3.2.1. generate null samples
R CMD BATCH '--args -n=10000 -g=sample_data/test_GENO -o=sample_data/test_NULL2 -f=1 -t=50' MARS_NULL.R
(expected run time for the test data: < 10 min)
3.2.2. run MARS on the null samples 
./MARS -z sample_data/test_NULL2 -x sample_data/test_GENO -n 338 -o sample_data/test_NULL2_output -a 0 -m 10000
(expected run time for the test data: < 2 hours)

3.3. Compute Pvalue
3.3.1. Compute pvalue to identify the multiloci association
R CMD BATCH '--args -a=sample_data/test_output_LRT -n=sample_data/test_NULL_output_LRT -o=sample_data/test_result -t=0.02797203' computePvalue.R
3.3.2. MARS for GWAS anlaysis 
R CMD BATCH '--args -a=sample_data/test_output_LRT -n=sample_data/test_NULL_output_LRT -o=sample_data/test_result_GWAS -f=0 -u=5e-6' computePvalue_GWAS.R
3.3.3. fastMARS for GWAS anlaysis
R CMD BATCH '--args -a=sample_data/test_output_LRT -n=sample_data/test_NULL2_output_LRT -o=sample_data/test_result_GWAS2 -f=1 -u=5e-6' computePvalue_GWAS.R
(expected run time for the test data: < 1 sec)

###########################################################################################################
