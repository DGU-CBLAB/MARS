g++ -std=c++14 ./mars.cpp ./PostCal.cpp Util.cpp -o MARS_Linux.o -w -Wall -g -I ./ -I ./armadillo/include/ -DARMA_DONT_USE_WRAPPER -llapack -lblas -lgslcblas -lgsl
g++ -std=c++11 ./MARS_alt.cpp -o MARS_alt.o -w -I "C:\boost_1_72_0" C:\boost_1_72_0\stage\lib\libboost_program_options-vc142-mt-s-x64-1_72.lib
