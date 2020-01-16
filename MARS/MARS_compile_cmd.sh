g++ -std=c++14 ./mars.cpp ./PostCal.cpp Util.cpp -o MARS_Linux.o -w -Wall -g -I ./ -I ./armadillo/include/ -DARMA_DONT_USE_WRAPPER -llapack -lblas -lgslcblas -lgsl
g++ -std=c++14 ./MARS_alt.cpp -w -o MARS_alt.o -I /boost_1_72_0/ /boost_1_72_0/stage/lib/libboost_program_options.a 
