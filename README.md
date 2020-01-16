# MARS
currently working on full refactoring of MARS using only C++

#

## How to Compile

Using Ubuntu 18.04.3 LTS
```
sudo apt-get update
sudo apt-get install build-essential
sudo apt-get install libgsl-dev libblas-dev liblapack-dev
```
Download and install <a href="https://www.boost.org/">Boost</a>

compile boost-static library
```
bootstrap.sh
./b2 link=static
```

1. Compile MARS.cpp
```
g++ -std=c++14 ./mars.cpp ./PostCal.cpp Util.cpp -o MARS_Linux.o -w -Wall -g -I ./ -I ./armadillo/include/ -DARMA_DONT_USE_WRAPPER -llapack -lblas -lgslcblas -lgsl
```

## Test MARS

MARS_alt.o
```
./MARS_alt.o --g sample_data/test_GENO --s sample_data/test_STAT --o sample_data/test_GENO50 --u sample_data/test_STAT50 --t 50
```