#!/bin/sh
rm -f openmptest.o
rm -f openmptest
make openmptest
export OMP_NUM_THREADS=4
time ./openmptest 

