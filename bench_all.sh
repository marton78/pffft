#!/bin/bash

FFTW=ON
if [ ! -z "$1" ]; then
  FFTW="$1"
fi

if [ ! -d build ]; then
  mkdir build
  cd build
else
  cd build
  make clean
  rm *.csv *.txt *.png
fi


make clean
cmake ../ -DUSE_FFTW=${FFTW} -DUSE_SIMD=OFF
make
echo -e "\n\nrunning 2nd pass without simd (==scalar) .."
time ctest -V

make clean
cmake ../ -DUSE_FFTW=${FFTW} -DUSE_SIMD=ON
make
echo -e "\n\nrunning 1st pass with simd .."
time ctest -V


echo "$@" >infos.txt

echo "" >>infos.txt
echo "gcc --version:" >>infos.txt
gcc --version &>>infos.txt

echo "" >>infos.txt
echo "lscpu:" >>infos.txt
lscpu >>infos.txt

echo "" >>infos.txt
echo "lsb_release -a" >>infos.txt
lsb_release -a &>>infos.txt

echo "" >>infos.txt
echo "cat /etc/*-release" >>infos.txt
cat /etc/*-release &>>infos.txt


echo "" >>infos.txt
echo "cat /proc/cpuinfo:" >>infos.txt
cat /proc/cpuinfo >>infos.txt


tar zcvf ../pffft_bench_${HOSTNAME}.tar.gz *.csv *.txt *.png
echo "all benchmark results in pffft_bench_${HOSTNAME}.tar.gz"
