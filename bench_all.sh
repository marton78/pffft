#!/bin/bash

FFTW="ON"
CMAKEOPT=""
# CMAKEOPT="-DUSE_NEON=ON"

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

echo "" >ToolChain.cmake
if [ -z "${GCC_WITH_CMAKE}" ]; then
  GCC_WITH_CMAKE="gcc"
else
  GCCPATH=$(basename "${GCC_WITH_CMAKE}")
  echo "SET(CMAKE_C_COMPILER     ${GCCPATH})" >>ToolChain.cmake
fi
if [ -z "${GPP_WITH_CMAKE}" ]; then
  GPP_WITH_CMAKE="g++"
else
  GPPPATH=$(basename "${GPP_WITH_CMAKE}")
  echo "SET(CMAKE_CXX_COMPILER   ${GPPPATH})" >>ToolChain.cmake
fi


#cmake -DCMAKE_TOOLCHAIN_FILE=ToolChain.cmake -DUSE_BENCH_FFTW=${FFTW} -DUSE_SIMD=OFF ${CMAKEOPT} ../
#make clean
#make
#echo -e "\n\nrunning without simd (==scalar) .."
#time ctest -V

cmake -DCMAKE_TOOLCHAIN_FILE=ToolChain.cmake -DUSE_BENCH_FFTW=${FFTW} -DUSE_SIMD=ON ${CMAKEOPT} ../
#make clean
make
echo -e "\n\nrunning with simd .."
time ctest -V


echo "$@" >infos.txt
echo "FFTW=${FFTW}" >>infos.txt
echo "CMAKEOPT=${CMAKEOPT}" >>infos.txt


echo "" >>infos.txt
echo "${GCC_WITH_CMAKE} --version:" >>infos.txt
${GCC_WITH_CMAKE} --version &>>infos.txt

echo "" >>infos.txt
echo "${GPP_WITH_CMAKE} --version:" >>infos.txt
${GPP_WITH_CMAKE} --version &>>infos.txt


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


tar zcvf ../pffft_bench_${GCCPATH}_${HOSTNAME}.tar.gz --exclude=CMakeCache.txt *.csv *.txt *.png
echo "all benchmark results in pffft_bench_${GCCPATH}_${HOSTNAME}.tar.gz"
