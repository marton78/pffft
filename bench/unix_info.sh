#!/bin/bash

lscpu > unix_lscpu.txt
cat /proc/cpuinfo > unix_cpuinfo.txt
lsb_release -a  > unix_lsb_release.txt
cp /etc/*-release ./
