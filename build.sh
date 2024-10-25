#!/bin/bash

cpu_cores=$(nproc --all)

FoBiS.py build -fc "gfortran-8 -cpp -fopenmp -O3 -ffree-line-length-none -frecursive" -j $cpu_cores
