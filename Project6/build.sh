#!/bin/bash

# Compiling each module
gfortran -Wall -fcheck=all -c nrtype.f90
gfortran -Wall -fcheck=all -c constants.f90
gfortran -Wall -fcheck=all -c input.f90
gfortran -Wall -fcheck=all -c functions.f90
gfortran -Wall -fcheck=all -c main.f90

# Linking to create the final executable
gfortran -Wall -o Project6 nrtype.o constants.o input.o functions.o main.o 

# Running the program
./Project6

# Running the Python script for graphing
 python graphs.py

rm constants.mod
rm constants.o
rm input.mod
rm input.o
rm main.o
rm nrtype.mod
rm nrtype.o
rm functions.o
rm functions.mod 
rm Project6


