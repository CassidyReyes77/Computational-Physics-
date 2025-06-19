#!/bin/bash

gfortran -Wall -c constants.f90 input.f90 integration.f90 Project2_Cassidy_Reyes.f90
gfortran -Wall -fcheck=bounds -o integration constants.o input.o integration.o Project2_Cassidy_Reyes.f90
./integration 
