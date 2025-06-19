#!/bin/bash
gfortran -Wall -o lj_potential Project1_part2_Cassidy_Reyes.f90
./lj_potential
python Project1_Part2_graphs.py 
