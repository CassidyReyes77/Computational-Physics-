# Compilation steps (keeping order intact)
gfortran -g -Wall -fcheck=all -c nrtype.f90
gfortran -g -Wall -fcheck=all -c constants.f90
gfortran -g -Wall -fcheck=all -c input_part1.f90
gfortran -g -Wall -fcheck=all -c rk4_part1.f90
gfortran -g -Wall -fcheck=all -c main_part1.f90

# Linking steps
gfortran -Wall -fcheck=bounds -o Part1 nrtype.o constants.o input_part1.o rk4_part1.o main_part1.o

./Part1
python graphs_part1.py 

# Cleanup 
rm constants.mod
rm constants.o
rm rk4_part1.o
rm runge_kutta.mod
rm user_input.mod
rm input_part1.o
rm nrtype.mod
rm nrtype.o
rm main_part1.o
rm Part1 
