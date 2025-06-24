# Compilation steps (keeping order intact)
gfortran  -Wall -fcheck=all -c nrtype.f90
gfortran  -Wall -fcheck=all -c constants.f90
gfortran  -Wall -fcheck=all -c input_part2.f90
gfortran  -Wall -fcheck=all -c rk4_part2.f90
gfortran  -Wall -fcheck=all -c main_part2.f90

# Linking steps
gfortran -Wall -fcheck=bounds -o Part2 nrtype.o constants.o input_part2.o rk4_part2.o main_part2.o

./Part2
python graphs_part2.py 

# Cleanup 
rm constants.mod
rm constants.o
rm rk4_part2.o
rm runge_kutta.mod
rm user_input.mod
rm input_part2.o
rm nrtype.mod
rm nrtype.o
rm main_part2.o
rm Part2
