# Compilation steps (keeping order intact)
gfortran -Wall -fcheck=all -c nrtype.f90
gfortran -Wall -fcheck=all -c nrutil.f90
gfortran -Wall -fcheck=all -c pythag.f90 

gfortran -Wall -fcheck=all -c constants.f90
gfortran -Wall -fcheck=all -c input.f90
gfortran -Wall -fcheck=all -c tred2.f90
gfortran -Wall -fcheck=all -c tqli.f90
gfortran -Wall -fcheck=all -c sort.f90
gfortran -Wall -fcheck=all -c functions.f90
gfortran -Wall -fcheck=all -c main.f90

# Linking steps
gfortran -Wall -fcheck=bounds -o code nrtype.o nrutil.o pythag.o constants.o input.o functions.o tred2.o tqli.o sort.o main.o

./code
 python graph1.py 

# Cleanup step
rm constants.mod constants.o
rm user_input.mod input.o
rm functions.mod functions.o
rm nrtype.mod nrtype.o
rm nrutil.mod nrutil.o
rm main.o
rm code
rm tred2.o tred2_.mod
rm tqli.o  tqli_.mod
rm sort.o  sort_.mod
rm pythag.mod pythag.o

