# Compilation steps (keeping order intact)
gfortran -Wall -fcheck=bounds -c nrtype.f90 nr.f90 nrutil.f90 constants.f90
gfortran -Wall -fcheck=bounds -c pythag.f90 matrices_arrays.f90
gfortran -Wall -fcheck=bounds -c ludcmp.f90 lubksb.f90 svdcmp.f90 svbksb.f90
gfortran -Wall -fcheck=bounds -c Reyes_Cassidy_Project3.f90

# Linking steps
gfortran -Wall -fcheck=bounds -o draft1 nrtype.o nr.o nrutil.o constants.o pythag.o matrices_arrays.o ludcmp.o lubksb.o svdcmp.o svbksb.o Reyes_Cassidy_Project3.o

./draft1

# Cleanup step
rm constants.mod constants.o lubksb.mod lubksb.o ludcmp.mod ludcmp.o
rm  svdcmp.mod svdcmp.o svbksb.o pythag.o
rm matrices_arrays.mod matrices_arrays.o Reyes_Cassidy_Project3.o
rm nr.mod nr.o nrtype.mod nrtype.o nrutil.mod nrutil.o
