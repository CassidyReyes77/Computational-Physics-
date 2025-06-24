module user_input
  use nrtype
  implicit none

  contains 
!-------------------------------------------------------------------------------
    subroutine input(m, v0, a, R, L, n, filename)
      real(DP)           :: m    !MeV
      real(DP)           :: v0   !MeV
      real(DP)           :: a    !fm
      real(DP)           :: R    !fm
      real(DP)           :: L    !fm
      integer            :: n
      character(len=100) :: filename

      write(*,*) "Input the mass of the particle in MeV"
     read(*,*) m
    !  m = 939.0 

      write(*,*) "Input the initial voltage in MeV"
     read(*,*) v0
    !  v0 = 50.0

      write(*,*) "Input the radius of the particle a in fm"
     read(*,*) a
    !  a = 0.2 

      write(*,*) "Input the radius, R, of the well used in the Woods-Saxon potential"
     read(*,*) R
    !  R = 1.0 

      write(*,*) "Input the size of the box L (this should be 2x R)"
     read(*,*) L
    !  L = 2.0

      write(*,*) "Input the number of grid points N "
     read(*,*) n
    !  n = 10

      write(*,*) "Input the name of the file that will be used to graph x vs. V"
      read(*,*) filename
     ! filename = "r_vs_energy.dat"

    end subroutine input
!-------------------------------------------------------------------------------
    end module user_input 

  
