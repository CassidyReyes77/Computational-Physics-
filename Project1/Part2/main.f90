! Programmer: Cassidy Reyes

! Purpose:
! In this program we will be using the leonard jones potential to model Van der Waals forces of Argon-Argon
! inter-atomic interactions. We'll use the first and second derivative of the leonard jones potential to find
! where the potential is at a minimum and where there may be inflection points which we'll print to the
! screen. We'll then make 3 graphs, V, V' and V" vs. r.

! Program Requirements:
! Need this program, Project1_Part2_graphs, and run.sh to compile and get graphs all at once 
  
! Subroutines in module subroutines_functions:
  ! declare_variables(v0, a, r1, r2, h)
  ! central_diff_d1(lj_potential, v0, a, r(i), h)
  ! central_diff_d2(lj_potential, v0, a, r(i), h)

! Functions in module subroutines_functions:
  ! lj_potential(v0, a, r)
 
! Variables:         Type          Units                 Description
  !  v0         | real(8)       | KJ/mol            | Lowest/initial potential for Argon-Argon interactions 
  !  i          | integer       | N/A               | an iterative integer 
  !  a          | real(8)       | angstroms         | Distance at which potential = 0
  !  r1         | real(8)       | angstroms         | Our lower boundary for r
  !  r2         | real(8)       | angstroms         | Our upper boundary for r 
  !  h          | real(8)       | N/A               | Step size or delta_r
  !  n          | integer       | N/A               | number of points between r1 ad r2 for step size h
  !  r          | real(8) array | angstroms         | An array to store all values of r between r1 and r2
  !  v          | real(8) array | KJ/mol            | An array to store values of the lj potential for all r 
  !  v1         | real(8) array | KJ/mol/angstrom   | An array to store values of the first derivative for all r
  !  v2         | real(8) array | KJ/mol/angstron^2 | An array to store values of the second derivative for all r
  ! v_min       | real(8)       | KJ/mol            | The minimum potential 
  ! r_min       | real(8)       | angstroms         | The minumum value of r that corresponds to v_min
  ! r_inflection| real(8)       | angstroms         | The value of r where v has an inflection point 


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
module subroutines_functions
  implicit none

contains
!-------------------------------------------------------------------------------
subroutine declare_variables(v0, a, r1, r2, h)
  implicit none
  real(8), intent(out) :: v0, a, h, r1, r2 

  v0 = 0.997000000000000
  a  = 3.400000000000000
  r1 = 3.200000000000000
  r2 = 8.000000000000000

  write(*,*) 'Input a step size h:'
  read(*,*) h
  
  end subroutine declare_variables 
!-------------------------------------------------------------------------------
real(8) function lj_potential(v0, a, r)
  implicit none
  real(8), intent(in) :: v0, a, r

  lj_potential = 4.0*v0*((a/r)**12.0 - (a/r)**6.0)

end function lj_potential
!-------------------------------------------------------------------------------
real(8) function central_diff_d1(f, v0, a, x, h)
  implicit none
  real(8), intent(in) ::  x, h, v0, a
  real(8), external :: f

 ! Use 2nd order central difference formula to calculate the first derivative
     central_diff_d1 = (f(v0, a, x + h) - f(v0, a, x - h))/(2.0*h)

end function central_diff_d1
!-------------------------------------------------------------------------------
real(8) function central_diff_d2(f, v0, a, x, h)
  implicit none
  real(8), intent(in) :: h, x, v0, a
  real(8), external :: f

  ! Use 4th order central difference formula to calculate the second derivative
     central_diff_d2 = (-f(v0, a, x+2.0*h) + 16.0*f(v0, a, x+h) - 30.0*f(v0, a, x) &
          + 16.0*f(v0, a, x-h)  - f(v0, a, x -2.0*h)) / (12.0*h*h)

end function central_diff_d2
!-------------------------------------------------------------------------------
end module subroutines_functions
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

program leonard_jones_potential
  use subroutines_functions 
  implicit none
  real(8) :: v0, a, r1, r2, h
  real(8) :: v_min, r_min, r_inflection
  integer :: i, n
  real(8), dimension(:), allocatable :: v, r, v1, v2

  call  declare_variables(v0, a, r1, r2, h)
  
! Calculate the number of points between r1 ad r2 for step size h
  n = int((r2 - r1) / h) + 1

  allocate(v(n), v1(n), v2(n), r(n))

! Open files 
  open(unit=10, file='LJ_original_function.dat', status='unknown')
  open(unit=20, file='LJ_derivative1.dat', status='unknown')
  open(unit=30, file='LJ_derivative2.dat', status='unknown')

 
! populate arrays for graphing v(r) vs r
  do i=1, n
     r(i) = r1 + (i-1) * h
     v(i) = lj_potential(v0, a, r(i))
     write(10,*) r(i), v(i)
  end do

! Initialize values 
  r_min  = -1.0
  r_inflection = -1.0
  
! create v1 and v2 arrays 
  do i=1, n
     v1(i) = central_diff_d1(lj_potential, v0, a, r(i), h)
     v2(i) = central_diff_d2(lj_potential, v0, a, r(i), h)
     write(20,*) r(i), v1(i)
     write(30,*) r(i), v2(i)

! Identify the minimum
    if (v1(i) > 0.0 .and. v1(i-1) < 0.0 .and. v2(i) > 0.0) then
        r_min = r(i)
    end if

! Identify inflection points by a sign change in the second derivative
    if ((v2(i) > 0.0 .and. v2(i-1) < 0.0) .or. (v2(i) < 0.0 .and. v2(i-1) > 0.0)) then
        write(*,*) 'Inflection point at v=', v(i), 'r=', r(i)
    end if     
  end do

  if (r_min /= -1.0) then
      v_min = lj_potential(v0, a, r_min)
      write(*,*) 'Minimum point at    v=', v_min, 'r=', r_min
  end if
  
  close(10)
  close(20)
  close(30)

  deallocate(v, v1, v2, r)

  end program leonard_jones_potential 
