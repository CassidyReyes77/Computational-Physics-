! Programmer: Cassidy Reyes

! Purpose:
  ! In this program we'll compare the exact 1st and 2nd derivatives of our
  ! given function to the 1st derivative calculated with the 2nd order central
  ! difference formula and the 2nd derivative calculated with the 4th order
  ! central difference formula. We'll do this by calclating 
  ! |f'_exact(x)-f'_approx(x)| and |f"_exact(x)-f"_approx(x)| graphing the
  ! result vs. delta_x or h(i) in two graphs.

! Graphs:
!Graphing was done in jupyter lab using the two output files made in this
! program and a python3 kernal. Code for graphing is in the file graphs2.py

  
! Subroutines in module:   
  ! input(h,x0)
  ! exact_derivative1(x0, h, exact_d1)
  ! exact_derivative2(x0, h, exact_d2)
  ! central_diff_d1(f, x0, h, i, central_d1)
  ! central_diff_d2(f, x0, h, i, central_d2)
  ! Functions in Module:
  ! f(x)
 
! Variables:          Type      Units  Description
  !  x0        | real(8)       | N/A | A user input point for calculation 
  !  i         | integer       | N/A | Integer value for iteration 
  ! central_d2 | real(8) array | N/A | Result of 4th order 2nd derivative central 
  !            |               |     | difference calculation at pt. x0 for each h
  ! exact_d1   | real(8)       | N/A | Exact 1st derivative at pt. x0 
  ! exact_d2   | real(8)       | N/A | Exact 2nd derivative at pt. x0  
  ! central_d1 | real(8) array | N/A | Result of 2nd order 1st derivative central 
  !            |               |     | difference calculation at pt. x0 for each h
  ! graph_d1   | real(8) array | N/A | Absolute difference between each 1st
  !            |               |     | derivative for each h value               
  ! graph_d2   | real(8) array | N/A | Absolute difference between each 2nd 
  !            |               |     | derivative for each h value    
  !    h       | real(8) array | N/A | User input step size (.1, .01, .001, .0001)

! Note:
  ! h(i) is treated as an array instead of a single value to simplify the steps
  !   involved in graphing and comparing h values. Instead of asking for a single h
  ! and compiling the program 4 times value, we ask the user for all h(i) and
  ! find all needed values in one compilation. In this program we also ask the user
  ! to input point x0. 

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------

Module derivative_subroutines_and_functions
    implicit none
    
  contains
    
!------------------------------------------------------------------------------
subroutine input(h,x0)
  implicit none
  real(8), intent(out) :: x0
  real(8), dimension(:), allocatable, intent(out) :: h

  allocate(h(4))

  ! Ask for user input 
  write(*,*) "Input 4 values for step sizes h1, h2, h3, h4: "
  read(*,*) h(1), h(2), h(3), h(4)
  write(*,*) "Input a point x0 for calculation: "
  read(*,*) x0

  write(*,*) "You entered: "
  write(*,*) h(1), h(2), h(3), h(4), x0

end subroutine input
!------------------------------------------------------------------------------
real(8) function f(x)
  implicit none
  real(8), intent(in) :: x

  f = x * sinh(1.0/x)

end function f
!------------------------------------------------------------------------------
subroutine exact_derivative1(x0,exact_d1)
  implicit none
  real(8), intent(in) :: x0
  real(8),dimension(:), allocatable :: exact_d1
  integer :: i 

 ! Calculate the exact first derivative at pt. x0. All values in the array
 ! should be equal since h isnt in the equation.
  do i=1, 4
     exact_d1(i) = sinh(1.0/x0) - cosh(1.0/x0)/x0
  end do

 ! Print results to screen. 
  write(*,*) "The exact first derivative at point x0 is: "
  write(*,*) exact_d1(1)
  
end subroutine exact_derivative1
!-------------------------------------------------------------------------------
subroutine exact_derivative2(x0, exact_d2)
  implicit none
  real(8), intent(in) :: x0
  real(8), dimension(:), allocatable :: exact_d2
  integer :: i

 ! Calculate the exact second derivative at pt. x0. All values in the array
 ! should be equal since h isnt in the equation. 
  do i=1, 4
     exact_d2(i) = sinh(1.0/x0) / x0**3.0
  end do
  
 ! Print result to screen. 
   write(*,*) "The exact second derivative at point x0 is: "
   write(*,*) exact_d2(1)
   
end subroutine exact_derivative2
!-------------------------------------------------------------------------------
subroutine central_diff_d1(f, x0, h, i, central_d1)
  implicit none
  real(8), dimension(:), allocatable, intent(in) :: h
  real(8), intent(in) ::  x0
  real(8), external :: f
  real(8), dimension(:), allocatable :: central_d1
  integer :: i

 ! Use 2nd order central difference formula to calculate the first derivative
 ! and print the result to the screen for each value of h
  
  write(*,*) "The 1st deriv. at pt. x0 using 2nd order central diff. is:"
  do i=1, 4
     central_d1(i) = (f(x0 + h(i)) - f(x0 - h(i)))/(2.0*h(i))
     write(*,*) 'h=', h(i),':', central_d1(i)
     end do 

end subroutine central_diff_d1
!-----------------------------------------------------------------------------------
subroutine central_diff_d2(f, x0, h, i, central_d2)
  implicit none
  real(8), dimension(:), allocatable, intent(in) :: h
  real(8), intent(in) ::  x0
  real(8), external :: f
  real(8), dimension(:), allocatable :: central_d2
  integer :: i

  ! Use 4th order central difference formula to calculate the second derivative
  ! and print the result to the screen for each value of h
  
  write(*,*) "The 2nd deriv. at pt. x0 using 4th order central diff. is:"
  do i=1, 4
     central_d2(i) = (-f(x0+2.0*h(i)) + 16.0*f(x0+h(i)) - 30.0*f(x0) &
          + 16.0*f(x0-h(i))  - f(x0 -2.0*h(i))) / (12.0*h(i)*h(i))
     write(*,*) 'h=', h(i),':', central_d2(i)
     end do 

end subroutine central_diff_d2
!-----------------------------------------------------------------------------------

end module derivative_subroutines_and_functions

!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------

program main_project1_part1
  use derivative_subroutines_and_functions

  implicit none
  real(8) :: x0
  integer :: i
  real(8), dimension(:), allocatable :: h, central_d2, central_d1
  real(8), dimension(:), allocatable ::graph_d1, graph_d2, exact_d1, exact_d2

  ! Call input function
  call input(h,x0)

  allocate(graph_d1(4),graph_d2(4), central_d1(4),central_d2(4))
  allocate(exact_d1(4), exact_d2(4))
  
 ! Open files for graphing 
  OPEN(UNIT=10,FILE='derivative1.dat',STATUS='UNKNOWN')
  OPEN(UNIT=20,FILE='derivative2.dat',STATUS='UNKNOWN')

 ! Call subroutines and calculate the functions exact and central difference derivs.
  call exact_derivative1(x0, exact_d1)
  call exact_derivative2(x0, exact_d2)
  call central_diff_d1(f, x0, h, i, central_d1)
  call central_diff_d2(f, x0, h, i, central_d2)


 ! Graph functions |f'_exact(x0)-f'_approx(x0)| and |f"_exact(x0)-f"_approx(x0)|
  do i=1, 4
     graph_d1(i) = abs(exact_d1(i) - central_d1(i))
     graph_d2(i) = abs(exact_d2(i) - central_d2(i))
     write(10,*) h(i), graph_d1(i)
     Write(20,*) h(i), graph_d2(i)
  end do

 ! Close files 
  close(10)
  close(20)

 ! Deallocate arrays 
  deallocate(central_d1, central_d2, graph_d1, graph_d2, h)
  deallocate( exact_d1, exact_d2)


  end program main_project1_part1
