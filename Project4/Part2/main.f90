! Purpose:
! The purpose of this program is to plot the orbit of a single object orbiting a primary given
! initial conditions of the system and a set of equations that map the position and velocity of
! such an object in time. This program solves the system of equations to map the orbit with the
! 4th order runge kutta method.
!
! Programmer: Cassidy Reyes
!---------------------------------------------------------------------------------------------------------------------

program orbit_advanced
  use nrtype
  use constants
  use runge_kutta
  implicit none
  
  real(DP)                            :: x1_0         !Initial position of m_1 in x direction 
  real(DP)                            :: y1_0         !Initial position of m_1 in y direction
  real(DP)                            :: x2_0         !Initial position of m_2 in x direction 
  real(DP)                            :: y2_0         !Initial position of m_2 in y direction 
  real(DP)                            :: vx1_0        !Initial velocity of m_1 in x direction
  real(DP)                            :: vy1_0        !Initial velocity of m_1 in y direction
  real(DP)                            :: vx2_0        !Initial velocity of m_2 in x direction
  real(DP)                            :: vy2_0        !Initial velocity of m_2 in y direction 
  real(DP)                            :: m_prim       !Mass of the object m_1 is orbiting around 
  real(DP)                            :: m_1          !Mass of object 1
  real(DP)                            :: m_2
  real(DP)                            :: h            !Step size 
  integer                             :: n_t          !Number of points or steps to evaluate 
  integer                             :: i            
  real(DP), dimension(:), allocatable :: x1           !Array that maps the position of m_1 in the x direction in time
  real(DP), dimension(:), allocatable :: y1           !Array that maps the position of m_1 in the y direction in time
  real(DP), dimension(:), allocatable :: x2           !Array that maps the position of m_2 in the x direction in time
  real(DP), dimension(:), allocatable :: y2           !Array that maps the position of m_2 in the y direction in time    
  real(DP), dimension(:), allocatable :: vx1          !Array that maps the velocity of m_1 in the x direction in time 
  real(DP), dimension(:), allocatable :: vy1          !Array that maps the velocity of m_1 in the y direction in time
  real(DP), dimension(:), allocatable :: vx2          !Array that maps the velocity of m_2 in the x direction in time 
  real(DP), dimension(:), allocatable :: vy2          !Array that maps the velocity of m_2 in the y direction in time 
  real(DP), dimension(:), allocatable :: t            !Time 
  real(DP), dimension(:), allocatable :: energy_total !Array that maps the total Energy of the system in time 
  real(DP), dimension(:), allocatable :: y0           !Array to store initial values needed to carry out runge kutta 
  real(DP), dimension(:), allocatable :: y4           !Array to store the results of the runge kutta method on the system
  character(len=100)                  :: filename2    !Name of the file that stores all of the data 


  call input_mass_and_stepsizes(m_prim,m_1,m_2,h,n_t)
  call input_initial_conditions(x1_0, y1_0, x2_0, y2_0, vx1_0, vy1_0, vx2_0, vy2_0)
  call input_file(filename2)
  
  allocate(x1(0:n_t), y1(0:n_t))
  allocate(x2(0:n_t), y2(0:n_t))
  allocate(vx1(0:n_t), vy1(0:n_t))
   allocate(vx2(0:n_t), vy2(0:n_t))
  allocate(t(0:n_t), energy_total(0:n_t))
  allocate(y0(n2), y4(n2))
    
  open(unit=10, file=trim(filename2), status='unknown', action='write')


  ! initialize arrays
  y0     = [x1_0, y1_0, x2_0, y2_0, vx1_0, vy1_0, vx2_0, vy2_0]
  x1(0)  = y0(1)
  y1(0)  = y0(2)
  x2(0)  = y0(3)
  y2(0)  = y0(4)
  vx1(0) = y0(5)
  vy1(0) = y0(6)
  vx2(0) = y0(7)
  vy2(0) = y0(8) 
  t(0)   = 0.0 
       
  part1: do i=0, n_t-1
     t(i) = real(i)*h
     
     call rk4_step(y0, h, m_prim, m_1, m_2, y4)
     x1(i)      = y4(1)
     y1(i)      = y4(2)
     x2(i)      = y4(3)
     y2(i)      = y4(4)
     vx1(i)     = y4(5)
     vy1(i)     = y4(6)
     vx2(i)     = y4(7)
     vy2(i)     = y4(8) 
     y0         = y4
     energy_total(i) = energy_eq_2(m_prim, m_1, m_2, x1(i), y1(i), x2(i), y2(i), vx1(i), vy1(i), vx2(i), vy2(i))

  !   if (m_2 == 0.0) then 
     write(10,*) t(i), energy_total(i), x1(i), y1(i), x2(i), y2(i) 
  end do part1
  
  write(*,*) 'Data written to file:', trim(filename2)

  deallocate(x1, y1, x2, y2, vx1, vy1, vx2, vy2, energy_total, t)
  deallocate(y0, y4)

  close(10)
    
  end program orbit_advanced
  
