! Purpose:
! The purpose of this program is to plot the orbit of a single object orbiting a primary given
! initial conditions of the system and a set of equations that map the position and velocity of
! such an object in time. This program solves the system of equations to map the orbit with the
! 4th order runge kutta method.
!
! Programmer: Cassidy Reyes
!---------------------------------------------------------------------------------------------------------------------

program orbitbasic
  use nrtype
  use constants
  use runge_kutta
  implicit none
  
  real(DP)                            :: x1_0         !Initial position of m_1 in x direction 
  real(DP)                            :: y1_0         !Initial position of m_1 in y direction 
  real(DP)                            :: vx1_0        !Initial velocity of m_1 in x direction
  real(DP)                            :: vy1_0        !Initial velocity of m_1 in y direction 
  real(DP)                            :: m_prim       !Mass of the object m_1 is orbiting around 
  real(DP)                            :: m_1          !Mass of the object 
  real(DP)                            :: h            !Step size 
  integer                             :: n_t          !Number of points or steps to evaluate 
  integer                             :: i            
  real(DP), dimension(:), allocatable :: x1           !Array that maps the position of m_1 in the x direction in time
  real(DP), dimension(:), allocatable :: y1           !Array that maps the position of m_1 in the y direction in time    
  real(DP), dimension(:), allocatable :: vx1          !Array that maps the velocity of m_1 in the x direction in time 
  real(DP), dimension(:), allocatable :: vy1          !Array that maps the velocity of m_1 in the x direction in time 
  real(DP), dimension(:), allocatable :: t            !Time 
  real(DP), dimension(:), allocatable :: energy1      !Array that maps the Energy of m_1 as it travels in time 
  real(DP), dimension(:), allocatable :: y0           !Array to store initial values needed to carry out runge kutta 
  real(DP), dimension(:), allocatable :: y4           !Array to store the results of the runge kutta method on the system
  character(len=100)                  :: filename1    !Name of the file that stores the data 

  call input1(x1_0,y1_0,vx1_0,vy1_0,m_prim,m_1,h,n_t,filename1)
  
  allocate(t(0:n_t), x1(0:n_t), y1(0:n_t))
  allocate(vx1(0:n_t), vy1(0:n_t), energy1(0:n_t))
  allocate(y0(n1), y4(n1))
    
  open(unit=10, file=trim(filename1), status='unknown', action='write')

  ! initialize arrays
  y0 = [x1_0, y1_0, vx1_0, vy1_0]
  x1(0)  = y0(1)
  y1(0)  = y0(2)
  vx1(0) = y0(3)
  vy1(0) = y0(4)
  t(0)   = 0.0
       
  part1: do i=0, n_t-1
     t(i) = real(i)*h
     
     call rk4_step(y0, h, m_prim, y4)
     x1(i)      = y4(1)
     y1(i)      = y4(2)
     vx1(i)     = y4(3)
     vy1(i)     = y4(4)
     y0         = y4
     energy1(i) = energy_eq_1(m_prim,m_1,x1(i),y1(i),vx1(i),vy1(i))
     
     write(10,*) t(i), energy1(i), x1(i), y1(i)
  end do part1
  
  write(*,*) 'Data for m_1 written to file:', trim(filename1)

  deallocate(x1, y1, vx1, vy1, energy1, t)
  deallocate(y0, y4)

  close(10)
    
  end program orbitbasic
  
