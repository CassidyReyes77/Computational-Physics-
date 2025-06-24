module user_input
  use nrtype
  implicit none
contains
!------------------------------------------------------------------------------
subroutine input1(x1_0,y1_0,vx1_0,vy1_0,m_prim,m_1,h,n_t,filename1)
  real(DP) :: x1_0, y1_0
  real(DP) :: vx1_0, vy1_0 
  real(DP) :: m_prim, m_1, h
  integer :: n_t
  character (len=100) :: filename1

  ! All commented values are for circular orbit case
  write(*,*) "Input the mass M of the primary (this is the body the object will orbit) : " 
  read(*,*)  m_prim
  ! m_prim = 1.0
  
  write(*,*) "Input the mass of the object m1 : "
  read(*,*)  m_1
  ! m_1 = 1.0
  
  write(*,*) "Input the initial x and y positions of object 1: " 
  read(*,*)  x1_0, y1_0
 !  x1_0 = 1.0
 !  y1_0 = 0.0
  
  write(*,*) "Input the intial velocity in the x and y direction of object 1:" 
  read(*,*)  vx1_0, vy1_0
 ! vx1_0 = 0.0 
 ! vy1_0 = 1.0
  
  write(*,*) "Input the time step size (dt): "
  read(*,*)  h
 ! dt = 0.001
  
  write(*,*) "Input the number of time steps N_t: " 
  read(*,*)  n_t
 ! n_t = 10000
  
  write(*,*) "Input the name of the output file to used for graphing mass 1: " 
  read(*,*) filename1
 ! filename1 = 'orbitbasic1.dat'
  
  END SUBROUTINE input1
!------------------------------------------------------------------------------
  end module user_input 
