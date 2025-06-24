module user_input
  use nrtype
  implicit none
contains
!------------------------------------------------------------------------------
subroutine input_mass_and_stepsizes(m_prim, m_1, m_2, h, n_t)
  real(DP) :: m_prim, m_1, m_2, h
  integer :: n_t

  ! All commented values are for circular orbit case
  write(*,*) "Input the mass M of the primary (this is the body the object will orbit) : "
  read(*,*)  m_prim
  !m_prim = 10.0
  
  write(*,*) "Input the mass of the object m1 : "
  read(*,*)  m_1
  ! m_1 = 1.0

  write(*,*) "Input the mass of the object m2 : "
  read(*,*)  m_2
  !m_2 = 0.0
  
  write(*,*) "Input the time step size (dt): "
  read(*,*)  h
  !dt = 0.001
  
  write(*,*) "Input the number of time steps N_t: " 
  read(*,*)  n_t
  !n_t = 10000

  END SUBROUTINE input_mass_and_stepsizes
!------------------------------------------------------------------------------
subroutine input_initial_conditions(x1_0, y1_0, x2_0, y2_0, vx1_0, vy1_0, vx2_0, vy2_0)
  real(DP) :: x1_0, y1_0
  real(DP) :: x2_0, y2_0
  real(DP) :: vx1_0, vy1_0
  real(DP) :: vx2_0, vy2_0 

  ! All commented values are for circular orbit of object 1 case
  write(*,*) "Input the initial x and y positions of object 1: " 
  read(*,*)  x1_0, y1_0
  !x1_0 = 1.0
  !y1_0 = 0.0

  write(*,*) "Input the initial x and y positions of object 2: " 
  read(*,*)  x2_0, y2_0
  !x2_0 = 0.0
  !y2_0 = 1.0
  
  write(*,*) "Input the intial velocity in the x and y direction of object 1:" 
  read(*,*)  vx1_0, vy1_0
  !vx1_0 = 0.0 
  !vy1_0 = sqrt(10.0)

  write(*,*) "Input the intial velocity in the x and y direction of object 2:"
  read(*,*)  vx2_0, vy2_0
  !vx2_0 = 0.0 
  !vy2_0 = 0.0
  
  END SUBROUTINE input_initial_conditions
  !--------------------------------------------------------------------------
subroutine input_file(filename2) 
  character (len=100) :: filename2

  write(*,*) "Input the name of the output file to used for graphing : " 
  read(*,*) filename2
  !filename2 = 'orbitadvanced1.dat'
  
  END SUBROUTINE input_file
  !-----------------------------------------------------------------------------------------
  end module user_input 
