module runge_kutta
  use nrtype
  use constants
  use user_input
  implicit none

contains
  !-----------------------------------------------------------------------------
  real(DP) function func(i, m_prim, m_1, m_2, y0) result(dfdt)
    real(DP), intent(in) :: m_prim, m_1, m_2
    real(DP), dimension(n2), intent(in) :: y0
    real(DP) :: r1, r2, r12
    integer :: i
 
    r1  = r(y0(1),y0(2))
    r2  = r(y0(3),y0(4))
    r12 = r((y0(1)-y0(3)), (y0(2)-y0(4)))
    
    select case(i) 
 case (1)
    dfdt = y0(5)                                !y0(5) = vx1
 case (2)
    dfdt = y0(6)                                !y0(6) = vy1
 case (3)
    dfdt = y0(7)                                !y0(7) = vx2
 case (4)                                    
    dfdt = y0(8)                                !y0(8) = vy2
 case (5)
    dfdt = velocity_eq_2(m_prim, m_2, r1, r12, y0(1), y0(3))     !y0(1) = x1
 case (6) 
    dfdt = velocity_eq_2(m_prim, m_2, r1, r12, y0(2), y0(4))     !y0(2) = y1
 case (7)
    dfdt = velocity_eq_2(m_prim, m_1, r2, r12, y0(3), y0(1))     !y0(3) = x2
 case (8) 
    dfdt = velocity_eq_2(m_prim, m_1, r2, r12, y0(4), y0(2))     !y0(4) = y2
 end select
 
end function func
!-----------------------------------------------------------------------------
  subroutine rk4_step(y0, h, m_prim, m_1, m_2, y4)                 !Fourth Order Runge Kutta solver
    real(DP), dimension(n2),       intent(in)  :: y0               !Array of initial values
    real(DP),                      intent(in)  :: h                !Step size
    real(DP),                      intent(in)  :: m_prim, m_1, m_2 !Mass  
    real(DP), dimension(size(y0)), intent(out) :: y4               !Array of solutions 
    real(DP), dimension(size(y0))              :: k1, k2, k3, k4   !Intermediate k values 
    integer :: i

    do i = 1, size(y0)
       k1(i) = h*func(i, m_prim, m_1, m_2, y0)
       k2(i) = h*func(i, m_prim, m_1, m_2, y0+k1(i)/2.)
       k3(i) = h*func(i, m_prim, m_1, m_2, y0+k2(i)/2.)
       k4(i) = h*func(i, m_prim, m_1, m_2, y0+k3(i))
       y4(i) = y0(i) + (k1(i) + 2.*k2(i) + 2.*k3(i) + k4(i)) / 6.0
    end do

    end subroutine rk4_step 
!-------------------------------------------------------------------------------
       real(DP) function velocity_eq_2(m_prim, m_2, r1, r12, x1, x2)     !Velocity equation for one component 
         real(DP), intent(in) :: m_prim, m_2
         real(DP), intent(in) :: r1, r12
         real(DP), intent(in) :: x1, x2
         real(DP)             :: velocity_eq_1

      velocity_eq_1 = -G*m_prim*x1/(r1**3.0)
      velocity_eq_2 =  velocity_eq_1 - G*m_2*(x1-x2)/(r12**3.0)

    end function velocity_eq_2
!-------------------------------------------------------------------------------
    real(DP) function energy_eq_2(m_prim, m_1, m_2, x1, y1, x2, y2, vx1, vy1, vx2, vy2)
      real(DP), intent(in) :: m_prim, m_1, m_2
      real(DP), intent(in) :: x1, y1, x2, y2  
      real(DP), intent(in) :: vx1, vy1, vx2, vy2
      real(DP)             :: r1, r2, r12
      real(DP)             :: ke1           !Kinetic energy of object 1
      real(DP)             :: ke2           !Kinetic energy of object 2
      real(DP)             :: pe            !Total Potential energy 


      ke1 = 0.5 * m_1 * (vx1*vx1 + vy1*vy1)
      ke2 = 0.5 * m_2 * (vx2*vx2 + vy2*vy2)
      
      r1 = r(x1, y1)
      r2 = r(x2, y2)
      r12 =r((x2-x1), (y2-y1))
      
      pe  = -G*m_1*m_prim/ r1 -G*m_2*m_prim/r2 -G*m_1*m_2/r12

      energy_eq_2 = ke1 + ke2 + pe   

    end function energy_eq_2
!-------------------------------------------------------------------------------
    real(DP) function r(x,y)
      real(DP), intent(in) :: x, y

      !if (x == 0.0 .and. y == 0.0) then
     !    r = tiny
     ! else
         r = sqrt(x*x + y*y)
     ! end if 

    end function r
!-------------------------------------------------------------------------------
    end module runge_kutta 
