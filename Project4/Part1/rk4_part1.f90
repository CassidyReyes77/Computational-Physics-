module runge_kutta
  use nrtype
  use constants
  use user_input
  implicit none

contains
  !-----------------------------------------------------------------------------
  real(DP) function func(i, m_prim, y0) result(dfdt)
    real(DP), intent(in) :: m_prim
    real(DP), dimension(n1) :: y0
    real(DP) :: r1
    integer :: i
 
    r1 = r(y0(1),y0(2))
    
    select case(i) 
 case (1)
    dfdt = y0(3)                                !y0(3) = vx
 case (2)
    dfdt = y0(4)                                !y0(4) = vy
 case (3)
    dfdt = velocity_eq_1(m_prim, r1, y0(1))     !y0(1) = x
 case (4) 
    dfdt = velocity_eq_1(m_prim, r1, y0(2))     !y0(2) = y
 end select
 
end function func
!-----------------------------------------------------------------------------
  subroutine rk4_step(y0, h, m, y4)                                !Fourth Order Runge Kutta solver
    real(DP), dimension(n1),       intent(in)  :: y0               !Array of initial values
    real(DP),                      intent(in)  :: h                !Step size
    real(DP),                      intent(in)  :: m                !Mass of primary object 
    real(DP), dimension(size(y0)), intent(out) :: y4               !Array of solutions 
    real(DP), dimension(size(y0))              :: k1, k2, k3, k4   !Intermediate k values 
    integer :: i

    do i = 1, size(y0)
       k1(i) = h*func(i, m, y0)
       k2(i) = h*func(i, m, y0+k1(i)/2.)
       k3(i) = h*func(i, m, y0+k2(i)/2.)
       k4(i) = h*func(i, m, y0+k3(i))
       y4(i) = y0(i) + (k1(i) + 2.*k2(i) + 2.*k3(i) + k4(i)) / 6.0
    end do

    end subroutine rk4_step 
!-------------------------------------------------------------------------------
       real(DP) function velocity_eq_1(m_prim, r1, x1)              !Velocity equation for one component 
      real(DP), intent(in) :: m_prim, x1, r1

      velocity_eq_1 = -G*m_prim*x1/(r1**3.0)

    end function velocity_eq_1
!-------------------------------------------------------------------------------
    real(DP) function energy_eq_1(m_prim,m_1,x,y,vx1,vy1)
      real(DP), intent(in) :: m_prim, m_1
      real(DP), intent(in) :: x, y         
      real(DP), intent(in) :: vx1, vy1      
      real(DP)             :: ke            !Kinetic energy 
      real(DP)             :: r1

      ke = 0.5 * m_1 * (vx1**2.0 + vy1**2.0)
      r1 = r(x, y)

      energy_eq_1 = ke - G*m_1*m_prim/r1    !Kinetic energy - Potential energy

    end function energy_eq_1
!-------------------------------------------------------------------------------
       real(DP) function r(x,y)
      real(DP), intent(in) :: x, y

         r = sqrt(x**2.0 + y**2.0)

    end function r
!-------------------------------------------------------------------------------
    end module runge_kutta 
