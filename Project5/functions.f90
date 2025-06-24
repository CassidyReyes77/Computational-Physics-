module functions
  use nrtype
  use user_input 
  implicit none

contains
!--------------------------------------------------------------------------------------
  real(DP) function potential(v0,x,R,a)
    real(DP), intent(in) :: v0
    real(DP), intent(in) :: x
    real(DP), intent(in) :: R
    Real(DP), intent(in) :: a 

    potential = v0 / (1 + exp((abs(x)-R)/a))

  end function potential
!--------------------------------------------------------------------------------------
  
  end module functions 
