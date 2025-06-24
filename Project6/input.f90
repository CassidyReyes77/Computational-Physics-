module input
  use nrtype
  implicit none

contains
!-------------------------------------------------------------------------------
  subroutine user_input(L, mu, sigma, dt, N, Nt)
    real(DP) :: L
    real(DP) :: mu
    real(DP) :: sigma, sigma2
    real(DP) :: dt
    integer  :: N
    integer  :: Nt
    

    write(*,*) 'Input the size of the box, L:'
    read(*,*) L
    ! L=

    write(*,*) 'Input the number of lattice points N'
    read(*,*) N

    write(*,*) 'Input the number of time steps Nt'
    read(*,*) Nt

    write(*,*) 'Input the time step dt'
    read(*,*) dt

    write(*,*) 'Input the center of the distribution mu'
    read(*,*) mu
    ! mu=

    write(*,*) 'Input the width of the initial distribution 2 * sigma'
    read(*,*) sigma2
    ! 2_sigma =
    sigma = 0.5*sigma2
   
  end subroutine user_input
  !-----------------------------------------------------------------------------
end module input


  

    

    
