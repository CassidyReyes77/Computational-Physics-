module functions
  use nrtype
  use constants
  use input
  implicit none

contains
  !-----------------------------------------------------------------------------
  subroutine diagonal_v(n, v, v_matrix)
    real(DP), dimension(0:n-1)  , intent(in) :: v
    real(DP), dimension(0:n-1,0:n-1)             :: v_matrix
    integer                              :: i, j
    integer                 , intent(in) :: n 

    do i=0, n-1
       v_matrix(i,i) = v(i)
    do j = 0, n-1
       if (i/=j) then
          v_matrix(i,j) = 0.0_DP
       end if
    end do
 end do

!  do i=0, n-1
!     v_matrix(i,i) = v(i)
!     ! Optionally set off-diagonal elements to 0.0_DP if they are not already initialized
!  end do

  end subroutine diagonal_v
!-------------------------------------------------------------------------------
  subroutine kintetic_energy(n, h, mass, t_matrix)
     real(DP), dimension(0:n-1,0:n-1)        :: t_matrix        !
     real(DP),                intent(in)     :: mass            !
     real(DP),                intent(in)     :: h               !
     integer                                 :: i
     integer ,                intent(in)     :: n
     
     do i=0, n-1
        t_matrix(i,i) =  2.0 * h_bar*h_bar/(2.0*mass*h*h)
        if (i>0) then
           t_matrix(i, i-1) = -1.0 * h_bar*h_bar/(2.0*mass*h*h)
        end if
        if (i < n-1) then 
           t_matrix(i, i+1) = -1.0 * h_bar*h_bar/(2.0*mass*h*h)
        end if
     end do

   end subroutine kintetic_energy
   !-------------------------------------------------------------------------------
   subroutine Hamiltonian(n, v_matrix, t_matrix, h_matrix)
     real(DP), dimension(0:n-1,0:n-1), intent(in) :: v_matrix        
     real(DP), dimension(0:n-1,0:n-1), intent(in) :: t_matrix        
     real(DP), dimension(0:n-1,0:n-1)             :: h_matrix        
     integer                              :: i, j
     integer,                  intent(in) :: n 
     
     write(*,*) 'h_matrix:'
     do i = 0, n-1
        do j = 0, n-1
           h_matrix(i,j) = v_matrix(i,j) + t_matrix(i,j)
           write(*,'(F10.4)', advance='no') h_matrix(i,j)
        end do
        write(*,*)  ! New line after each row
     end do

   end subroutine Hamiltonian
   !----------------------------------------------------------------------------
   subroutine initialize_psi(n, x, mu, sigma, psi)
     real(DP),                   intent(in) :: mu              
     real(DP),                   intent(in) :: sigma           
     complex(DP), dimension(0:n-1)          :: psi
     real(DP), dimension(0:n-1), intent(in) :: x   
     integer                                :: i
     integer ,                   intent(in) :: n 

     do i=0, n-1
        psi(i) = exp(-((x(i) - mu)**2) / (2*sigma**2)) / (sigma*sqrt(2*pi))
     end do 

   end subroutine initialize_psi
   !----------------------------------------------------------------------------
   subroutine cn_setup(h_matrix, cn_plus_matrix, cn_minus_matrix, n, dt)
     real   (DP),                         intent(in) :: dt
     real   (DP), dimension(0:n-1,0:n-1), intent(in) :: h_matrix        !
     complex(DP), dimension(0:n-1,0:n-1)             :: cn_plus_matrix  !
     complex(DP), dimension(0:n-1,0:n-1)             :: cn_minus_matrix !
     integer,                             intent(in) :: n               !
     integer                                         :: i, j            !

     do i = 0, n-1
        do j = 0, n-1
           if (i == j) then
              cn_plus_matrix(i,j) = 1.0/(1.0_DP + i_complex * (dt / (2.0_DP * h_bar)) * h_matrix(i,j))
              cn_minus_matrix(i,j) = 1.0_DP - i_complex * (dt / (2.0_DP * h_bar)) * h_matrix(i,j)
           else
              cn_plus_matrix(i,j) = 1.0/(i_complex * (dt / (2.0_DP * h_bar)) * h_matrix(i,j))
              cn_minus_matrix(i,j) = -i_complex * (dt / (2.0_DP * h_bar)) * h_matrix(i,j)
           end if
        end do
     end do

   end subroutine cn_setup
   !-------------------------------------------------------------------------------
   subroutine crank_nicholson(cn_plus_matrix, cn_minus_matrix, psi_next, psi, n)
   complex(DP), dimension(0:n-1,0:n-1), intent(in) :: cn_plus_matrix  !
   complex(DP), dimension(0:n-1,0:n-1), intent(in) :: cn_minus_matrix !
   complex(DP), dimension(0:n-1)                   :: psi_next        !
   complex(DP), dimension(0:n-1),       intent(in) :: psi             !
   integer                                         :: i, j            !
   integer,                             intent(in) :: n               !
   
    ! do i=0, n-1
    !    do j=0, n-1
    !       psi_next(i) = cn_minus_matrix(i,j) * cn_plus_matrix(i,j) * psi(i)
    !    end do
   ! end do

   do i = 0, n-1
      psi_next(i) = complex(0.0_DP, 0.0_DP)  ! Initialize to zero
      do j = 0, n-1
         psi_next(i) = psi_next(i) + cn_minus_matrix(i,j) * cn_plus_matrix(i,j) * psi(j)
      end do
   end do
   
     end subroutine crank_nicholson 
!--------------------------------------------------------------------------------------
subroutine standard_deviation_x(psi, x, n, dx, std_dev)
    use nrtype
    implicit none

    complex(DP), dimension(0:n-1), intent(in) :: psi                 ! Wave function array
    real   (DP), dimension(0:n-1), intent(in) :: x                   ! Position array
    integer    ,                   intent(in) :: n                   ! Number of points
    real   (DP),                   intent(in) :: dx                  ! Spacing between points
    real   (DP)                               :: std_dev             ! Standard deviation
    real   (DP)                               :: expectation_x       ! Expectation value of x
    real   (DP)                               :: expectation_x2      ! Expectation value of x squared
    real   (DP)                               :: probability_density
    integer                                   :: j

    expectation_x = 0.0_DP
    expectation_x2 = 0.0_DP

    ! Calculate expectation values
    do j = 0, n-1
        probability_density = (abs(psi(j)))**2
        expectation_x = expectation_x + x(j) * probability_density * dx
        expectation_x2 = expectation_x2 + x(j)**2 * probability_density * dx
    end do

    ! Calculate standard deviation
    std_dev = sqrt(expectation_x2 - expectation_x**2)
end subroutine standard_deviation_x

 end module functions

 
