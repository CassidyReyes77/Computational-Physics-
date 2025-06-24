program main
  use nrtype
  use constants
  use input
  use functions
  implicit none

   real   (DP)                              :: L               !
   real   (DP)                              :: mu              !
   real   (DP)                              :: sigma           !
   real   (DP)                              :: dt              !
   real   (DP)                              :: dx
   real   (DP)                              :: standard_deviation
   real   (DP), allocatable, dimension(:)   :: x               !
   real   (DP), allocatable, dimension(:)   :: v               !
   complex(DP), allocatable, dimension(:)   :: psi             ! Wave function array
   complex(DP), allocatable, dimension(:)   :: psi_next        !
   real   (DP), allocatable, dimension(:,:) :: v_matrix        !
   real   (DP), allocatable, dimension(:,:) :: t_matrix        !
   real   (DP), allocatable, dimension(:,:) :: h_matrix        !
   complex(DP), allocatable, dimension(:,:) :: cn_plus_matrix  !
   complex(DP), allocatable, dimension(:,:) :: cn_minus_matrix !
   integer                                  :: n               !
   integer                                  :: nt              !
   integer                                  :: i,j             !

   call user_input(L, mu, sigma, dt, n, nt)

   allocate(x(0:n-1), v(0:n-1), psi(0:n-1), psi_next(0:n-1))
   allocate( v_matrix(0:n-1,0:n-1), t_matrix(0:n-1,0:n-1), h_matrix(0:n-1,0:n-1))
   allocate(cn_plus_matrix(0:n-1,0:n-1), cn_minus_matrix(0:n-1,0:n-1))

   open(unit=10, file='standard_deviation.dat', status='replace')

   dx = 2.0*L/(n-1)

   x_vs_potential : do i=0, n-1
      x(i) = -L + real(i) * dx
      v(i) = 0.0
   end do x_vs_potential
   
! This routine takes the potential array 'v' and makes a diagonal matrix 
   call diagonal_v(n, v, v_matrix)

! Print potential energy matrix
!   write(*,*) 'Potential Energy Matrix v_matrix:'
!   do i = 0, n-1
!       do j = 0, n-1
!           write(*,'(F10.4)', advance='no') v_matrix(i,j)
!       end do
!       write(*,*)  ! New line after each row
!   end do   

! Calculate the kinetic energy by applying the second derivative central difference formula to schrodingers equation 
   call kintetic_energy(n, dx, m, t_matrix)

! Adds the potential and kinetic energy to get the Hamiltonian matrix 
   call Hamiltonian(n, v_matrix, t_matrix, h_matrix)

! Print Hamiltonian matrix
   write(*,*) 'Hamiltonian Matrix h_matrix:'
   do i = 0, n-1
       do j = 0, n-1
           write(*,'(F10.4)', advance='no') h_matrix(i,j)
       end do
       write(*,*)  ! New line after each row
   end do

! Initialize psi 
   call initialize_psi(n, x, mu, sigma, psi)
   
 write(*,*) 'psi:'
    do i = 0, n-1
        write(*,'(A, F10.4, A, F10.4, A)') '(', real(psi(i)), ', ', aimag(psi(i)), 'i)'
    end do
     
! Calculate the 
   call cn_setup(h_matrix, cn_plus_matrix, cn_minus_matrix, n, dt)

! Main loop for time evolution
   do i = 0, nt-1
      call crank_nicholson(cn_plus_matrix, cn_minus_matrix, psi_next, psi, n)

      ! Compute necessary quantities
      call standard_deviation_x(psi_next, x, n, dx, standard_deviation) ! Implement this

      ! Write to files
      write(10, *) i*dt, standard_deviation

    ! Swap psi and psi_next for next iteration 
     end do
  
deallocate(x, v, v_matrix, t_matrix, h_matrix, cn_plus_matrix, cn_minus_matrix, psi_next)
close(10) 
     
 end program main
 
