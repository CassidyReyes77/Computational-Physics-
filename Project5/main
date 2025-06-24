program Woods_saxon_potential
  use nrtype
  use nrutil
  use tqli_
  use tred2_
  use sort_
  use constants 
  use user_input
  use functions
  implicit none

  real(DP)                              :: m          !MeV
  real(DP)                              :: v0         !MeV
  real(DP)                              :: a          !fm
  real(DP)                              :: R          !fm
  real(DP)                              :: r1         !fm
  real(DP)                              :: r2         !fm
  real(DP)                              :: L          !fm
  real(DP)                              :: L2         !fm
  real(DP)                              :: h          !step size
  real(DP)                              :: h2         !step size
  real(DP), allocatable, dimension(:)   :: x          !poistion?
  real(DP), allocatable, dimension(:)   :: x2         !poistion?
  real(DP), allocatable, dimension(:)   :: radius     !changing radius
  real(DP), allocatable, dimension(:)   :: d          !poistion
  real(DP), allocatable, dimension(:)   :: e          !poistion
  real(DP), allocatable, dimension(:)   :: f          !poistion
  real(DP), allocatable, dimension(:)   :: g          !poistion
  real(DP), allocatable, dimension(:)   :: E0         !Ground state energy
  real(DP), allocatable, dimension(:)   :: E1         !First energy level
  real(DP), allocatable, dimension(:)   :: E2         !Second energy level
  real(DP), allocatable, dimension(:)   :: E3         !Third energy level
  real(DP), allocatable, dimension(:)   :: v          !1D potential array
  real(DP), allocatable, dimension(:)   :: v2         !1D potential array
  real(DP), allocatable, dimension(:,:) :: v_matrix   !2D potential, diagonal matrix
  real(DP), allocatable, dimension(:,:) :: t_matrix   !2D Kinetic energy
  real(DP), allocatable, dimension(:,:) :: h_matrix   !2D Total energy
  real(DP), allocatable, dimension(:,:) :: v_matrix2  !2D potential, diagonal matrix
  real(DP), allocatable, dimension(:,:) :: t_matrix2  !2D Kinetic energy
  real(DP), allocatable, dimension(:,:) :: h_matrix2  !2D Total energy
  integer                               :: n          !number of grid points 
  integer                               :: i,j,k,ii,jj
  character(len=100)                    :: filename   !Name of the file that stores the r vs. energy data 

  call input(m, v0, a, R, L, n, filename)
  
  allocate(x(0:n), v(0:n), v_matrix(0:n,0:n))
  allocate(t_matrix(0:n,0:n), h_matrix(0:n,0:n))
  allocate(x2(0:n), v2(0:n), v_matrix2(0:n,0:n))
  allocate(t_matrix2(0:n,0:n), h_matrix2(0:n,0:n))
  allocate(d(0:n), e(0:n), f(0:n), g(0:n), radius(0:n))
  allocate(E0(0:n), E1(0:n), E2(0:n), E3(0:n))

  open(unit=10, file=trim(filename), status='unknown', action='write')
  
  h = 2.0*L/(n-1) 

  x_vs_potential : do i=0, n-1
     x(i) = -L + real(i) * h
     v(i) = potential(v0, x(i), R, a)
  end do x_vs_potential 
  
  diagonal_v : do i=0, n-1
     v_matrix(i,i) = v(i)
     do j = 0, n-1
        if (i/=j) then
           v_matrix(i,j) = 0.0_DP
        end if
     end do
  end do diagonal_v

  kinetic_energy : do i=0, n-1
     t_matrix(i,i) =  2.0 * h_bar*h_bar/(2.0*m*h*h)
     if (i>0) then
        t_matrix(i, i-1) = -1.0 * h_bar*h_bar/(2.0*m*h*h)
     end if
     if (i < n-1) then 
        t_matrix(i, i+1) = -1.0 * h_bar*h_bar/(2.0*m*h*h)
     end if
  end do kinetic_energy

  ! Construct and print h_matrix
  write(*,*) 'h_matrix:'
  do i = 0, n-1
     do j = 0, n-1
        h_matrix(i,j) = v_matrix(i,j) + t_matrix(i,j)
        write(*,'(F10.4)', advance='no') h_matrix(i,j)
     end do
     write(*,*)  ! New line after each row
  end do

  !Solve tridiagoal eigenvalue problem 
  call tred2(h_matrix, d, e)
  call tqli(d, e, h_matrix)

  ! Assuming 'd' is your array of eigenvalues
  ! Sort the array 'd'
  call sort(d)  ! 'n' is the size of your array

  ! Print the first three energy levels
  print *, 'The first four energy levels are:'
  print *, 'E0 = ', d(0)
  print *, 'E1 = ', d(1)
  print *, 'E2 = ', d(2)
  print *, 'E3 = ', d(3)
!-------------------------------------------------------------------------------
  !Re-do calculations and write radius vs. energy to another file
  r1   = 2.0
  r2   = 10.0
  r_vs_energy : do i= 0, n-1
     radius(i) = r1 + real(i) * (r2-r1)/(n-1)
     L2 = 2.0*radius(i)
     h2 = 2.0*L2/(n-1)
     
     recalculate_x_v: do j =0, n-1    
        x2(j) = -L2 + real(j) * h2
        v2(j) = potential(v0, x2(j), radius(j), a)  
     end do recalculate_x_v

     re_diagonal_v : do k=0, n-1
        v_matrix2(k,k) = v2(k)
        do ii = 0, n-1
           if (k/=ii) then
              v_matrix2(k,ii) = 0.0_DP
           end if
        end do
      end do re_diagonal_v

     kinetic_energy2 : do jj=0, n-1
        t_matrix2(jj,jj) =  2.0 * h_bar*h_bar/(2.0*m*h2*h2)
        if (jj>0) then
           t_matrix2(jj, jj-1) = -1.0 * h_bar*h_bar/(2.0*m*h2*h2)
        end if
        if (jj < n-1) then 
           t_matrix2(jj, jj+1) = -1.0 * h_bar*h_bar/(2.0*m*h2*h2)
        end if
     end do kinetic_energy2

     h_matrix2 = v_matrix2 + t_matrix2

     call tred2(h_matrix2, f, g)
     call tqli(f, g, h_matrix2)
     call sort(f)

     E0(i) = f(0)
     E1(i) = f(1)
     E2(i) = f(2)
     E3(i) = f(3)      
     
  end do r_vs_energy

  do i= 0, n-1
     write(10,*) radius(i), E0(i), E1(i), E2(i), E3(i)
  end do

  write(*,*) 'Data for r vs. energy written to file: ', trim(filename)

  deallocate(x,v,v_matrix,t_matrix,h_matrix,d,e,f,g,radius)
  deallocate(x2,v2,v_matrix2,t_matrix2,h_matrix2)
  deallocate(E0, E1, E2, E3)

  close(10)

  end program Woods_saxon_potential
