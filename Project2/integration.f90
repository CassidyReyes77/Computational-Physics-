module integration
  use constants
  implicit none

contains
!----------------------------------------------------------------------------
    real(dbl) function d_flux(r_n, pi, x0, y0, z0, x, y, z)
      real(dbl), intent(in) :: r_n, pi, x0, y0, z0
      real(dbl), intent(in) :: x, y, z
      real(dbl) :: dist_squared

       dist_squared = (x-x0)**2.0 + (y-y0)**2.0 + (z-z0)**2.0

! If statement to help control singularities near machine precision limit
! If statement prevents the flux from blowing up to infinity for small dist_squared
       if (dist_squared < 1.0d-10) then
          d_flux = 0.0
       else
          d_flux = r_n / (4.0*pi*dist_squared)
       end if
     
  end function d_flux
!-------------------------------------------------------------------------------- 
  real(dbl) function boole_array(df, h, n)
    
    real(dbl), dimension(:), allocatable :: df
    real(dbl), intent(in) :: h
    integer, intent(in) :: n
    integer :: i
    real(dbl) :: sum 

    sum = 0.0

! Take the sum of the function evaluations for df 
    do i=4, n, 4
       sum = sum + (7.0*df(i) + 32.0*df(i-1) + 12.0*df(i-2) + &
            32.0*df(i-3) + 7.0*df(i-4))
    end do

! Calculate the final boole result 
    boole_array = (2.0*h/45.0)*sum
    
  end function boole_array
!---------------------------------------------------------------------------------
  subroutine integration_3d(d_flux,nx,ny,nz,r_n,pi,x0,y0,z0,hx,hy,hz,dz,g,f,flux_boole)
    implicit none
    integer, intent(in) :: nx, ny, nz
    real(dbl), intent(in) :: r_n, pi, x0, y0, z0, hx, hy, hz
    real(dbl), intent(out) :: flux_boole
    real(dbl), dimension(:,:,:), allocatable :: dz
    real(dbl), dimension(:,:), allocatable :: g
    real(dbl), dimension(:), allocatable :: f, temp_array1, temp_array2
    real(8), external :: d_flux
    integer :: i, j, k

    allocate(dz(0:nx,0:ny,0:nz))
    allocate(temp_array1(0:nz), temp_array2(0:ny))
    allocate(g(0:nx,0:ny),f(0:nx))

    dz = 0.0D0
    g = 0.0D0
    f = 0.0D0
    
    do i=0, nx-1
      do j=0, ny-1
        ! Compute the dz integrand
        do k=0, nz-1
            dz(i,j,k) = d_flux(r_n,pi,x0,y0,z0,k*hz, k*hz, k*hz)
        end do
        ! Transfer an xy slice into a temp array and integrate using Boole
        temp_array1(:) = dz(i,j,:)
        g(i,j) = boole_array(temp_array1, hz, nz)
      end do
      ! Transfer pipes in the x direction into a temp array and integrate using Boole
      temp_array2(:) = g(i,:)
      f(i) = boole_array(temp_array2, hy, ny)
    end do
    ! Integrate f(i) using Boole to get the final result
    flux_boole = boole_array(f, hx, nx)
    
  end subroutine integration_3d
!----------------------------------------------------------------------------
subroutine monte_carlo_3d(d_flux,flux_monte, d, w, h, r_n, x0, y0, z0, nx, ny, nz)

  real(dbl) :: r, x, y, z, avg, sum_flux, sum_square, weight, temp_flux
  real(dbl) :: x_a, x_b, y_a, y_b, z_a, z_b, alpha, dist, variance, uncertainty
  real(dbl), intent(in) :: d, w, h, r_n, x0, y0, z0
  real(dbl), intent(out) :: flux_monte
  integer, intent(in) :: nx, ny, nz
  integer :: i, N
  real(8), external :: d_flux

  sum_flux = 0.0
  sum_square = 0.0

  ! Calculate a new N based on nx, ny, nz inputs 
  N = nx*ny*nz

  x_a = 0.0
  x_b = d
  y_a = 0.0
  y_b = w
  z_a = 0.0
  z_b = h
  
 alpha = 2000.0

  call random_init(repeatable = .false., image_distinct = .true.)
  
  do i=0, N-1
     
! Get a random number r that lies between 0 and 1
     call random_number(r)
     
! Use the bounds of the reactor to generate random x, y, z values based on random r
     x   = x_a + (x_b - x_a)*r
     y   = y_a + (y_b - y_a)*r
     z   = z_a + (z_b - z_a)*r

! Apply a weight function to get more accurate results
     dist = sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)
     weight = exp(-alpha*dist)

! Calculate the function for the given random values and take the sum of all fluxes 
     temp_flux = d_flux(r_n, pi, x0, y0, z0, x, y, z) 
     sum_flux = sum_flux + temp_flux
     
! Take the sum of the square of the function for calculating the uncertainty 
     sum_square = sum_square + temp_flux**2.0/N
  end do

! Calculate the average and the total flux 
  avg = sum_flux/N
  flux_monte = (x_b-x_a)*(y_b-y_a)*(z_b-z_a)*avg

! Calculate the uncertainty 
  variance = sum_square/N - avg**2.0
  uncertainty = variance/N * (x_b-x_a)**2.0*(y_b-y_a)**2.0*(z_b-z_a)**2.0

  end subroutine monte_carlo_3d
!----------------------------------------------------------------------------
end module integration
