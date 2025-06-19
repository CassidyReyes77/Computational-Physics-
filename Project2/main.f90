! Programmer: Cassidy Reyes

! Purpose:
! The purpose of this program is to compare booles method with the monte carlo method
! in calculating the triple integral that gives us the total flux at the detector a
! distance (x0,y0,z0) away from a plutonium reactor of height H, width W, and depth D. 
! 
! Modules
!------------------------------
! constants
! input
! integration
  
! Integration module:
!-------------------------------------------------
!   d_flux(r_n, pi, x0, y0, z0, x, y, z)
!   boole_array(df, h, n)
!   integration_3d(d_flux,nx,ny,nz,r_n,pi,x0,y0,z0,hx,hy,hz,dz,g,f,flux_boole)
!   monte_carlo_3d(d_flux,flux_monte, d, w, h, r_n, x0, y0, z0, nx, ny, nz)
 
! Variables:       Type              Units         Description
  !  nx        | integer          |    N/A       | Number of gridpoints in x direction
  !  ny        | integer          |    N/A       | Number of gridpoints in y direction 
  !  nz        | integer          |    N/A       | Number of gridpoints in z direction 
  !  hx        | real(8)          |    cm        | Step size in the x direction 
  !  hy        | real(8)          |    cm        | Step size in the y direction 
  !  hz        | real(8)          |    cm        | Step size in the z direction 
  !  dz        | real(8) 3D array | cm^-5s^-1    | Calculated d_flux integrand for all space 
  !  g         | real(8) 2D array | cm^-4s^-1    | xy_slices in the z direction
  !  f         | real(8) 1D array | cm^-3s^-1    | x pipes in the y direction 
  !  h         | real(8)          |    cm        | Height of the reactor            
  !  w         | real(8)          |    cm        | Width of the reactor 
  !  d         | real(8)          |    cm        | Depth of the reactor  
  !  r_n_input | real(8)          | neutrons/g/s | User input neutron emission rate of the reactor 
  !  r_n       | real(8)          | s*cm^-3      | Neutron emission rate adjusted with the density of plutonium
  !  density   | real(8)          | g*cm^-3      | The density of plutonium 
  !  x0        | real(8)          |    cm        | x location of detector relative to the reactor
  !  y0        | real(8)          |    cm        | y location of detector relative to the reactor 
  !  z0        | real(8)          |    cm        | z location of detector relative to the reactor 
  !  x         | real(8)          |    cm        | x distance from a point in space to the detector 
  !  y         | real(8)          |    cm        | y distance from a point in space to the detector 
  !  z         | real(8)          |    cm        | z distance from a point in space to the detector 
  ! flux_boole | real(8)          | cm^-2s^-1    | Calculated flux using booles method in 3D
  ! flux_monte | real(8) array    | cm^-2s^-1    | Calculated flux using monte carlos method in 3D
  ! pt_flux    |                  | cm^-2s^-1    | Flux at the detector due to a point source reactor 
  ! uncertainty| real(8) array    | cm^-2s^-1    | Calculated uncertainty of monte carlos method 
  !            |                  |              |            
!-----------------------------------------------------------------------------------------------------------
program boole
  use constants
  use integration
  use input
  implicit none

  integer :: nx, ny, nz
  real(dbl), dimension(:,:,:), allocatable :: dz
  real(dbl), dimension(:,:), allocatable :: g
  real(dbl), dimension(:), allocatable :: f
  real(dbl) :: hx, hy, hz, y0, x0, z0, r_n, d, w, h
  real(dbl) :: flux_boole, flux_monte, pt_flux
  real(dbl) :: r_n_input, density, uncertainty 

   print *,"3D Boole vs 3D Monte Carlo Integration "
   print *,"==========================="

   call input_variables(r_n_input,r_n,density,h,w,d,x0,y0,z0,nx,ny,nz)

! Calculate the step size for each x, y, z direction based on the dimensions of the reactor and number of pts
   hx = d/nx
   hy = w/ny
   hz = h/nz

   call integration_3d(d_flux,nx,ny,nz,r_n,pi,x0,y0,z0,hx,hy,hz,dz,g,f,flux_boole)
   call monte_carlo_3d(d_flux,flux_monte, d, w, h, r_n, x0, y0, z0, nx, ny, nz)

! Calculate the flux at the detector a distance (x0,y0,z0) away from a point reactor source
! As r or (x0,y0,z0) gets very large for our original reactor, the flux should approach this value
   pt_flux =  d_flux(r_n, pi, x0, y0, z0, d/2.0, w/2.0, h/2.0)

! Print results to screen    
   print *,"The boole flux is", flux_boole, "cm^-2s^-1"
   print *, "The monte carlo flux is", flux_monte, "cm^-2s^-1"
   print *, "The statistical uncertainty of the flux found with the monte carlo method is:", uncertainty  
   print *, "The flux for very large r =", pt_flux, "cm^-2s^-1"
end program boole 
