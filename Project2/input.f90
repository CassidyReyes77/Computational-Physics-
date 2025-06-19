module input
  use constants
  implicit none

contains
!---------------------------------------------------------------------------------  
  subroutine input_variables(r_n_input, r_n, density, h, w, d, x0, y0, z0, nx, ny, nz)
  implicit none
  real(dbl) :: r_n_input,r_n, density, h, w, d, x0, y0, z0
  integer, intent(out) :: nx, ny, nz

  write(*,*) 'Input the neutron emission rate of the reactor (neutrons/s/g).'
  read (*,*) r_n_input
  write(*,*) 'Input the depth, width, and height in that order.'
  read (*,*) d, w, h
  write(*,*) 'Input the validated location of the detector (x0,y0,z0).'
  read (*,*) x0, y0, z0
  write(*,*) 'Input the number of gridpoints in each direction (Nx,Ny,Nz) that is a multiple of 4.'
  read (*,*) nx, ny, nz

! Error message for user when nx, ny, or nz isnt a multiple of 4
  if (mod(nx,4) /= 0 .or. mod(ny,4) /= 0 .or. mod(nz,4) /= 0) then
    write(*,*) "Error: nx, ny, and nz must all be multiples of 4"
  end if

! Unit conversion from neutrons/s/g to s/cm^3 using the density of plutonium 
  density = 19.84
  r_n = r_n_input * density

end subroutine input_variables
!------------------------------------------------------------------------------

end module input
