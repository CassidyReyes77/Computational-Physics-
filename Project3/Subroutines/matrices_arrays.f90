module matrices_arrays
  use constants
  implicit none

contains
!-------------------------------------------------------------------------------
  subroutine calculate_A_error(A, g, BE_error, row_number)
    integer, intent(in) :: row_number
    real(dbl), dimension(row_number,l), intent(in) :: g
    real(dbl), dimension(l,l), intent(out) :: A
    real(dbl), dimension(row_number), intent(in) :: BE_error
    integer :: j, k

    A = 0.0
    do j=1, l
       do k=1, l
          A(j,k) = sum(g(:,j)*g(:,k)*1.0/(BE_error(1:row_number)**2.0))
       end do
    end do

  end subroutine calculate_A_error
  !-------------------------------------------------------------------------------
    subroutine calculate_A(A, g)
    real(dbl), dimension(m,l), intent(in) :: g
    real(dbl), dimension(l,l), intent(out) :: A
    integer :: j, k

    A = 0.0
    do j=1, l
       do k=1, l
          A(j,k) = sum(g(:,j)*g(:,k))
       end do
    end do

  end subroutine calculate_A
  !-------------------------------------------------------------------------------
  subroutine calculate_single_g(g_row, n_i, z_i)
    real(dbl), dimension(l) :: g_row
    real(dbl) :: n_i, z_i

       g_row(1) = n_i + z_i
       g_row(2) = (n_i + z_i)**(2.0/3.0)
       g_row(3) = (n_i - z_i)**2.0/(n_i + z_i)
       g_row(4) = z_i * (z_i - 1.0)/((n_i+z_i)**(1.0/3.0))

       if (n_i + z_i > 0.0) then
          g_row(5) = 1.0/((n_i + z_i)**(3.0/4.0))
       else
          write(*,*) "Error: Invalid input for g(5). n(i)+z(i) must be >0"
          stop
       end if
       
       ! Kronecker delta 
       if (mod(n_i, 2.0) == 0 .and. mod(z_i, 2.0) == 0) then
          g_row(5) = abs(g_row(5))
       elseif (mod(n_i, 2.0) == 1.0 .and. mod(z_i, 2.0) == 1.0) then
          g_row(5) = -abs(g_row(5))
       elseif (mod(n_i + z_i, 2.0) == 1.0) then
          g_row(5) = 0.0
       else
          g_row(5) = g_row(5)
       end if

end subroutine calculate_single_g
!------------------------------------------------------------------------------
subroutine calculate_g(g,n,z,row_number)
  integer, intent(in) :: row_number
  real(dbl), dimension(row_number,l) ::g
  real(dbl), dimension(l) :: g_row
  real(dbl), dimension(row_number) :: n, z
  integer :: i

  g = 0.0
  do i=1, row_number 
     call calculate_single_g(g_row,n(i),z(i))
     g(i,:) = g_row
  end do
  
  end subroutine calculate_g
!------------------------------------------------------------------------------
  SUBROUTINE calculate_b_error(b, g, BE_experimental, BE_error, row_number)
  integer, intent(in) :: row_number
  real(dbl), dimension(row_number, l), intent(in) :: g
  real(dbl), dimension(row_number), intent(in) :: BE_experimental, BE_error
  real(dbl), dimension(l), intent(out) :: b

  integer :: i, j
  real(dbl) :: summ

  b=0.0
  do j = 1, l
     summ = 0.0
     do i = 1, row_number
        summ = summ + g(i,j) * BE_experimental(i) * 1.0/(BE_error(i)**2.0)
     end do
     b(j) = summ
  end do
  
end subroutine calculate_b_error
!------------------------------------------------------------------------------
SUBROUTINE calculate_b(b, g, BE_experimental)
  real(dbl), dimension(m, l), intent(in) :: g
  real(dbl), dimension(m), intent(in) :: BE_experimental
  real(dbl), dimension(l), intent(out) :: b  
  INTEGER :: i, j
  REAL(dbl) :: summ

  b=0.0
  do j = 1, l
     summ = 0.0
     do i = 1, m
        summ = summ + g(i,j) * BE_experimental(i)
     end do
     b(j) = summ
  end do
  
end subroutine calculate_b
!------------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  subroutine uncertainty_coefficients(inverse, diagonal, uncertainty)
    real(dbl), dimension(l,l), intent(in) :: inverse   
    real(dbl), dimension(l) :: diagonal
    real(dbl), dimension(l), intent(out) :: uncertainty
    integer :: i 

    do i = 1, l
       diagonal(i) = inverse(i,i)
       uncertainty(i) = sqrt(diagonal(i))
    end do

  end subroutine uncertainty_coefficients

  !-------------------------------------------------------------------------------------------
subroutine uncertainty_single_BE(n_input,z_input,uncertainty_c_alpha,g_row,c_alpha,uncertainty_BE)
  real(dbl), intent(in) :: n_input, z_input
  real(dbl) :: h, derivative_sum
  real(dbl), dimension(l), intent(in) :: g_row, c_alpha, uncertainty_c_alpha
  real(dbl), dimension(l) :: BE_plus, BE_minus
  real(dbl), intent(out) :: uncertainty_BE
  
  ! Take the derivative of the function BE
  h = 1000.0
  call BE_equation(n_input, z_input, g_row, c_alpha + h, BE_plus)
  call BE_equation(n_input, z_input, g_row, c_alpha - h, BE_minus)
  derivative_sum = sum(BE_plus(1:l) - BE_minus(1:l) / (2.0*h))
  uncertainty_BE = sqrt(derivative_sum**2.0 * sum(uncertainty_c_alpha(1:l)))

end subroutine uncertainty_single_BE
!-------------------------------------------------------------------------------
subroutine BE_equation(n_input, z_input, g_row, c_alpha, BE)
  real(dbl) :: n_input, z_input 
  real(dbl), dimension(l) :: g_row, c_alpha, BE
  integer :: i

  call calculate_single_g(g_row, n_input, z_input) 
  do i=1, l    
     BE(i) = c_alpha(i) * g_row(i)
  end do

  end subroutine BE_equation 
!------------------------------------------------------------------------------ 

end module matrices_arrays
