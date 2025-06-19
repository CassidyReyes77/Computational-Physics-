! Programmer: Cassidy Reyes
!
! Purpose:
! Find the coefficients of the nuclear semi-empirical mass formula by using experimental data from a file the
! user provides, and using methods such as lu decomposition and single variable decomposition to solve the
! system of linear equations.
!
! In part 1 we read in the first 800 rows of data, construct matrices based on the
! nuclear semi empirical mass formula and input data, perform lu decompisition on the lxl A matrix, apply back
! substitution to the result of the lu decmposition to get the soltution vector c_alpha1, print results to the
! screen, and finally we ask the user for z and n to calculate the binding energy based on our calcuted solution
! c_alpha1.
!
! In part 2 we read all rows of data, construct matrices based on the semi empirical mass formula and
! and the input data (including the error this time), perform single variable decomposition on a lxl A matrix, 
! apply back substitution to the result to get the solution vector c_alpha2, find the uncertainty in each 
! coefficient, we ask the user for z and n to calculate the binding energy based on our calculated solution, and
! finally we find the uncertainty of our calculated binding energy and print to the screen. 

! Notes:
! I'm not sure that my A matrix is correct for part 2 because my solution vector isn't correct or close to
! correct (its off by at least an order of 10 for each coefficient) and also because the calculated binding
! energy isn't entirely correct either (it's of the same order but not close to the correct value). Also I'm not
! sure that I calculated the uncertainty correctly because the values for both seem ridiculously large. 
! 
 
!              Variables:       Type              Units         Description
  !  z1                    | integer          |    N/A       | atomic number part 1
  !  n1                    | integer          |    N/A       | neutron number part 1 
  !  z2                    | integer          |    N/A       | atomic number part 2
  !  n2                    | real(8)          |    cm        | neutron number part 2
  !  z_input1              | real(8)          |    cm        | User input atomic number part 1 
  !  n_input1              | real(8)          |    cm        | User input neutron number part 1
  !  z_input2              | real(8) 3D array | cm^-5s^-1    | User input atomic number part 2
  !  n_input2              | real(8) 2D array | cm^-4s^-1    | User input neutron number part 2
  !  BE_experimental_kev_1 | real(8) 1D array | cm^-3s^-1    | Experimental binding energy from the file before
  !                        |                  |              | conversion part 1
  !  BE_error_kev_1        | real(8)          |    cm        | Experimental binding energy error from file
  !                        |                  |              | before conversion part 1  
  !  BE_experimental_1     | real(8)          |    cm        | Experimental binding energy after converting from
  !                        |                  |              | kev to eV  
  !  BE_experimental_kev_2 | real(8)          |    cm        | Experimental binding energy from the file before
  !                        |                  |              | conversion part 2
  !  BE_error_kev_2        | real(8)          | neutrons/g/s | Experimental binding energy error from file
  !                        |                  |              | before conversion part 2 
  !  BE_experimental_2     | real(8)          | s*cm^-3      | Experimental binding energy after converting from
  !                        |                  |              | kev to eV 
  !  BE_error_2            | real(8)          | g*cm^-3      | Experimental bindingenergy error after converting
  !                        |                  |              | from kev to eV  
  !  BE_calculated_1       | real(8)          |    cm        | The binding energy found by multiplying the 
  !                        |                  |              | solution by the the g array (found with the user 
  !                        |                  |              | input z and n) for part 1
  !  BE_calculated_2       | real(8)          |    cm        | The binding energy found by multiplying the 
  !                        |                  |              | solution by the the g array (found with the user 
  !                        |                  |              | input z and n) for part 2
  !  w_min                 | real(8)          |    cm        | definition for really small 
  !  w_max                 | real(8)          |    cm        | The max singular value 
  !  g1                    | real(8)          |    cm        | g array that stores the values of the portion of 
  !                        |                  |              | the semi emipirical mass equation that's a
  !                        |                  |              | function of z,n for all combinations of z,n part1
  !  g2                    | real(8)          |    cm        | g array that stores the values of the portion of 
  !                        |                  |              | the semi emipirical mass equation that's a
  !                        |                  |              | function of z,n for all combinations of z,n part2
  !  b1                    | real(8)          | cm^-2s^-1    | right hand side vector in equation to be solved
  !                        |                  |              | that's equal to g array * binding energy part 1
  !  b2                    | real(8) array    | cm^-2s^-1    | right hand side vector in equation to be solved
  !                        |                  |              | that's equal to g array * binding energy part 1
  !  c_alpha1              |                  | cm^-2s^-1    | Solution vector that holds coeffients of the
  !                        |                  |              | semi empirical mass formula part 1
  !  c_alpha2              | real(8) array    | cm^-2s^-1    | Solution vector that holds coeffients of the
  !                        |                  |              | semi empirical mass formula part 2
  !  g_row                 | integer          |    N/A       | The portion of the semi empirical mass equation
  !                        |                  |              | that's a function of z and n (1 iteration) 
  !  column                | integer          |    N/A       | A column of A matrix part 2's inverse 
  !  diagonal              | integer          |    N/A       | The diagonal element of A inverse 
  !  uncertainty_c_alpha   | real(8)          |    cm        | Array to store uncertainty of each coeffiecient
  !                        |                  |              | found in part 2
  !  w                     | real(8)          |    cm        | singular values 
  !  A1_pt1                | real(8)          |    cm        | A matrix for part 1
  !  A1_pt2                | real(8) 3D array | cm^-5s^-1    | A matrix for part 2
  !  A2_pt2                | real(8) 2D array | cm^-4s^-1    | Saved A matrix for part 2
  !  A3_pt2                | real(8) 1D array | cm^-3s^-1    | Saved A matrix for part 2
  !  matrix_A1             | real(8)          |    cm        | A maxtrix part 1           
  !  matrix_A2             | real(8)          |    cm        | A matrix part 2 
  !  inverse               | real(8)          |    cm        | inverse matrix of A matrix in part 2
  !  identity              | real(8)          | neutrons/g/s | An identity matrix 
  !  filename              | real(8)          | s*cm^-3      | The file name of the data to be read in 
  !  row_number            | real(8)          | g*cm^-3      | The number of rows in the data file 
  !  indx                  | real(8)          |    cm        | array to store pivoted values of A matrix part 1
  !  lndx                  | real(8)          |    cm        | array to store pivoted values of A matrix part 2
  
!-----------------------------------------------------------------------------------------------------------

program main
    use constants
    use nrtype
    use nrutil
    use matrices_arrays
    use ludcmp
    use lubksb
    use svdcmp
    use svbksb
    implicit none
    
    real(dbl) :: n_input_1, z_input_1, n_input_2, z_input_2
    real(dbl) :: BE_calculated_1, BE_calculated_2, uncertainty_BE
    real(dbl) :: w_min, w_max
    real(dbl), dimension(l,l) :: A1_pt1, A1_pt2, A2_pt2, A3_pt2
    real(dbl), dimension(l,l) :: matrix_A1, matrix_A2, inverse, identity
    real(dbl), dimension(m,l) :: g1, g2 
    real(dbl), dimension(m) :: z1, n1, z2, n2
    real(dbl), dimension(m) :: BE_experimental_kev_1, BE_error_kev_1, BE_experimental_1
    real(dbl), dimension(m) :: BE_experimental_kev_2, BE_error_kev_2, BE_experimental_2, BE_error_2  
    real(dbl), dimension(l) ::  c_alpha1, c_alpha2, b1, b2
    real(dbl), dimension(l) :: g_row, column, diagonal, uncertainty_c_alpha, w  
    integer(I4B) :: i, j, row_number
    character (len=100) :: file_name
    real(dbl) :: d, q 
    integer(I4B), dimension(l) :: indx, lndx
    
    write(*,*) "********************Part 1****************************************************"
    
    ! Ask for file name and read file 
    write(*,*) "Input the name of the data file to be used in this program:"
    read(*,*) file_name
    open(10, file=trim(file_name))
    read(10,*)
    do i=1, m
       read(10,*) z1(i), n1(i), BE_experimental_kev_1(i), BE_error_kev_1(i)
       BE_experimental_1(i) =  BE_experimental_kev_1(i) * 1000.0
    end do

    ! Calculate matrices for part 1
    call calculate_g(g1, n1, z1, m)
    call calculate_b(b1, g1, BE_experimental_1)
    call calculate_A(matrix_A1, g1)

    ! Print A matrix for part 1 and save for later 
    write(*,*) "A matrix part 1:"
    do i =1, l
       write(*,'(9E14.6)') (matrix_A1(i,j), j=1, l)
    end do
    A1_pt1(1:l,1:l) = matrix_A1(1:l,1:l)
    
    ! Do lu decomposition for part 1
    call ludcmp_(A1_pt1(1:l,1:l), indx(1:l), d)

    ! solve system of equations for part 1 and print solution
    c_alpha1(1:l) = b1(1:l)
    call lubksb_(A1_pt1(1:l,1:l),indx(1:l),c_alpha1(1:l))   
      
    ! Print solution and test results for part 1
    write(*,*) "Solution vector/c_alpha part 1 in eV:"
    write(*,'(9E14.6)') c_alpha1(1:l)

    ! Ask user for N and Z of nuclei part 1
    print *, "What are the values Z and N of the nuclei?"
    read(*,*) z_input_1, n_input_1

    ! Calculate the binding energy based on user input N and Z part 1
    call calculate_single_g(g_row, n_input_1, z_input_1)
    BE_calculated_1 = dot_product(c_alpha1(1:l), g_row(1:l))
    write(*,*) "The binding energy for part 1 calculated using the solution vector c_alpha is:"
    write(*,'(E14.6, A3)')  BE_calculated_1, "eV"
!------------------------------------------------------------------------------------------------
    ! Part 2
    write(*,*) "********************Part 2****************************************************"
    read(10,*) row_number 
    do i=1, row_number 
       read(10,*) z2(i), n2(i), BE_experimental_kev_2(i), BE_error_kev_2(i)
       BE_experimental_2(i) =  BE_experimental_kev_2(i) * 1000.0
       BE_error_2(i) =  BE_error_kev_2(i) * 1000.0
    end do

    call calculate_g(g2, n2, z2, row_number)
    call calculate_b_error(b2, g2, BE_experimental_2, BE_error_2, row_number)
    call calculate_A_error(matrix_A2, g2, BE_error_2, row_number)

    ! Print A matrix for part 2 and save matrices for later 
    write(*,*) "A matrix part 2:"
    do i=1, l
       write(*,'(9E14.6)') (matrix_A2(i,j), j=1, l)
    end do
    A1_pt2(1:l,1:l) = matrix_A2(1:l,1:l)
    A2_pt2(1:l,1:l) = matrix_A2(1:l,1:l)
    
    ! Decompose A matrix and find the max singular value 
    w = 0.0_dp
    call svdcmp_dp(A1_pt2(1:l,1:l), w(1:l), A3_pt2(1:l,1:l))
    w_max = max(maxval(w(1:l)), 0.0_sp)

    ! Define small and zero the "small" singular values
    w_min = w_max*(1.0e-6)
    where (w(1:l) < w_min) w(1:l) = 0.0  

    ! backsubstitute for each right-hand side vector and solve the system 
    do j=1,m
       c_alpha2(1:l)=b2(1:l)
       call svbksb_dp(A1_pt2(1:l,1:l),w(1:l),A3_pt2(1:l,1:l),b2(1:l),c_alpha2(1:l))
    end do

      ! Make identity matrix
    identity = 0.0
    do i = 1, l
       identity(i,i) = 1.0
    end do

    ! Get A inverse
    inverse = 0.0
    call ludcmp_(A2_pt2(1:l,1:l), lndx(1:l), q)
    do j = 1, l
       column = identity(:,j)
       call lubksb_(A2_pt2(:,:), lndx(:), column)
       inverse(:,j) = column
    end do

    ! Extract diagonal elements of the inverse of matrix_A and take the square root
    ! to get the uncertainties in the coefficients
    call uncertainty_coefficients(inverse, diagonal, uncertainty_c_alpha)
    
    ! Print solution, uncertainty, and test results for part 2 
    write(*,*) 'Solution vector/c_alpha part 2 in eV:'
    write(*,'(9E14.6)') (c_alpha2(i), i=1,l)
    write(*,*) "Uncertainty of the coeffiencts in eV:"
    write(*,'(9E14.6)') (uncertainty_c_alpha(i), i=1,l)

    ! Ask user for N and Z of nuclei part 2
    print *, "What are the values Z and N of the nuclei?"
    read(*,*) z_input_2, n_input_2

    ! Calculate the binding energy based on user input N and Z part 1
    call calculate_single_g(g_row, n_input_2, z_input_2)
    BE_calculated_2 = dot_product(c_alpha2(1:l), g_row(1:l))
    write(*,*) "The binding energy for part 2 based on the solution vector c_alpha is:"
    write(*,'(E14.6, A3)')  BE_calculated_2, "eV"

    ! Print the uncertainty of the binding energy
    call uncertainty_single_BE(n_input_2,z_input_2,uncertainty_c_alpha,g_row,c_alpha2,uncertainty_BE)
      write(*,*) "The uncertainty of the calculated binding energing is:"
      write(*,'(E14.6, A3)')  uncertainty_BE, "eV" 

       
    close (10)

  end program main
