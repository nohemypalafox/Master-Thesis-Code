
!modulo para escribir la forma analitica de f(x)
module vec_to_fun_mod
  implicit none
  private

  public :: str_func, in_strin

contains

  !this function writes the associated function to the integer vector 
  !vector v
  function str_func(v) result(f_string)
    integer, intent(in), dimension(:) :: v
    integer :: n, poss, i
    character(len = 10) :: tmp !all the functions of the used grammar 
    !has less than 7 characters.

    character(len = 10*size(v)) :: f_stringaux
    character(len = :), allocatable :: f_string

    poss = 0
    f_stringaux =''

    n = size(v)

    do i=1,n
       if (v(i) .eq. 1) then
          tmp = '1.0'
       elseif (v(i) .eq. 2) then
          tmp = '2.0'
       elseif (v(i) .eq. 3) then
          tmp = '3.0'
       elseif (v(i) .eq. 4) then
          tmp = '4.0'
       elseif (v(i) .eq. 5) then
          tmp = '5.0'
       elseif (v(i) .eq. 6) then
          tmp = '6.0'
       elseif (v(i) .eq. 7) then
          tmp = '7.0'
       elseif (v(i) .eq. 8) then
          tmp = '8.0'
       elseif (v(i) .eq. 9) then
          tmp = '9.0'
       elseif (v(i) .eq. 10)then
          tmp = '10.0'
       elseif (v(i) .eq. 11) then
          tmp = 'x'
       elseif (v(i) .eq. 12) then
          tmp = '()+()'
       elseif (v(i) .eq. 13) then
          tmp = '()*()'
       elseif (v(i) .eq. 14) then
          tmp = '()/()'
       elseif (v(i) .eq. 15) then
          tmp = '()**()'
       elseif (v(i) .eq. 16) then
          tmp = '-1.0*()'
       elseif (v(i) .eq. 17) then
          tmp = 'sin()'
       elseif (v(i) .eq. 18) then
          tmp = 'cos()'
       elseif (v(i) .eq. 19) then
          tmp = 'exp()'
       elseif (v(i) .eq. 20) then
          tmp = 'log()'
       elseif (v(i) .eq. 21) then
          tmp = 'sqrt()'
       elseif (v(i) .eq. 22) then
          tmp = 'tan()'
       elseif (v(i) .eq. 23) then    
          tmp = 'asin()'
       elseif (v(i) .eq. 24) then    
          tmp = 'acos()'
       elseif (v(i) .eq. 25) then    
          tmp = 'atan()'
       elseif (v(i) .eq. 26) then    
          tmp = 'erf()'
       elseif (v(i) .eq. 27) then    
          tmp = 'sinh()'
       elseif (v(i) .eq. 28) then    
          tmp = 'cosh()'
       elseif (v(i) .eq. 29) then    
          tmp = 'tanh()'
       endif

       poss = index(f_stringaux,'()')
       f_stringaux = trim(in_strin(f_stringaux,tmp,poss+1))

    enddo

    allocate(f_string, source = trim(f_stringaux))
  end function str_func

  !this function insert a string in another at position poss
  function in_strin(strin,str_in,poss)  result(str)
    character(len = *), intent(in) :: strin, str_in
    integer, intent(in) :: poss

    character(len = len_trim(strin) + len_trim(str_in)) :: str
    character(len = poss - 1) :: part1
    character(len = len_trim(strin) - poss +1) :: part2
    character(len = len_trim(strin)) :: strinaux

    strinaux = trim(strin)
    part1 = strinaux(:poss-1)
    part2 = strinaux(poss:)
    str = trim(trim(part1) // trim(str_in) // trim(part2))

  end function in_strin
end module vec_to_fun_mod
