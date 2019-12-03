
!Parser module
module parser_mod
  implicit none

  private

  public :: nan_hack, nargs, operf, operf2, parsv, filterx

contains

  !this function extract a subset of x which corresponds to a valid
  !function
  function filterx(x,xval) result(f_r)
    integer, intent(in), dimension(:) :: x
    real(8), intent(in) ::  xval
    integer, allocatable, dimension(:) :: f_r
    integer :: k, cont
    real(8) :: zero

    zero = tiny(zero)
    zero = 10.d0*zero !a little margin of 1 order of magnitude

    cont = 0
    do k = 1, size(x)
       cont = cont + 1 
       if( abs(parsv(x,xval) - parsv(x(1:k),xval)) .le. zero ) exit
    end do

    allocate(f_r(cont))
    f_r = x(1:cont)
  end function filterx

  !parsv(x,xval) evaluates the function associated to x (or associated to
  !a subset of x) in xval.
  !
  function parsv(x,xval) result(f_r)
    integer, intent(in), dimension(:) :: x
    real(8), intent(in) ::  xval
    real(8)  :: f_r
    integer :: length, k
    real(8), dimension(size(x)) :: auxv

    length = size(x)

    auxv = nan_hack() 
    f_r = auxv(1)

    !evaluating the "functions" of zero arguments. We have chosen the x 
    !variable to have zero arguments. If x is defined as the identity 
    !function, this identity function would have 1 argument
    do k = 1, length
       if(x(k) .ge. 1 .and. x(k) .le. 10) auxv(k) = 1.0d0*x(k)
       if(x(k) .eq. 11) auxv(k) = xval
    end do

    f_r = auxv(1)

    ! evaluating functions of 1 argument and 2 arguments
    do k = length - 1, 1, -1
       if(.not. isnan( auxv(k+1) ) .and. nargs(x(k)) .eq. 0) then
          cycle
       elseif(.not. isnan(auxv(k+1)) .and. nargs(x(k)) .eq. 1) then
          auxv(k) = operf(x(k),auxv(k+1))
          if( isnan(auxv(k)) ) return
       elseif(.not. isnan( auxv(k+1) ) .and. nargs(x(k)).eq. 2)then
          auxv(k) = operf2(x(k),auxv(k+1:length))
          if( isnan(auxv(k)) ) return
       end if
    end do

    f_r = auxv(1)
  end function parsv


  ! function that evaluates the functions of two arguments,
  !i.e., +,*,/,^. vecval is an array with real numbers, and possibly
  !nan components this function operates on the first two real no-nan
  !components of the vecval array and then transforms them into nan
  !components 
  function operf2(x,vecval) result(f_r)
    integer, intent(in) :: x
    real(8), intent(inout), dimension(:) :: vecval
    real(8) :: f_r
    real(8), dimension(2) :: auxiliar
    integer :: j, cont, length

    length = size(vecval)
    auxiliar = nan_hack()

    cont = 0
    do j = 1, length
       if( .not. isnan(vecval(j)) ) then 
          cont = cont + 1
          auxiliar(cont) = vecval(j)
          vecval(j) = nan_hack()
       end if
       if(cont .eq. 2) exit
    end do

    if(x .eq. 12) f_r = auxiliar(1) + auxiliar(2) 
    if(x .eq. 13) f_r = auxiliar(1) * auxiliar(2)
    if(x .eq. 14) f_r = auxiliar(1) / auxiliar(2)
    if(x .eq. 15) f_r = auxiliar(1) ** auxiliar(2)

  end function operf2

  !function that evaluates the functions of one argument at the point 
  !xval changing it to NaN 
  function operf(x,xval) result(f_r)
    implicit none
    integer, intent(in) :: x
    real(8), intent(inout) :: xval
    real(8) ::  f_r

    if(x .eq. 16) f_r = -xval
    if(x .eq. 17) f_r = sin(xval)
    if(x .eq. 18) f_r = cos(xval)
    if(x .eq. 19) f_r = exp(xval)
    if(x .eq. 20) f_r = log(xval)
    if(x .eq. 21) f_r = sqrt(xval)
    if(x .eq. 22) f_r = tan(xval)
    if(x .eq. 23) f_r = asin(xval)
    if(x .eq. 24) f_r = acos(xval)
    if(x .eq. 25) f_r = atan(xval)
    if(x .eq. 26) f_r = erf(xval)
    if(x .eq. 27) f_r = sinh(xval)
    if(x .eq. 28) f_r = cosh(xval)
    if(x .eq. 29) f_r = tanh(xval)

    xval = nan_hack()

  end function operf


  ! funcion that calculates the number of arguments
  function nargs(x) result(f_r)
    integer, intent(in) :: x
    integer :: f_r

    if(x .ge.  1 .and. x .le. 11) f_r = 0
    if(x .ge. 12 .and. x .le. 15) f_r = 2
    if(x .ge. 16 .and. x .le. 29) f_r = 1

  end function nargs

  !A NaN implementation. If the compiler supports IEEE arithmetic use:
  !use, intrinsic :: ieee_arithmetic
  !real(8) :: NaN
  !NaN = ieee_value(0d0, ieee_quiet_nan)
  function nan_hack() result(f_r)
    real(8) :: f_r
    character(len=3) :: nanstr = "nan"

    read(nanstr,*) f_r
  end function nan_hack

end module parser_mod

