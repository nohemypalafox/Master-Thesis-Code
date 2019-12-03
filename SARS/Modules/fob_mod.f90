
!modulo para las funciones objetivo
module fob_mod
  implicit none
  
  private
  
  public :: fob, Fvec
  
  
  contains
  
  function fob(r) result(f_r)
  real(8), intent(in), dimension(:) :: r
  real(8) :: x, y
  real(8) :: f_r
  
  x = r(1)
  y = r(2)
  
  f_r = x**2 + y**2
  
  end function fob
  
  
    function Fvec(t, xvec) result(vecf_r)
        real(8), intent(in) :: t
        real(8), intent(in), dimension(:) :: xvec
        real(8), dimension(size(xvec)) :: vecf_r
        
        real(8) :: f1, f2, x1, x2
        
        x1 = xvec(1)
        x2 = xvec(2)
        
        f1 = x2
        f2 = 5.0d0 * t * x1 - 2d0 * x2 + exp(- 2.0d0 * t)
        
        vecf_r = [f1, f2]
    end function Fvec

end module fob_mod
