module rand_mod
  implicit none
  
  private
  
  public :: rand01, irand, randperm

  contains
  !genera una permutacion aleatoria del 1 al m de tamaño s
  function randperm(m,s) result(rs)
  integer, intent(in) :: m, s
  integer, dimension(m) :: x
  integer, dimension(s) :: rs
  integer :: i, r
  
  do i = 1, m
    x(i) = i
  end do
  
  do i = 1, s
    r = irand(m - i + 1)
    rs(i) = x(r)
    x(r) = x(m - i + 1)
  end do
  end function randperm
  
  !para generar aleatorio entero en [1,m]
  function irand(m) result(f_r) 
  integer, intent(in) :: m
  integer :: f_r
  
  f_r = floor(m*rand01() + 1) 
  end function irand

  !genera un aleatorio en [0,1) 
  function rand01() result(f_r)
  real(8) :: f_r
  
    call random_number(f_r)
  end function rand01

end module rand_mod
