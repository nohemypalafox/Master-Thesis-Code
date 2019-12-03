module RK4_mod
    implicit none
    
    private
    
    public :: RK4
    
    contains 
    
    !Fvec = rhs de un sistema de ED1
    function RK4(Fvec,t0,tf,x0,npts) result(matf_r)
    
    interface
    function Fvec(t,x) result(Fv_r)
      real(8), intent(in) :: t
      real(8), intent(in), dimension(:) :: x
      real(8), dimension(size(x)) :: Fv_r
    end function Fvec
    
    
    end interface
    
    real(8), intent(in) :: t0, tf
    real(8), intent(in), dimension(:) :: x0
    integer, intent(in) :: npts
    real(8), dimension(npts + 1, size(x0) + 1) :: matf_r
    
    real(8), dimension(npts + 1) :: tvec
    real(8), dimension(npts + 1, size(x0)) :: Xmat 
    real(8) :: h, ti, timas1
    real(8), dimension(size(x0)) :: xi, ximas1, k1vec, k2vec,k3vec,k4vec
    integer :: i, ncol
    
    h = (tf -t0)/npts
  
    ncol = size(x0) + 1
    
    tvec(1) = t0
    
    Xmat(1,:) = x0
    do i= 1, npts
        ti = tvec(i)
        
        xi = Xmat(i,:)
        
        k1vec = Fvec(ti,xi)
        k2vec = Fvec(ti + 0.5d0*h,xi + 0.5d0*h*k1vec)
        k3vec = Fvec(ti + 0.5d0*h,xi + 0.5d0*h*k2vec)
        k4vec = Fvec(ti + h,xi + h*k3vec)
        
        timas1 = ti + h
        ximas1 = xi + 1d0/6d0*(k1vec + 2d0*k2vec + 2d0*k3vec + k4vec)*h
        tvec(i+1) = timas1
        Xmat(i+1,:) = ximas1
    end do
    
    matf_r(:,1) = tvec
    matf_r(:, 2:ncol) = Xmat
    
    end function RK4
    
end module RK4_mod
