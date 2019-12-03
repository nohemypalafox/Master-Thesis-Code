!Evolucion diferencial
module DE_mod
  use rand_mod
  implicit none

  private

  public :: DE_Method

contains

  !JD = "J" <--- v<riante jiter; JD = "D"  <-- variante dither
  subroutine DE_Method(m,gmax,bL,bU,F,Cr,JD,fob,mejor,minimo,tiempo)
    integer, intent(in) :: m, gmax
    real(8), intent(in), dimension(:) :: bL, bU
    real(8), intent(in) :: F, Cr
    character(len = 1), intent(in) :: JD
    real(8), intent(out), dimension(size(bL)) :: mejor
    real(8), intent(out) :: minimo, tiempo

    interface 
       function fob(x) result(f_r)
         real(8), intent(in), dimension(:) :: x
         real(8) :: f_r
       end function fob
    end interface

    real(8), dimension(m, size(bL)) :: Xmat, Vmat, Umat, Smat
    integer :: i
    real(8) :: t1, t2

    call cpu_time(t1)

    Xmat = X0(m,bL,bU)

    do i=1, gmax
       Vmat = V(Xmat, bL, bU, F, JD)
       Umat = U(Xmat,Vmat,Cr)
       Smat = Sop(Xmat,Umat,fob)
       Xmat = Smat
    end do

    call Best(Xmat, fob, mejor, minimo)

    call cpu_time(t2)

    tiempo = t2-t1

  end subroutine DE_Method


  ! operador que extrae el mejor individuo y el minimo
  subroutine Best(X, fob, mejor, minimo)
    real(8), intent(in), dimension(:,:) :: X
    interface 
       function fob(x) result(f_r)
         real(8), intent(in), dimension(:) :: x
         real(8) :: f_r
       end function fob
    end interface

    real(8), intent(out), dimension(size(X,2)) :: mejor
    real(8), intent(out) :: minimo

    real(8), dimension(size(X,1)) :: vecvals
    real(8), dimension(size(X,2)) :: xi
    integer :: i, m, ibest
    integer, dimension(1) :: ibestvec

    m = size(X,1)

    do i = 1, m
       xi = X(i,:)
       vecvals(i) = fob(xi)
    end do

    minimo = minval(vecvals)
    ibestvec = minloc(vecvals) 
    ibest = ibestvec(1)
    mejor = X(ibest,:)

  end subroutine Best

  ! operador de seleccion
  function Sop(X, U, fob) result(Smat)
    real(8), intent(in), dimension(:,:) :: X, U

    interface 
       function fob(x) result(f_r)
         real(8), intent(in), dimension(:) :: x
         real(8) :: f_r
       end function fob
    end interface

    real(8), dimension(size(X,1),size(X,2)) :: Smat

    real(8), dimension(size(X,2)) :: xvec, uvec
    integer :: i, m

    m = size(X,1)
    Smat = U

    do i = 1, m
       xvec = X(i,:)
       uvec = U(i,:)

       if(fob(xvec) <= fob(uvec)) Smat(i,:) = xvec

    end do

  end function Sop

  !Mutacion
  function V(X, bL, bU, F, JD) result(Vmat)
    real(8), intent(in), dimension(:,:) :: X
    real(8), intent(in), dimension(size(X,2)) :: bL, bU
    real(8), intent(in) :: F
    character(len = 1), intent(in) :: JD

    real(8), dimension(size(X,1),size(X,2)) :: Vmat
    integer :: m, n, i, j, r1, r2, r3
    integer, dimension(3) :: r
    real(8) :: rd

    m = size(X,1) 
    n = size(X,2)

    do i = 1, m
       r = randperm(m,3)
       r1 = r(1)
       r2 = r(2)
       r3 = r(3)

       rd = rand01()
       do j = 1, n
          if (JD == "J" .or. JD == "j") &
               Vmat(i,j) = X(r1,j) + F*rand01()*(X(r2,j) - X(r3,j)) !jitter

          if (JD == "D" .or. JD == "d") &
               Vmat(i,j) = X(r1,j) + F*rd*(X(r2,j) - X(r3,j)) !dither

          if(Vmat(i,j) < bL(j) .or. Vmat(i,j) > bU(j)) &
               Vmat(i,j) = bL(j) + rand01()*(bU(j)-bL(j))
       end do

    end do
  end function V

  !Cruza
  !Cr <--- crossover probability
  function U(X, V, Cr) result(Umat)
    real(8), intent(in), dimension(:,:) :: X, V
    real(8), intent(in) :: Cr
    real(8), dimension(size(X,1),size(X,2)) :: Umat

    integer :: i, j, m, n, jrand
    real(8) :: r

    m = size(X,1)
    n = size(X,2)

    Umat = X
    do i = 1, m
       do j = 1, n
          r = rand01()
          jrand = irand(n)
          if(r <= Cr .or. j == jrand) Umat(i,j) = V(i,j)
       end do
    end do
  end function U

  ! para generar poblacion inicial
  function X0(m, bL, bU) result(f_r)
    integer, intent(in) :: m
    real(8), intent(in), dimension(:) :: bL, bU
    real(8), dimension(m,size(bL)) :: f_r
    integer :: i, j, dimbL

    dimbL = size(bL)
    do i = 1,m
       do j = 1, dimbL
          f_r(i,j) = bL(j) + rand01()*(bU(j)-bL(j))
       end do
    end do
  end function X0

end module DE_mod
