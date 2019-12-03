  ! para compilar
  ! ifort -heap-arrays 1024 SARS_Optimiza.f90
  ! la opcion '-heap-arrays 1024' es necesaria por las grandes dimensiones
  ! de los arreglos. gfortran no lo requiere
  ! en mi computadora (Intel(R) Core(TM) i5-4460  CPU @ 3.20GHz)
  ! tarda como hora y cuarto
  
  !mimimimimimimimimimimimimimimimimimimimimimi
  include "../../modules/RK4_mod.f90"
  include "../../modules/rand_mod.f90"
  include "../../modules/DE_mod.f90"


  ! modulo parametros
  ! 
  !
  ! lambda          = tasa de reclutamiento
  ! beta, epsE,     = coeficientes de transmision desde la clase S al resto.
  ! epsQ, epsJ
  ! mu              = tasa de mortalidad (natural)
  ! p               = tasa de reclutamiento de asintomaticos
  ! k1              = tasa de transferancia entre clases E y I
  ! k2              = tasa de transferancia entre clases Q y J
  ! d1, d2          = tasas de mortalidad inducidas por la enfermedad
  ! sigma1, sigma2  = tasas de recuperacion
  module par_mod
    implicit none
    real(8), parameter :: t0 = 0d0, tf = 360d0
    integer, parameter :: npts = 360, dimf = 6
    real(8), parameter :: S0 = 12d6, E0 = 1565d0, Q0 = 292d0, I0 = 695d0
    real(8), parameter :: J0 = 326d0, R0 = 20d0
    real(8), parameter, dimension(dimf) :: x0 = [S0, E0, Q0, I0, J0, R0]

    real(8), parameter :: beta = 0.2d0, epsE = 0.3d0, epsQ = 0d0, &
         epsJ = 0.1d0, mu = 0.000034d0, p = 0d0, k1 = 0.1d0,      &
         k2 = 0.125d0, d1 = 0.0079d0, d2 = 0.0068d0,              &
         sigma1 = 0.0337d0, sigma2 = 0.0386d0

    real(8), parameter :: B1 = 1d0, B2 = 1d0, B3 = 1d0, B4 = 1d0, &
         C1 = 300d0, C2 = 600d0

  end module par_mod


  ! modulo para fob
  module fob_mod
    use par_mod
    use RK4_mod

    implicit none
    private
    public :: fob

  contains

    function fob(x) result(fob_r)
      real(8), intent(in), dimension(:) :: x
      real(8) :: fob_r
      integer :: i
      real(8), dimension(size(x)/2) :: r1
      real(8), dimension(size(x) - size(x)/2) :: r2

      real(8), dimension(npts + 1) :: ct, cE, cQ, cI, cJ, cu1, cu2
      real(8), dimension(npts + 1) :: integ1, integ2, integ
      real(8), dimension(npts + 1, 2) :: matinteg
      real(8), dimension(npts + 1, dimf + 1) :: Amat

      integer :: dimx, dimr1, dimr2

      dimx = size(x)
      dimr1 = dimx/2
      dimr2 = dimx - dimr1
      r1 = x(1:dimr1) 
      r2 = x(dimr1+1:dimr2)

      Amat = RK4(Fvec,t0,tf,x0,npts)

      ct = Amat(:,1)
      cE = Amat(:,3)
      cQ = Amat(:,4)
      cI = Amat(:,5)
      cJ = Amat(:,6)

      integ1 = B1*cE + B2*cQ + B3*cI + B4*cJ

      do i =1, npts+1
         cu1(i) = u1(ct(i),r1)
         cu2(i) = u2(ct(i),r2)
      end do

      integ2 = C1/2d0*cu1**2 + C2/2d0*cu2**2
      integ = integ1 + integ2

      matinteg(:,1) = ct
      matinteg(:,2) = integ

      fob_r = trapz(matinteg)

    contains

      function Fvec(t, xvec) result(vecf_r)
        real(8), intent(in) :: t
        real(8), intent(in), dimension(:) :: xvec
        real(8), dimension(size(xvec)) :: vecf_r
        real(8) :: fS, fE, fQ, fI, fJ, fR, xS, xE, xQ, xI, xJ, xR        

        real(8) :: Npob, lambda

        xS = xvec(1)
        xE = xvec(2)
        xQ = xvec(3)
        xI = xvec(4)
        xJ = xvec(5)
        xR = xvec(6)

        Npob = sum(xvec)
        lambda = mu*Npob

        fS = lambda - (1d0 / Npob) * xS * (beta*xI + epsE*beta*xE + &
             epsQ*beta*xQ + epsJ*beta*xJ) - mu * xS

        fE = p + (1d0 / Npob) * xS * (beta*xI + epsE*beta*xE + &
             epsQ*beta*xQ + epsJ*beta*xJ) - (u1(t,r1) + k1 + mu) * xE

        fQ = u1(t,r1)*xE - (k2 + mu) * xQ
        fI = k1 * xE - (u2(t,r2) + d1 + sigma1 + mu) * xI
        fJ = u2(t,r2)*xI + k2 * xQ - (d2 + sigma2 + mu) * xJ
        fR = sigma1 * xI + sigma2 * xJ - mu * xR

        vecf_r = [fS, fE, fQ, fI, fJ, fR]

      end function Fvec

    end function fob


    !control 1 (cuarentena)
    function u1(t,r) result(u1_r)
      real(8), intent(in) :: t
      real(8), intent(in), dimension(:) :: r
      real(8) :: u1_r
      real(8) :: dt, eps
      integer :: idt

      eps = epsilon(eps)

      dt = (tf - t0)/dble(size(r))
      idt = floor(t/dt) + 1

      if(abs(t-tf)<=eps) idt = size(r)

      u1_r = r(idt) 

    end function u1

    ! control 2 (aislamiento)
    function u2(t,r) result(u2_r)
      real(8), intent(in) :: t
      real(8), intent(in), dimension(:) :: r
      real(8) :: u2_r
      real(8) :: dt, eps
      integer :: idt

      eps = epsilon(eps)
      dt = (tf - t0)/size(r)
      idt = floor(t/dt) + 1

      if(abs(t-tf)<=eps) idt = size(r)      

      u2_r = r(idt) 

    end function u2

    ! trapz (dt constante)
    function trapz(A) result(f_r)
      real(8), intent(in), dimension(:,:) :: A
      real(8) :: f_r
      real(8) :: dt
      real(8), dimension(size(A,1)) :: t, y
      real(8), dimension(size(A,1)-2) :: ycut
      integer :: ny


      t = A(:,1)
      y = A(:,2)
      ny = size(A,1)

      dt = t(2)-t(1)

      ycut = y(2:(ny-1))

      f_r = dt/2d0*(y(1) + y(ny)) + dt*sum(ycut)
    end function trapz

  end module fob_mod



  !---------                    Main Program                   ---------!
  !---------!---------!---------!---------!---------!---------!---------!
  Program prueba
    use RK4_mod
    use fob_mod
    use DE_mod
    !use par_mod

    implicit none

    integer, parameter :: dimx = 180*4, m = 1000, gmax = 10000
    real(8), parameter, dimension(dimx) :: bL = 0.05d0, bU = 0.5d0
    real(8), parameter :: F = 1d0, Cr = 0.9d0
    character(len = 1), parameter :: JD = "D" 
    real(8), dimension(dimx) :: mejor
    real(8) :: minimo, tiempo

    real(8), dimension(dimx/2) :: u1vec
    real(8), dimension(dimx - dimx/2) :: u2vec
    integer :: i


!    call random_seed()
    
    
    call DE_Method(m,gmax,bL,bU,F,Cr,JD,fob,mejor,minimo,tiempo)

    u1vec = mejor(1:dimx/2)
    u2vec = mejor(dimx/2+1:dimx)

    print*, minimo
    print*, 'tiempo(s) =', tiempo
    
    open(unit = 11, file = "u1.dat", status = "unknown")
    open(unit = 12, file = "u2.dat", status = "unknown")
    
    do i = 1, dimx/2
       write(11,"(f0.4)") u1vec(i)
    end do
    close(unit = 11)

    do i = 1, dimx - dimx/2
       write(12,"(f0.4)") u2vec(i)
    end do
    close(unit = 12)

  end program prueba





