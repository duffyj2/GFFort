module shared_data
implicit none
  ! Mathematical parameters
  complex(8), parameter :: im = (0.0d0,1.0d0)
  real(8), parameter :: pi = acos(-1.0d0)
  ! Material parameters
  real(8), parameter :: U = 10.0d0
  real(8), parameter :: n0 = 1.0d0
  real(8), parameter :: E_F = 0.0d0
  real(8), parameter :: eta = 1.0d-4
  ! Variables
  real(8) :: ex_split, mag_m
  integer :: nE, DZ

end module shared_data


program main
use shared_data
implicit none
  real(8) :: band_shift
  real(8) :: mag_temp
  real(8) :: n_up, n_down
  real(8), parameter :: tol = 1.0d-4
  ! RootFinder
  real(8) :: a, b, eps, eps_fxn, h
  integer :: ifail
  
  nE = 8
  DZ = 6
  
  h = 0.1d0
  eps = 1.0d-4
  eps_fxn = 0.0d0
  ifail = 1

  mag_m = 0.5d0
  band_shift = 0.0d0
  ex_split = U*mag_m
  do
    mag_temp = mag_m
    call C05AGF(band_shift,h,eps,eps_fxn,f,a,b,ifail)
    mag_m = n_up - n_down
    if ( abs(mag_m-mag_temp) .le. tol ) exit
    ex_split = U*mag_m
    write(*,*) mag_m, band_shift
  end do
  
  write(*,*) mag_m, band_shift
  
  contains
  
    function f(x)
      real(8) :: f
      real(8), intent(in) :: x
      
      call occupation(n_up,n_down,x)
      f = n_up + n_down - n0
    end function

end program main




subroutine occupation(n_up, n_down, band_shift)
use GF_AGNR
use shared_data
implicit none
  real(8), intent(out) :: n_up, n_down
  real(8), intent(in) :: band_shift
  real(8) :: V_up, V_down
  ! Integration parameters
  integer, parameter :: lw=800, liw=lw/4
  real(8) :: w(lw)
  integer :: iw(liw)
  real(8) :: abserr, lim, epsabs, epsrel
  real(8) :: res_up, res_down
  integer :: ifail, inf
  external d01amf

  V_up = band_shift - ex_split/2.0d0
  V_down = band_shift + ex_split/2.0d0
  
  epsabs = 0.0d0
  epsrel = 1.0d-4
  lim = 1.0d-4
  inf = 1		! The range is from "lim" to infinity
  ifail = 1		! You could set this to -1 to suppress error messages
  
  call d01amf(g_up,lim,inf,epsabs,epsrel,res_up,abserr,w,lw,iw,liw,ifail)
  call d01amf(g_down,lim,inf,epsabs,epsrel,res_down,abserr,w,lw,iw,liw,ifail)
  
  n_up = 1.0d0/2.0d0 + res_up/pi 
  n_down = 1.0d0/2.0d0 + res_down/pi 
  
  contains
  
    function g_up(E)
      real(8) :: g_up
      real(8), intent(in) :: E
      complex(8) :: g
      
      g = GF(nE,0,DZ,DZ,'bb',im*E)
      g_up = real( g/(1.0d0-g*V_up) )
    end function g_up
    
    function g_down(E)
      real(8) :: g_down
      real(8), intent(in) :: E
      complex(8) :: g
      
      g = GF(nE,0,DZ,DZ,'bb',im*E)
      g_down = real( g/(1.0d0-g*V_down) )
    end function g_down
  
end subroutine occupation