module Routines
implicit none

  contains
  
  function Dyson1(g,V)
    complex(8) :: Dyson1
    complex(8), intent(in) :: g
    real(8), intent(in) :: V
    
    Dyson1 = g/(1-g*V)
  end function Dyson1
  
  
  function Dyson(g,V)
  use NumRoutines
    complex(8), dimension(:,:), intent(in) :: g
    real(8), dimension(:,:), intent(in) :: V
    complex(8), dimension(size(g,1),size(V,2)) :: Dyson

    Dyson = matmul(inv(Imx(size(g,1))-matmul(g,V)),g)
  end function Dyson
  
  
  function gBulkSubsMx(m,n,s,E)
  use GFMod
    ! Returns the GFMod matrix for two atomic sites in bulk graphene
    complex(8), dimension(2,2) :: gBulkSubsMx
    complex(8), intent(in) :: E
    integer, intent(in) :: m,n,s
    
    gBulkSubsMx(1,1) = gBulk_kZ(0,0,0,E)
    gBulkSubsMx(2,2) = gBulkSubsMx(1,1)
    gBulkSubsMx(1,2) = gBulk_kZ(m,n,s,E)
    gBulkSubsMx(2,1) = gBulkSubsMx(1,2)
  end function gBulkSubsMx
  
  
  function gGNRSubsMx(nE,m1,n1,m2,n2,s,E)
  use GFMod
    complex(8), dimension(2,2) :: gGNRSubsMx
    complex(8), intent(in) :: E
    integer, intent(in) :: nE,m1,n1,m2,n2,s
    
    gGNRSubsMx(1,1) = gRib_Arm(nE,m1,n1,m1,n1,0,E)
    gGNRSubsMx(2,2) = gRib_Arm(nE,m2,n2,m2,n2,0,E)
    gGNRSubsMx(1,2) = gRib_Arm(nE,m1,n1,m2,n2,s,E)
    gGNRSubsMx(2,1) = gGNRSubsMx(1,2)
  end function gGNRSubsMx
  
  
  function gGNRTopMx(nE,m1,n1,m2,n2,s,E)
  use GFMod
  use shared_data
    complex(8), dimension(2,2) :: gGNRTopMx
    complex(8), intent(in) :: E
    integer, intent(in) :: nE,m1,n1,m2,n2,s 
    complex(8), dimension(4,4) :: g, gtemp
    real(8), dimension(4,4) :: V
    complex(8) :: g_imp

    ! Introduce the connecting GFs
    g(1,1) = gRib_Arm(nE,m1,n1,m1,n1,0,E)
    g(2,2) = gRib_Arm(nE,m2,n2,m2,n2,0,E)
    g(1,2) = gRib_Arm(nE,m1,n1,m2,n2,s,E)
    g(2,1) = g(1,2)
    
    !Introduce the impurity GFs
    g_imp = 1/(E-eps_imp)
    g(3,3) = g_imp
    g(4,4) = g(3,3)
    
    ! The peturbation connects the impurities to the lattice
    data V(3,1),V(1,3),V(2,4),V(4,2) /4*tau/
    
    gtemp = Dyson(g,V)
    gGNRTopMx = gtemp(3:4,3:4)
  end function gGNRTopMx

  
end module Routines


module Dynamic
implicit none

  contains
  
  function XHFBulk1(Vup,Vdown,w)
  use Routines
  use GFMod
  use NumRoutines
  use shared_data
    !Calculates the Hartree-Fock spin susceptibility in Bulk Graphene
    complex(8) :: XHFBulk1
    real(8), intent(in) :: Vup, Vdown
    real(8), intent(in) ::  w
    complex(8) :: r1, r2
	  
    r1 = zqags(spin_sus_int3,wf-w,wf)
    r2 = zqagi(spin_sus_int12,eta,1)
    XHFBulk1 = r1+r2
    
      contains
      
      function spin_sus_int12(y)
	complex(8) :: spin_sus_int12
	real(8), intent(in) :: y
	spin_sus_int12 = hbar/(2*pi)*( GF(wf + im*y,Vup)*GF(wf + w + im*y,Vdown) + GF(wf - im*y,Vdown)*GF(wf-w- im*y,Vup) )
      end function spin_sus_int12
      
      function spin_sus_int3(w_dum)
	complex(8) :: spin_sus_int3
	real(8), intent(in) :: w_dum
	spin_sus_int3 = -im*hbar/(2*pi)*GF(w_dum - im*eta,Vup)*GF(w + w_dum + im*eta,Vdown)
      end function spin_sus_int3
      
      function GF(E,V)
	complex(8) :: GF
	complex(8), intent(in) :: E
	real(8), intent(in) :: V
	complex(8) :: g
	g = gBulk_kZ(0,0,0,E)
	GF = Dyson1(g,V)
      end function GF
    
  end function XHFBulk1
  
  
  subroutine SCBulkSubs1()
  use shared_data
  use Routines
  use GFMod
    real(8) :: mag_m, tolerance, delta
    real(8) :: Vup, Vdown
    real(8) :: n0

    mag_m = 0.8d0
    tolerance = dtol
    delta = 0
    
    print *, mag_m, tolerance, delta
    
    contains
    
    function FZero()
      real(8) :: FZero
      ! Dummy arguments
      real(8) :: rup, rdown
      real(8) :: nup, ndown
      real(8) :: abserr, epsabs, epsrel
      integer :: ier, last, neval
      integer, parameter :: limit = 100
      integer, parameter :: lenw = limit*4
      integer, dimension(limit) :: iwork
      real(8), dimension(lenw) :: work

      epsabs = 0.0d0
      epsrel = tolerance
      
      call dqagi (GFup, eta, 1, epsabs, epsrel, rup, abserr, neval, ier, limit, lenw, last, iwork, work)
      call dqagi (GFup, eta, 1, epsabs, epsrel, rdown, abserr, neval, ier, limit, lenw, last, iwork, work)
      
      nup =  1.0d0/2.0d0 + rup/pi
      ndown =  1.0d0/2.0d0 + rdown/pi
      
      FZero = n0 - nup - ndown
      
      end function FZero
      
      function GFup(y)
	complex(8) :: GFup
	real(8), intent(in) :: y
	complex(8) :: g
	g = gBulk_kZ(0,0,0,EF+im*y)
	GFup =  real(Dyson1(g,Vup))
      end function GFup
      
      function GFdown(y)
	complex(8) :: GFdown
	real(8), intent(in) :: y
	complex(8) :: g
	g = gBulk_kZ(0,0,0,EF+im*y)
	GFdown =  real(Dyson1(g,Vdown))
      end function GFdown

  end subroutine SCBulkSubs1

end module Dynamic


program main
use Dynamic
use shared_data
implicit none

  call SCBulkSubs1()
  
end program main


! def SC1(GF,n0=1.0):
!   """Calculate Vup/Vdown for a single impurity in graphene."""
!   # This returns values for the up/down spin that are only separated by a sign. There is a symmetry here that you are not exploiting.
!   mag_m = 0.8
!   tolerance = dtol
!   delta = 0.0
!   
!   def n_occ(V):
!     integral = quad(GF, eta, np.inf, args=V, epsabs=0.0, epsrel=dtol, limit=200 )
!     return 1.0/2.0 + integral[0]/pi
!     
!   def FZero(delta):
!     ex_split = U*mag_m
!     Vdown = delta + (ex_split + hw0)/2.0
!     Vup = delta - (ex_split + hw0)/2.0
!     return n0 - n_occ(Vup) - n_occ(Vdown)
! 
!   while True:
!     mag_temp = mag_m
!     delta = newton(FZero, delta, tol=dtol, maxiter=50)
!     ex_split = U*mag_m
!     Vdown = delta + (ex_split + hw0)/2.0
!     Vup = delta - (ex_split + hw0)/2.0
!     mag_m = n_occ(Vup) - n_occ(Vdown)
!     if abs(mag_m - mag_temp) <= tolerance:
!       break
! 
!   return Vup, Vdown

! def SCBulkSubs1(n0=1.0):
!   def GF(y,V):
!     g = gBulk_kZ(0,0,0,EF+1j*y)
!     return Dyson1(g,V).real
!   return SC1(GF,n0)