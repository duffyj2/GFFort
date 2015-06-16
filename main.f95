module Routines
implicit none

  contains
  
  function Dyson1(g,V)
    complex(8) :: Dyson1
    complex(8), intent(in) :: g
    real(8), intent(in) :: V
    
    Dyson1 = g/(1.0-g*V)
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
    g_imp = 1.0/(E-eps_imp)
    g(3,3) = g_imp
    g(4,4) = g(3,3)
    
    ! The peturbation connects the impurities to the lattice
    data V(3,1),V(1,3),V(2,4),V(4,2) /4*tau/
    
    gtemp = Dyson(g,V)
    gGNRTopMx = gtemp(3:4,3:4)
  end function gGNRTopMx

  
end module Routines


program main
use shared_data
use GFMod
use Routines
use NumRoutines
implicit none
  complex(8) :: E
  integer :: m,n,s
  !
  real(8) :: Vup, Vdown
  real(8) :: wf, w, w_dum
  real(8) :: y
  
  m = 3
  n = 3
  s = 1
 
  E = (1.2d0,eta)


  Vup = -1.2d0
  Vdown = -Vup
  wf = 0.0d0
  w = 1.0d0
  y = 1.0d0
  
!   print *, spin_sus_int12(w)  
!   spin_sus_int3,wf-w,wf
  print *, complexIntegrate(spin_sus_int3,wf-w,wf)
  
  contains
  
    function spin_sus_int12(w)
      complex(8) :: spin_sus_int12
      real(8), intent(in) :: w
      spin_sus_int12 = hbar/(2.0*pi)*( GF(wf + im*y,Vup)*GF(wf + w + im*y,Vdown) + GF(wf - im*y,Vdown)*GF(wf-w- im*y,Vup) )
    end function spin_sus_int12
    
    function spin_sus_int3(w_dum)
      complex(8) :: spin_sus_int3
      real(8), intent(in) :: w_dum
      spin_sus_int3 = -im*hbar/(2.0*pi)*GF(w_dum - im*eta,Vup)*GF(w + w_dum + im*eta,Vdown)
    end function spin_sus_int3
    
    function GF(E,V)
      complex(8) :: GF
      complex(8), intent(in) :: E
      real(8), intent(in) :: V
      complex(8) :: g
      g = gBulk_kZ(0,0,0,E)
      GF = Dyson1(g,V)
    end function GF


end program main