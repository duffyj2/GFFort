! The GFs in FORTRAN. Hopefully all modern and stuff.

module shared_data
implicit none
! Contains various parameters used throughout the work
  save
  ! Mathematical parameters
  complex(8), parameter :: im = (0.0d0,1.0d0)
  real(8), parameter :: pi = acos(-1.0d0)
  real(8), parameter :: hbar = 1.0
  ! Material parameters
  real(8), parameter :: t = -1.0d0
  real(8), parameter :: eps_imp = 0.0
  real(8), parameter :: tau = -1.0
  ! Accuracy parameters
  real(8), parameter :: eta = 1.0d-4
  real(8), parameter :: dtol = 1.0d-4
end module


module GFMod
implicit none

  contains 

  function gBulk_kZ(m,n,s,E)
  use shared_data
  use NumRoutines
    complex(8) :: gBulk_kZ	! Should this be an I/O variable?
    ! Input variables
    complex(8), intent(in) :: E
    integer, intent(in) :: m, n, s
    ! Dummy variables
    real(8) :: lim1, lim2
    
    lim1 = -pi/2.0d0
    lim2 = pi/2.0d0
    
    gBulk_kZ = complexIntegrate(gbulk_kz_int,lim1,lim2)
    
    contains

    function gbulk_kz_int(kZ)
      ! Input arguments
      complex(8) :: gbulk_kz_int
      real(8), intent(in) :: kZ
      ! Dummy arguments
      complex(8) :: f, ft 
      complex(8) :: q
      complex(8) :: Const, Den
      integer :: sig

      q = acos( (E**2 - t**2 - 4.0d0*t**2 *cos(kZ)**2)/(4.0d0*t**2 *cos(kZ) ) )

      if (aimag(q) < 0.0d0) q = -q

      Const = im/(4*pi*t**2)
      Den = cos(kZ)*sin(q)

      if (s == 0) then
	sig = sign(1,m+n)
	gbulk_kz_int = Const*E*exp( im*(sig*q*(m+n) + kZ*(m-n) ) )/ Den 
      else if (s == 1) then 
	sig = sign(1,m+n)
	f = t*( 1.0d0 + 2.0d0*cos(kZ)*exp(sig*im*q) )
	gbulk_kz_int = Const*f*exp( im*(sig*q*(m+n) + kZ*(m-n) ) )/Den  
      else if (s == -1) then
	sig = sign(1,m+n-1)
	ft = t*( 1.0d0 + 2.0d0*cos(kZ)*exp(-sig*im*q) )
	gbulk_kz_int = Const*ft*exp( im*(sig*q*(m+n) + kZ*(m-n) ) )/Den 
      else
	print *, "Sublattice error in gbulk_kz_int"
      end if
	
    end function gbulk_kz_int 
    
  end function gBulk_kZ
  

  function gRib_Arm(nE,m1,n1,m2,n2,s,E)
  use shared_data
    complex(8) :: gRib_Arm
    integer :: j
    ! Other parameters
    complex(8) :: E
    integer :: m1, n1, m2, n2, s
    integer :: nE

    gRib_Arm = 0.0d0
    if ( mod(nE,2) .eq. 0 ) then
      do j = 1, nE-1
	if ( j .eq. nE/2 ) cycle		! Avoid singularities
	gRib_Arm = gRib_Arm + g_term(pi*j/nE)
      end do
      if ( m2+n2-m1-n1 .eq. 0 ) then
	gRib_Arm = gRib_Arm + limit_term(pi/2)
      end if
    else 
      do j = 1, nE-1
	gRib_Arm = gRib_Arm + g_term(pi*j/nE)
      end do
    end if

    contains

      function g_term (ky)
	complex(8) :: g_term
	real(8), intent(in) :: ky

	complex(8) :: f, ft 
	complex(8) :: q
	complex(8) :: Const, Den
	integer :: sig

	q = acos( (E**2 - t**2 - 4.0d0*t**2 *cos(ky)**2)/(4.0d0*t**2 *cos(ky) ) )
	if (aimag(q) < 0.0d0) q = -q

	Const = im/(2.0d0*nE*t**2)
	Den = cos(ky)*sin(q)

	if (s == 0) then
	  sig = sign(1,m2+n2-m1-n1)
	  g_term = Const*E*exp( im*sig*q*(m2+n2-m1-n1) )*sin(ky*(m2-n2))*sin(ky*(m1-n1))/ Den 
	else if (s == 1) then 
	  sig = sign(1,m2+n2-m1-n1)
	  f = 1.0d0 + 2.0d0*cos(ky)*exp(sig*im*q)
	  g_term = Const*t*f*exp( im*sig*q*(m2+n2-m1-n1) )*sin(ky*(m2-n2))*sin(ky*(m1-n1))/ Den 
	else if (s == -1) then
	  sig = sign(1,m2+n2-m1-n1-1)
	  ft = 1.0d0 + 2.0d0*cos(ky)*exp(-sig*im*q)
	  g_term = Const*t*ft*exp( im*sig*q*(m2+n2-m1-n1) )*sin(ky*(m2-n2))*sin(ky*(m1-n1))/ Den 
	else
	  write(*,*) 'Sublattice error in gRib_Arm'
	end if

      end function g_term 
      
      function limit_term (ky)
	complex(8) :: limit_term
	real(8), intent(in) :: ky
	complex(8) :: N_ab

	if (s == 0) then
	  N_ab = E
	else if ( (s == 1) .or. (s == -1) ) then
	  N_ab = t
	else
	  write(*,*) 'Sublattice error in gRib_Arm'
	end if
	
	limit_term = 2.0d0*N_ab*sin(ky*(m2-n2))*sin(ky*(m1-n1))/( nE*( E**2 - t**2 ) )

      end function limit_term 

  end function gRib_Arm


  function gtube_arm(nC,m,n,s,E)
  use shared_data
    ! Calculates the GF for an armchair nanotube
    ! Input arguments
    complex(8) :: gtube_arm
    complex(8), intent(in) :: E
    integer, intent(in) :: m,n,s
    integer, intent(in) :: nC
    ! Dummy arguments
    complex(8) :: qp, qm
    complex(8) :: const
    complex(8) :: fm, fp, ftm, ftp
    integer :: k
    integer :: sig
    
    gtube_arm = 0.0d0    
    do k = 0, nC-1	! possible

      qp = acos( 0.5d0*(-cos(pi*k/nC) + sqrt( (E**2/t**2) - (sin(pi*k/nC))**2 ) )  )
      qm = acos( 0.5d0*(-cos(pi*k/nC) - sqrt( (E**2/t**2) - (sin(pi*k/nC))**2 ) )  ) 

      if (aimag(qp) < 0.0) qp = -qp
      if (aimag(qm) < 0.0) qm = -qm

      sig = sign(1,m-n)		! Really fucking hope this is correct. Doesn't really match up with the other ones
      const = im/(4.0d0*t**2)
      if (s == 0) then
	gtube_arm = gtube_arm + const*( E*exp( im*( pi*k/nC *(m+n) + sig*qp*(m-n) )  ) / ( sin(2.0d0*qp) + sin(qp)*cos(pi*k/nC)  ) &
	  + E*exp( im*(  pi*k/nC *(m+n) + sig*qm*(m-n) )  ) / ( sin(2.0d0*qm) + sin(qm)*cos(pi*k/nC)  )  ) 
      else if (s == 1) then
	fp = t*( 1.0d0 + 2.0d0*cos(qp)*exp(im*pi*k/nC)  )
	fm = t*( 1.0d0 + 2.0d0*cos(qm)*exp(im*pi*k/nC)  )

	gtube_arm = gtube_arm + const*( fp*exp( im*( pi*k/nC *(m+n) + sig*qp*(m-n) )  ) / ( sin(2.0d0*qp) + sin(qp)*cos(pi*k/nC)  ) &
	+ fm*exp( im*( pi*k/nC *(m+n) + sig*qm*(m-n) )  ) / ( sin(2.0d0*qm) + sin(qm)*cos(pi*k/nC)  )  ) 
      else if (s == -1) then
	ftp = t*( 1.0d0 + 2.0d0*cos(qp)*exp(-im*pi*k/nC)  )
	ftm = t*( 1.0d0 + 2.0d0*cos(qm)*exp(-im*pi*k/nC)  )

	gtube_arm = gtube_arm + const*( exp( im*( pi*k/nC *(m+n) + sig*qp*(m-n) )  )*ftp / ( sin(2.0d0*qp) + sin(qp)*cos(pi*k/nC)  ) &
	+ exp( im*( pi*k/nC *(m+n) + sig*qm*(m-n) )  )*ftm / ( sin(2.0d0*qm) + sin(qm)*cos(pi*k/nC)  )  ) 
      else
	print *, "Sublattice error in gTube_Arm"
      end if
    end do
    
    gtube_arm = gtube_arm/nC
      
  end function gtube_arm


  function gsi_kz_int(m1,n1,m2,n2,s,E,kZ)
  use shared_data
    !The Semi-Infinite Graphene Green's function
    !The kZ integration is performed last
    ! Input arguments
    complex(8) :: gsi_kz_int
    complex(8), intent(in) :: E
    real(8), intent(in) :: kZ
    integer, intent(in) :: m1,n1,m2,n2,s
    ! Dummy arguments
    complex(8) :: q, f, ft
    complex(8) :: const, den
    integer :: sig
    
    q = acos( (E**2 - t**2 - 4.0d0*t**2 *cos(kZ)**2)/(4.0d0*t**2 *cos(kZ) ) )
    if (aimag(q) < 0.0d0) q = -q
    
    const = im/(2.0d0*pi*t**2)
    den = cos(kZ)*sin(q)
    
    if (s == 0) then
      sig = sign(1,m2+n2-m1-n1)
      gsi_kz_int = const*E*exp( im*sig*q*(m2+n2-m1-n1) )*sin( kZ*(m1-n1) )*sin( kZ*(m2-n2) )/den
    else if (s == 1) then
      sig = sign(1,m2+n2-m1-n1)
      f = t*( 1.0d0 + 2.0d0*cos(kZ)*exp(im*sig*q) )
      gsi_kz_int = const*f*exp( im*sig*q*(m2+n2-m1-n1) )*sin( kZ*(m1-n1) )*sin( kZ*(m2-n2) )/den
    else if (s == -1) then
      sig = sign(1,m2+n2-m1-n1-1)
      ft = t*( 1.0d0 + 2.0d0*cos(kZ)*exp(-im*sig*q) )
      gsi_kz_int = const*ft*exp( im*sig*q*(m2+n2-m1-n1) )*sin( kZ*(m1-n1) )*sin( kZ*(m2-n2) )/den
    else
      print *, "Sublattice error in gSI_kZ"
    end if

  end function gsi_kz_int

end module GFMod



module NumRoutines
implicit none
! Currently just a module for integration

  contains

  function complexIntegrate(f,lim1,lim2)
  use shared_data
  ! Integrates a complex function, uses the default tolerance
  ! Probably should include an optional parameter for the tolerance
  ! Also, doesn't return the error
  ! Might be better as a subroutine
    complex(8) :: complexIntegrate
    complex(8), external :: f
    real(8), intent(in) :: lim1,lim2
    ! Dummy
    real(8) :: rRe, rIm
    real(8) :: epsabs, epsrel, abserr
    integer :: ier, last, neval
    integer, parameter :: limit = 500
    integer, parameter :: lenw = limit*4
    integer, dimension(limit) :: iwork
    real(8), dimension(lenw) :: work

    epsabs = 0.0d0
    epsrel = dtol
    
    call dqags (fRe, lim1, lim2, epsabs, epsrel, rRe, abserr, neval, ier, limit, lenw, last, iwork, work)
    call dqags (fIm, lim1, lim2, epsabs, epsrel, rIm, abserr, neval, ier, limit, lenw, last, iwork, work)
        
    complexIntegrate = cmplx(rRe,rIm,8)
    
    contains
    
    function fRe(x)
      real(8) :: fRe
      real(8), intent(in) ::  x
      fRe = real( f(x) )
    end function fRe
    
    function fIm(x)
      real(8) :: fIm
      real(8), intent(in) ::  x
      fIm = aimag( f(x) )
    end function fIm
    
  end function
  
  
  function complexIntegrateInf(f,bound,inf)
  use shared_data
  ! Integrates a complex function over infinite ranges, uses the default tolerance
  ! Probably should include an optional parameter for the tolerance
  ! Also, doesn't return the error
    complex(8) :: complexIntegrateInf
    complex(8), external :: f
    real(8), intent(in) :: bound
    integer, intent(in) :: inf
    ! Dummy
    real(8) :: rRe, rIm
    real(8) :: abserr, epsabs, epsrel
    integer :: ier, last, neval
    integer, parameter :: limit = 100
    integer, parameter :: lenw = limit*4
    integer, dimension(limit) :: iwork
    real(8), dimension(lenw) :: work

    epsabs = 0.0d0
    epsrel = dtol
    
    call dqagi (fRe, bound, inf, epsabs, epsrel, rRe, abserr, neval, ier, limit, lenw, last, iwork, work)
    call dqagi (fRe, bound, inf, epsabs, epsrel, rIm, abserr, neval, ier, limit, lenw, last, iwork, work)
        
    complexIntegrateInf = cmplx(rRe,rIm,8)
    
    contains
    
    function fRe(x)
      real(8) :: fRe
      real(8), intent(in) ::  x
      fRe = real( f(x) )
    end function fRe
    
    function fIm(x)
      real(8) :: fIm
      real(8), intent(in) ::  x
      fIm = aimag( f(x) )
    end function fIm
    
  end function
  
  
  function inv(m)
  use, intrinsic :: iso_fortran_env, only : error_unit
    complex(8), dimension(:,:) :: m
    complex(8), dimension(size(m,1),size(m,1)) :: inv
    complex(8),dimension(2*size(m,1),2*size(m,1)) :: q
    complex(8), dimension(64*size(m,1)) :: work
    integer, dimension(2*size(m,1)) :: ipiv
    integer :: info, lwork, lda, i, j, n
    n = size(m,1)
    lwork = 64*n ! blocksize * n, where blocksize is a machine dependent variable 16 
    lda   = 2*n

    q = 0
    q(1:n,1:n) = m

    call zgetrf(n,n,q,lda,ipiv,info)
    if (info.eq.0) then
    call zgetri(n,q,lda,ipiv,work,lwork,info)
    else
      ! I'd prefer to write this to the error output
      write (error_unit,*) 'Matrix is numerically singular!'
    end if

    inv = q(1:n,1:n)

  end function
  
  function Imx(n)
  implicit none
    integer, intent(in) :: n
    real(8), dimension(n,n) :: Imx
    integer :: i
    
    Imx = 0.0d0 
    do i = 1, n
      Imx(i,i) = 1
    end do
    
  end function Imx
    
end module 



