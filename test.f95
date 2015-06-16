program main
use shared_data
use GFMod
use Routines
use NumRoutines
implicit none
  complex(8) :: E
  complex(8), dimension(2,2) :: g
  integer :: m,n,s
  
  data m,n,s /6,1,0/

  E = (-3.0d0,eta)
  print *, gBulk_kZ(m,n,s,E)


end program main