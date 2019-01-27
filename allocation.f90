subroutine allocation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine allocates memory space using read-dependent variables
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use kai
use system
use chainsdat
use solver
use results
use mkinsol
implicit none

allocate(Xu(-Xulimit:Xulimit))
allocate(xflag(2*dimz))
allocate(xh(dimz))
allocate(psi(dimz))
allocate(xpos(dimz))
allocate(xneg(dimz))
allocate(xHplus(dimz))
allocate(xOHmin(dimz))
allocate(qtot(dimz))
allocate(pro(cuantas))
allocate(avpol(dimz))
allocate(fdis(dimz))
allocate(pz(cuantas,long))
allocate(pp(2*dimz))
end subroutine
