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

allocate(Xu(-Xulimit:Xulimit)) !vdw
allocate(xflag(3*dimz))
allocate(xh(dimz)) !xsolvent
allocate(psi(dimz)) 
allocate(xpos(dimz))
allocate(xneg(dimz))
allocate(xHplus(dimz))
allocate(xOHmin(dimz))
allocate(qtot(dimz))
allocate(proA(cuantas)) ! It could be equal to both polyms, but may be depends if we took or not pkaA=!pkaB ;proA=proB=pro
allocate(proB(cuantas)) ! It could be equal to both polyms, but may be depends if we took or not pkaA=!pkaB ;proA=proB=pro
allocate(avpolA(dimz)) !	
allocate(avpolB(dimz)) !
allocate(xna(dimz))
allocate(xnb(dimz))
allocate(eta(dimz))
allocate(m(dimz))	
allocate(fdisANC(dimz)) !	fraction not charge pol-A
allocate(fdisBNC(dimz)) !	fraction not charge pol-B
allocate(fdisAas(dimz)) !	fraction associate pol-A
allocate(fdisBas(dimz)) !	fraction associate pol-B
allocate(pzA(cuantas,long)) !This will be related with the position of the segment, of the polyms, so it be different for each one. 
allocate(pzB(cuantas,long)) !
allocate(pp(3*dimz))
end subroutine
