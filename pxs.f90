subroutine pxs(in1,il)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Puts chains into lattice
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use system
use chainsdat
implicit none
    
real*8 in1(long,3)  ! Posicion de cada segmento
integer i, il,j, k, ix, iy, iz, ii

do j=1,long
   pz(il,j)=int((in1(j,1))/delta)+1
     if(pz(il, j).gt.dimz) then
       print*,'Increase dimz', il, j, pz(il, j)
       call endall
     endif         
enddo        
return
end subroutine
