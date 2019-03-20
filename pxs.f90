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
integer flag

flag=0
do j=1,long
   pzA(newcuantas+1,j)=int((in1(j,1))/delta)+1
     if(pzA(newcuantas+1,j).gt.dimz) then
		flag=1
     endif         
	pzB(newcuantas+1,j)=int((float(dimz)*delta-in1(j,1))/delta)+1
     if(pzB(newcuantas+1, j).lt.0) then
       flag=1
     endif  	
enddo       
if(flag.eq.0) then
	newcuantas=newcuantas+1
endif
return
end subroutine
