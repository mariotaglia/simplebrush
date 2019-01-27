
subroutine kai


!#####################################################################
!
! This program calculates the poor-sv coefficient to use in planar 1D
!
! (no PBC in z direction)
!
!
!
!#####################################################################

use const
use kai
use system
use chainsdat

implicit none
integer seed
real*8 xmin,xmax,ymin,ymax,zmin,zmax
integer MCsteps ! number of integration points per lattice site
real*8 R,theta,z
real*8 rn
integer i, ii, is, js, a, b
real*8 rands
real*8 x1,x2,y1, y2, z1, z2, vect
integer iR, ix,iy,iz, itheta
integer j
real*8 radio
real*8 cutoff

print*,'Kai calculation'

cutoff = (float(Xulimit)+0.5)*delta ! limits of integration box
Xu = 0.0 ! vector Xu

seed = 1010
MCsteps = 60*Xulimit

! Defines the integration box

ymax = cutoff
ymin = -cutoff

xmax = cutoff
xmin = -cutoff

zmax = cutoff + 0.5*delta    ! segment is located at the center of a lattice site
zmin = -cutoff + 0.5*delta
      
do ix = 1, MCsteps
do iy = 1, MCsteps
do iz = 1, MCsteps

! segment is at (x1,y1,z1) and point to integrate at (x2,y2,z2)

     z1 = 0.5*delta           ! segment is located at the center of a lattice site
     y1 = 0.0
     z1 = 0.0

     x2 = xmin + (xmax-xmin)*dfloat(ix-1)/dfloat(MCsteps-1)
     y2 = ymin + (ymax-ymin)*dfloat(iy-1)/dfloat(MCsteps-1)
     z2 = zmin + (zmax-zmin)*dfloat(iz-1)/dfloat(MCsteps-1)

     vect = sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2) ! vector 
     j = int(anint(z2/delta))                          ! j has the cell with the point to integrate
      
     if(j.le.dimz) then
         if(vect.le.cutoff) then ! inside cutoff sphere
         if(vect.ge.lseg) then   ! outside segment sphere
              Xu(j) = Xu(j) + ((lseg/vect)**6.0) 
         endif
         endif
     endif
         
enddo!iz
enddo!iy
enddo!ix
      
Xu = Xu /(MCsteps**3)*(2.0*cutoff)**3 ! divide by number of points and multiply by box volume

end
