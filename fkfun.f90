subroutine fkfun(x,f,ier)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! User provided routine for kinsol
! x is the input vector
! f is the output vector, kinsol will change x in order to get f = 0
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use results
use chainsdat
use solver
use system
use molecules
use bulk
use const
use kai
implicit none
 
integer ntot
real*8 x(2*dimz),f(2*dimz)
real*8 protemp, protemp1
integer i,j, k, ix, iy, iz, ii, ax, ay, az, temp, iiZ
real*8 temp2


real*8 xpot(dimz)
real*8 psi2(0:dimz+1) ! psi plus boundaries at z=0 and dimz+1
real*8 xtotal(-Xulimit:dimz+Xulimit) ! xtotal for poor solvent
 
! Kinsol
integer*4 ier

! Recovers xh and psi from x
ntot = dimz 

do iz=1,dimz
   xh(iz)=x(iz)
   psi(iz)=x(iz+ntot)
enddo
psi2(1:dimz) = psi(1:dimZ)

! Electrostatic potential boundary conditions
psi2(dimz+1) = 0.0 ! bulk
psi2(0) = psi2(1)  ! wall, no charge

! Volume fractions of ions and fdis

do iz=1,dimz

   xpos(iz) = expmupos*(xh(iz)**vsalt)*dexp(-psi2(iz)*zpos) ! ion plus volume fraction
   xneg(iz) = expmuneg*(xh(iz)**vsalt)*dexp(-psi2(iz)*zneg) ! ion neg volume fraction
   xHplus(iz) = expmuHplus*(xh(iz))*dexp(-psi2(iz))         ! H+ volume fraction
   xOHmin(iz) = expmuOHmin*(xh(iz))*dexp(+psi2(iz))         ! OH-  volume fraction
   fdis(iz) = 1.0 /(1.0 + xHplus(iz)/(K0*xh(iz)))
   
enddo

! Calculation of xtotal

do iz=1,dimz
  xtotal(iz) = 1.0 - xpos(iz) - xneg(iz) - xh(iz) - xHplus(iz) - xOHmin(iz) ! xtotal is everything but solvent and ions
enddo

xtotal(dimz+1:dimz+Xulimit) = 0.0 ! xtotal in bulk = 0.0
xtotal(-Xulimit:0) = 0.0 ! xtotal in surface = 0.0

! Calculation of xpot

do iz = 1, dimz
  xpot(iz) = xh(iz)**vpol*dexp(-psi2(iz)*zpol)/fdis(iz)

  do iiZ = -Xulimit, Xulimit
    xpot(iz) = xpot(iZ)*dexp(xtotal(iz+iiz)*Xu(iiZ)*st/(vpol*vsol))
  enddo
enddo

! Calculation of pro

q = 0.0
avpol = 0.0
  
do i=1,cuantas

  pro(i)=1.0

  do j=1,long
    az = pz(i, j)         
    pro(i) = pro(i) * xpot(az)
  enddo
  q=q+pro(i)
  do j = 1,long
    az = pz(i, j)
    avpol(az) = avpol(az) + pro(i)*sigma*vsol*vpol/delta
  enddo

enddo ! i

do iz=1, dimz            ! norm avpol
    avpol(iz)=avpol(iz)/q
enddo

!         temp2 = 0.0
!         do iz = 1, dimz
!         temp2 = temp2 + avpol(iz)
!         enddo
!         print*, temp2, sigma, vpol




! Construct f vector

! calculation of qtot (charge in units of |e|/vsol)

do iz=1,dimz
 qtot(iz) = (zpos*xpos(iz)+zneg*xneg(iz))/vsalt + avpol(iz)*zpol/vpol*fdis(iz) + xHplus(iz)-xOHmin(iz) 
enddo

! first block of f,  packing constraint

do iz=1,dimz
 f(iz)=avpol(iz) + xh(iz) + xneg(iz) + xpos(iz) + xHplus(iz) + xOHmin(iz) - 1.000000d0
enddo

! second block of f, Poisson Equation

do iz=1,dimz
  f(iz+ntot)= psi2(iz+1) -2*psi2(iz) + psi2(iz-1) + qtot(iz)*constq
  f(iz+ntot)=f(iz+ntot)/(-2.0)
enddo
 
iter = iter + 1
norma = 0.0

do i = 1, 2*ntot
norma = norma +(f(i))**2    
enddo

print*, iter, norma, q
ier = 0.0
   
return
end subroutine
