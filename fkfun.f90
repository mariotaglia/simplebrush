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
real*8 x(3*dimz),f(3*dimz)
real*8 protemp, protemp1
integer i,j, k, ix, iy, iz, ii, ax, ay, az, temp, iiZ
real*8 temp2
real*8 vpair
real*8 xpotA(dimz)
real*8 xpotB(dimz)
real*8 psi2(0:dimz+1) ! psi plus boundaries at z=0 and dimz+1
real*8 xtotal(-Xulimit:dimz+Xulimit) ! xtotal for poor solvent

 
! Kinsol
integer*4 ier


! Recovers xh and psi from x
ntot = dimz 

do iz=1,dimz
   xh(iz)=x(iz)
   psi(iz)=x(iz+ntot)	
	xna(iz)=x(iz+ntot+ntot)	!!!! <n_a>
enddo
psi2(1:dimz) = psi(1:dimZ)
 	!!!! <n_a>

! Electrostatic potential boundary conditions
psi2(dimz+1) = psi2(dimz) ! wall, no charge
psi2(0) = psi2(1)  ! wall, no charge

! Volume fractions of ions and fdis

do iz=1,dimz

   xpos(iz) = expmupos*(xh(iz)**vsalt)*dexp(-psi2(iz)*zpos) ! ion plus volume fraction
   xneg(iz) = expmuneg*(xh(iz)**vsalt)*dexp(-psi2(iz)*zneg) ! ion neg volume fraction
   xHplus(iz) = expmuHplus*(xh(iz))*dexp(-psi2(iz))         ! H+ volume fraction
   xOHmin(iz) = expmuOHmin*(xh(iz))*dexp(+psi2(iz))         ! OH-  volume fraction
	xnb(iz)=1.0-xna(iz)-xh(iz)-xpos(iz)-xneg(iz)-xHplus(iz)-xOHmin(iz)
   fdisAas(iz)=0.0d0
	fdisBas(iz)=0.0d0
	if ((0.0 .lt. xna(iz)).AND.(0.0 .lt. xnb(iz))) THEN		!g 
		eta(iz)=xna(iz)/xnb(iz)
		M(iz)=( 1.0+ (xOHmin(iz))/(K0B*xh(iz)) )*( 1.0+ (xHplus(iz))/(K0A*xh(iz)) )/(K0Eo*vpol*xna(iz)) 							   !!gabi: vpair = vpol!!
 	 	fdisAas(iz) = -(1.0+eta(iz)+eta(iz)*M(iz))/(2.0*eta(iz))+(	1.0/eta(iz)+(	(1.0+eta(iz)+eta(iz)*M(iz))/(2.0*eta(iz))	)**2	)**0.5 !!!!!!!!!!!!!!!!!!!!!!!
   	fdisBas(iz) = -(1.0+eta(iz)+eta(iz)*M(iz))/(2.0)+eta(iz)*(	1.0/eta(iz)+(	(1.0+eta(iz)+eta(iz)*M(iz))/(2.0*eta(iz))	)**2	)**0.5 !!!!!!!!!!!!!!!!!!!!!!!
		!print*, 'fdisAas(iz) y fdisBas(iz):', fdisAas(iz), fdisBas(iz)
		!fdisAas(iz) =0.0
		!fdisBas(iz) =0.0
	endif
	fdisANC(iz) = (1.0-fdisAas(iz) )/(1.0 + (K0A*xh(iz))/(xHplus(iz)))						   !g
   fdisBNC(iz) = (1.0-fdisBas(iz) )/(1.0 + (K0B*xh(iz))/(xOHmin(iz)))						   !g
  	KK0check(iz)=-dlog10( (Na/1.0d24)*fdisBas(iz)/(	(1.0-fdisAas(iz)-fdisANC(iz))*(1.0-fdisBas(iz)-fdisBNC(iz))*vpol*xna(iz) )	)-pKeo
	KKaAcheckplus(iz)= -dlog10( (xHplus(iz)/xh(iz))*((1-fdisANC(iz)-fdisAas(iz))/fdisANC(iz))*(xsolbulk*1.0d24/(Na*vsol)))-pKaA !! esto era para chequear pkaA
	kkaBcheckmin(iz)=	 (xOhmin(iz)/xh(iz))*(1.0-fdisBas(iz)-fdisBNC(iz))/fdisBNC(iz)-K0B
	print*, 'KKcheck,fdisAas,fdisBas:', KK0check(iz),	KKaAcheckplus(iz),kkaBcheckmin(iz), fdisAas(iz),fdisBas(iz)
enddo

! Calculation of xtotal (KaA*vsol/xsolbulk)*(Na/1.0d24)!

do iz=1,dimz
  xtotal(iz) = 1.0 - xpos(iz) - xneg(iz) - xh(iz) - xHplus(iz) - xOHmin(iz) ! xtotal is everything but solvent and ions
enddo

xtotal(dimz+1:dimz+Xulimit) = 0.0 ! xtotal in bulk = 0.0
xtotal(-Xulimit:0) = 0.0 ! xtotal in surface = 0.0

! Calculation of xpot; es la suma de P 

do iz = 1, dimz
  xpotA(iz) = xh(iz)**vpol*dexp(-psi2(iz)*zpolA)/(1.0-fdisAas(iz)-fdisANC(iz)) 
  xpotB(iz) = xh(iz)**vpol*dexp(-psi2(iz)*zpolB)/(1.0-fdisBas(iz)-fdisBNC(iz))

  do iiZ = -Xulimit, Xulimit
    xpotA(iz) = xpotA(iZ)*dexp(xtotal(iz+iiz)*Xu(iiZ)*st/(vpol*vsol))
    xpotB(iz) = xpotB(iZ)*dexp(xtotal(iz+iiz)*Xu(iiZ)*st/(vpol*vsol))

  enddo
enddo

! Calculation of pro

qA = 0.0
avpolA = 0.0
qB = 0.0
avpolB = 0.0
  
do i=1,newcuantas

  proA(i)=1.0
  proB(i)=1.0

  do j=1,long
    az = pzA(i, j)         
    proA(i) = proA(i) * xpotA(az)
  enddo
  qA=qA+proA(i)
  do j = 1,long
    az = pzA(i, j)
    avpolA(az) = avpolA(az) + proA(i)*sigmaA*vsol*vpol/delta
  enddo
  do j=1,long
    az = pzB(i, j)         
    proB(i) = proB(i) * xpotB(az)
  enddo
  qB=qB+proB(i)
  do j = 1,long
    az = pzB(i, j)
    avpolB(az) = avpolB(az) + proB(i)*sigmaB*vsol*vpol/delta
  enddo

enddo ! i

do iz=1, dimz            ! norm avpol
    avpolA(iz)=avpolA(iz)/qA
enddo
do iz=1, dimz            ! norm avpol
    avpolB(iz)=avpolB(iz)/qB
enddo
!         temp2 = 0.0
!         do iz = 1, dimz
!         temp2 = temp2 + avpol(iz)
!         enddo
!         print*, temp2, sigma, vpol




! Construct f vector

! calculation of qtot (charge in units of |e|/vsol)

do iz=1,dimz
 qtot(iz) = (zpos*xpos(iz)+zneg*xneg(iz))/vsalt + avpolA(iz)*zpolA/vpol*(1.0-fdisANC(iz)-fdisAas(iz))&
+ avpolB(iz)*zpolB/vpol*(1.0-fdisBNC(iz)-fdisBas(iz)) + xHplus(iz)-xOHmin(iz) 
enddo

! first block of f,  packing constraint

do iz=1,dimz
 f(iz)=avpolA(iz) + avpolB(iz) + xh(iz) + xneg(iz) + xpos(iz) + xHplus(iz) + xOHmin(iz) - 1.000000d0
enddo

! second block of f, Poisson Equation

do iz=1,dimz
  f(iz+ntot)= psi2(iz+1) -2*psi2(iz) + psi2(iz-1) + qtot(iz)*constq
  f(iz+ntot)=f(iz+ntot)/(-2.0)
enddo

! third block of f, ASSOCIATE RELATION
 
do iz=1,dimz
	f(iz+ntot+ntot)=avpolA(iz)-xna(iz) !		!
enddo

iter = iter + 1
norma = 0.0

do i = 1, 3*ntot
norma = norma +(f(i))**2    
enddo

print*, iter, norma, qA, qB
ier = 0.0
   
return
end subroutine
