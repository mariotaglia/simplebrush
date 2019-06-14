subroutine fe(cc, ccc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! this routine calculates the free energy of the system and chemical potential 
! of chains
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use results
use chainsdat
use system
use kai
use bulk
use molecules
implicit none

integer cc, ccc

real*8 Free_energy, F_Mix_s, F_Mix_pos
real*8 F_Mix_neg, F_Mix_Hplus
real*8 Free_energy2, sumpi, sumrho, sumel, sum, mupolA, mupolB, pilat, sumas,diffener,sumex
real*8 F_Mix_OHmin, F_Conf, F_Eq, F_vdW, F_eps, F_electro
integer i, iz, iiz
real*8 xtotal(-Xulimit:dimz+Xulimit) ! xtotal for poor solvent


! Calculation of xtotal

do iz=1,dimz
  xtotal(iz) = 1.0 - xpos(iz) - xneg(iz) - xh(iz) - xHplus(iz) - xOHmin(iz) ! xtotal is everything but solvent and ions
enddo

xtotal(dimz+1:dimz+Xulimit) = 0.0 ! xtotal in bulk = 0.0
xtotal(-Xulimit:0) = 0.0 ! xtotal in surface = 0.0



Free_Energy = 0.0
Free_Energy2 = 0.0

! 1. solvent entropy

F_Mix_s = 0.0 

do iz = 1, dimz
  F_Mix_s = F_Mix_s + xh(iz)*(dlog(xh(iz))-1.0)
  F_Mix_s = F_Mix_s - xsolbulk*(dlog(xsolbulk)-1.0)
enddo      

F_Mix_s = F_Mix_s * delta/vsol
Free_Energy = Free_Energy + F_Mix_s

! 2. cations entropy

F_Mix_pos = 0.0 

do iz = 1, dimz
  F_Mix_pos = F_Mix_pos + xpos(iz)*(dlog(xpos(iz)/vsalt)-1.0 - dlog(expmupos) + dlog(vsalt))
  F_Mix_pos = F_Mix_pos - xposbulk*(dlog(xposbulk/vsalt)-1.0 - dlog(expmupos) + dlog(vsalt))
enddo

F_Mix_pos = F_Mix_pos * delta/vsol/vsalt
Free_Energy = Free_Energy + F_Mix_pos

! 3. anions entropy

F_Mix_neg = 0.0

do iz = 1, dimz
  F_Mix_neg = F_Mix_neg + xneg(iz)*(dlog(xneg(iz)/vsalt)-1.0 - dlog(expmuneg) + dlog(vsalt))
  F_Mix_neg = F_Mix_neg - xnegbulk*(dlog(xnegbulk/vsalt)-1.0 - dlog(expmuneg) + dlog(vsalt))
enddo 

F_Mix_neg = F_Mix_neg * delta/vsol/vsalt
Free_Energy = Free_Energy + F_Mix_neg

! 4. proton entropy

F_Mix_Hplus = 0.0

do iz = 1, dimz
  F_Mix_Hplus = F_Mix_Hplus + xHplus(iz)*(dlog(xHplus(iz))-1.0 - dlog(expmuHplus))
  F_Mix_Hplus = F_Mix_Hplus - xHplusbulk*(dlog(xHplusbulk)-1.0 - dlog(expmuHplus))
enddo
      
F_Mix_Hplus = F_Mix_Hplus * delta/vsol
Free_Energy = Free_Energy + F_Mix_Hplus

! 5. hydroxyl ions entropy

F_Mix_OHmin = 0.0

do iz = 1, dimz
  F_Mix_OHmin = F_Mix_OHmin + xOHmin(iz)*(dlog(xOHmin(iz))-1.0 - dlog(expmuOHmin))
  F_Mix_OHmin = F_Mix_OHmin - xOHminbulk*(dlog(xOHminbulk)-1.0 - dlog(expmuOHmin))
enddo
      
F_Mix_OHmin = F_Mix_OHmin * delta/vsol
Free_Energy = Free_Energy + F_Mix_OHmin

! 6. Conformational entropy of polymer 

F_Conf = 0.0

do i = 1, newcuantas
  F_Conf = F_Conf + (proA(i)/qA)*dlog((proA(i))/qA)*sigmaA
enddo

do i = 1, newcuantas
  F_Conf = F_Conf + (proB(i)/qB)*dlog((proB(i))/qB)*sigmaB
enddo
Free_Energy = Free_Energy + F_Conf

! 7. Chemical Equilibrium
															!!G_ cambie todo esto va a estar mal pero me crasheaba sino !!
F_Eq = 0.0 
	
do iz  = 1, dimz
  F_Eq = F_Eq + (1.0-fdisANC(iz)-fdisAas(iz)-fdisANa(iz))*dlog(1.0-fdisANC(iz)-fdisAas(iz)-fdisANa(iz))*avpolA(iz)/vpol
  F_Eq = F_Eq + (fdisANC(iz))*dlog(fdisANC(iz))*avpolA(iz)/vpol
  F_Eq = F_Eq + (fdisANa(iz))*dlog(fdisANa(iz))*avpolA(iz)/vpol

  F_Eq = F_Eq + (fdisANC(iz))*dlog(K0A)*avpolA(iz)/vpol
  F_Eq = F_Eq + (fdisANC(iz))*(-dlog(expmuHplus))*avpolA(iz)/vpol
  F_Eq = F_Eq + (fdisANa(iz))*(dlog(K0ANa))*avpolA(iz)/vpol
  F_Eq = F_Eq + (fdisANa(iz))*(-dlog(expmupos))*avpolA(iz)/vpol

  F_Eq = F_Eq + (fdisAas(iz))*(-dlog(K0Eo))*avpolA(iz)/vpol
  
	if (1.0d-10 < fdisAas(iz)) then 
 	 F_Eq = F_Eq + (fdisAas(iz))*dlog(fdisAas(iz))*avpolA(iz)/vpol
 		if (1.0d-10 < avpolA(iz))then 
  			F_Eq = F_Eq + (-fdisAas(iz))*(dlog(avpolA(iz)*fdisAas(iz))-1.0)*avpolA(iz)/vpol ! usando que Vpol =Vab
		endif
	endif

enddo
do iz  = 1, dimz
  F_Eq = F_Eq + (1.0-fdisBNC(iz)-fdisBas(iz)-fdisBCl(iz))*dlog(1.0-fdisBNC(iz)-fdisBas(iz)-fdisBCl(iz))*avpolB(iz)/vpol


  F_Eq = F_Eq + (fdisBNC(iz))*dlog(fdisBNC(iz))*avpolB(iz)/vpol
  F_Eq = F_Eq + (fdisBCl(iz))*dlog(fdisBCl(iz))*avpolB(iz)/vpol

  F_Eq = F_Eq + (fdisBNC(iz))*dlog(K0B)*avpolB(iz)/vpol
  F_Eq = F_Eq + (fdisBNC(iz))*(-dlog(expmuOHmin))*avpolB(iz)/vpol
  F_Eq = F_Eq + (fdisBCl(iz))*(dlog(K0BCl))*avpolB(iz)/vpol
  F_Eq = F_Eq + (fdisBCl(iz))*(-dlog(expmuneg))*avpolB(iz)/vpol

  if ((1.0d-10 < fdisBas(iz))) then 
	  F_Eq = F_Eq +( (fdisBas(iz))*dlog(fdisBas(iz)) )*avpolB(iz)/vpol
	endif
enddo

F_eq = F_eq *delta/vsol
Free_Energy = Free_Energy + F_Eq
															!!G_ cambie todo esto va a estar mal pero me crasheaba sino !!
! 8.vdW 

F_vdW = 0.0

do iz = 1, dimz
 do iiz = -Xulimit, Xulimit
   F_vdW = F_vdW - 0.5*delta*Xu(iiZ)*xtotal(iz)*xtotal(iz+iiZ)*st/(vpol*vpol*vsol*vsol)
 enddo
enddo

Free_Energy = Free_Energy + F_vdW

! 9. Electrostatic ! VER ESTO...

F_electro = 0.0    

do iz  = 1, dimz
  F_electro = F_electro + delta*psi(iz)*qtot(iz)/2.0/vsol
enddo

Free_Energy = Free_Energy + F_electro

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Free energy at extrema 

Free_Energy2 = 0.0

sumpi = 0.0
sumrho=0.0
sumel=0.0
sumas=0.0

do i=1,dimz
            
   sumpi = sumpi+dlog(xh(i))     
   sumpi = sumpi-dlog(xsolbulk)     
   sumrho = sumrho + ( - xh(i) -xHplus(i) -xOHmin(i)-(xpos(i)+xneg(i))/vsalt)! sum over  rho_i i=+,-,si
   sumrho = sumrho - ( - xsolbulk -xHplusbulk -xOHminbulk-(xposbulk+xnegbulk)/vsalt)! sum over  rho_i i=+,-,si
   sumel = sumel - qtot(i)*psi(i)/2.0 ! electrostatic part free energy	
	sumas = sumas +  avpolA(i)*fdisAas(i)/vpol


enddo         
        
sumpi = (delta/vsol)*sumpi
sumrho = (delta/vsol)*sumrho
sumel = (delta/vsol)*sumel
sumas = (delta/vsol)*sumas

sum = sumpi + sumrho + sumel +sumas

Free_Energy2 = -sigmaA*dlog(qA)-sigmaB*dlog(qB) + sum -F_vdW 

! Chemical potential chains
mupolA = -dlog(qA)
mupolB = -dlog(qB)
 
! Pilat
pilat = sigmaA*mupolA+sigmaB*mupolB - Free_energy

! Save to disk

write(301,*)cc, ccc, Free_energy
write(302,*)cc, ccc, F_Mix_s 
write(303,*)cc, ccc, F_Mix_pos
write(304,*)cc, ccc, F_Mix_neg
write(305,*)cc, ccc, F_Mix_Hplus
write(306,*)cc, ccc, F_Mix_OHmin
write(307,*)cc, ccc, F_Conf
write(308,*)cc, ccc, F_Eq
write(309,*)cc, ccc, F_vdW
write(311,*)cc, ccc, F_electro
write(312,*)cc, ccc, Free_energy2
write(313,*)cc, ccc, mupolA
write(314,*)cc, ccc, pilat
write(315,*)cc, ccc, mupolB
! print
diffener= Free_energy- Free_energy2

print*, 'Free energy:', Free_energy, Free_energy2,diffener

end subroutine
 
