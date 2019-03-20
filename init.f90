subroutine init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
! This subroutine initializes input-dependent variables
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use molecules
use const
use bulk
use system
implicit none
real*8 KaA
real*8 KaB
real*8 cHplus, cOHmin
real*8 pOHbulk
real*8 xsalt
real*8 Kw

zpos = 1.0      ! charge of cation
zneg = -1.0     ! charge of anion
!zpol = -1.0     ! charge of polyelectrolyte segment 
zpolA = -1.0      ! charge of polyelectrolyte segment A
zpolB = 1.0      ! charge of polyelectrolyte segment B
vsol = 0.030    ! volume of solvent molecule
vsalt=((4.0/3.0)*pi*(0.27)**3)/vsol  ! volume salt in units of vsol 0.27=radius salt  
!vpol= 0.095/vsol!                     ! volume polymer segment in units of vsol 
vpol= 0.095/vsol!                     ! volume polymer segment in units of vsol 

constq=delta*delta*4.0*pi*lb/vsol   ! multiplicative factor in poisson eq  
pKw = 14.0 ! -log10(Kw)
kW=10**(-pKw)
!Esta primer parte esta relacionada con el sistema sin el polimero por lo que quedaría igual
KaA=10**(-pKaA)
KaB=10**(-pKaB)
cHplus = 10**(-pHbulk)    ! concentration H+ in bulk
xHplusbulk = (cHplus*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol
pOHbulk= pKw -pHbulk
cOHmin = 10**(-pOHbulk)   ! concentration OH- in bulk
xOHminbulk = (cOHmin*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol  

xsalt=(csalt*Na/(1.0d24))*(vsalt*vsol)   ! volume fraction salt,csalt in mol/l 

if(pHbulk.le.7) then  ! pH<= 7
   xposbulk=xsalt/zpos
   xnegbulk= -xsalt/zneg +(xHplusbulk -xOHminbulk) *vsalt ! NaCl+ HCl  
else                  ! pH >7 
   xposbulk=xsalt/zpos +(xOHminbulk -xHplusbulk) *vsalt ! NaCl+ NaOH   
   xnegbulk=-xsalt/zneg
endif

xsolbulk=1.0 -xHplusbulk -xOHminbulk - xnegbulk -xposbulk 

K0A = (KaA*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant 
K0B = (Kw/KaB*vsol/xsolbulk)*(Na/1.0d24) 
expmupos = xposbulk /xsolbulk**vsalt
expmuneg = xnegbulk /xsolbulk**vsalt
expmuHplus = xHplusbulk /xsolbulk   ! vsol = vHplus 
expmuOHmin = xOHminbulk /xsolbulk   ! vsol = vOHmin 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Open common files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(unit=301, file='F_tot.dat')
open(unit=302, file='F_mixs.dat')
open(unit=303, file='F_mixpos.dat')
open(unit=304, file='F_mixneg.dat')
open(unit=305, file='F_mixH.dat')
open(unit=306, file='F_mixOH.dat')
open(unit=307, file='F_conf.dat')
open(unit=308, file='F_eq.dat')
open(unit=309, file='F_vdW.dat')
open(unit=311, file='F_electro.dat')
open(unit=312, file='F_tot2.dat')
open(unit=313, file='mupolA.dat')
open(unit=314, file='pilat.dat')
open(unit=315, file='mupolB.dat')

end subroutine
