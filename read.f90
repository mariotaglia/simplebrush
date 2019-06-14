subroutine readinput

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!  This routine reads variables from fort.8
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use chainsdat
use system
use solver
use kai
implicit none

character basura 

! read starts here, not that read is performed sequentially! 
! dimzA=dimzB=dimZ
! deltaA=deltaB=delta
! cuantasA=cuantasB=cuantas
! longA=longB=long 
! lsegA=lsegB=lseg
! sigmaA=sigmaB=sigma
! pkaA=!pkaB. Me parece que no necesariamente son iguales para los dos polimeros 
! stA=stB=st Uso el mismo 

read(8,*), basura
read(8,*), dimz    ! number of lattice sites

read(8,*), basura
read(8,*), delta   ! size of lattice site in nm

read(8,*), basura
read(8,*), cuantas ! number of polymer conformations

read(8,*), basura
read(8,*), long    ! lenght of polyelectrolyte conformations

read(8,*), basura
read(8,*), lseg    ! lenght of a polyelectrolyte segment in nm

read(8, *), basura
read(8, *), sigmaA  ! surface coverage of the chainsA	

read(8, *), basura
read(8, *), sigmaB  ! surface coverage of the chainsB	

read(8, *), basura
read(8, *), csalt  ! salt concentration in bulk (Molar)

read(8, *), basura
read(8, *), pKaA    ! pKaA of weak polyacid segments

read(8, *), basura
read(8, *), pKaANa    ! pKaA of weak polyacid segments

read(8, *), basura
read(8, *), pKaB    ! pKaB of weak polyacid segments

read(8, *), basura
read(8, *), pKaBCl    ! pKaB of weak polyacid segments

read(8, *), basura
read(8, *), pHbulk ! bulk pH

read(8, *), basura
read(8, *), st     ! polymer-polymer attraction strenght in kBT

read(8, *), basura					!!!!!!!!!!!!!!!!!!!!!!!
read(8, *), pKEo     ! Interation 	!!!!!!!!!!!!!!!!!!!!!!!

read(8, *), basura
read(8, *), Xulimit  ! cutoff for porr sv interaction in lattice sites

read(8, *), basura
read(8, *), infile ! read input from file?

end subroutine
