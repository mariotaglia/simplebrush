subroutine read

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
read(8, *), sigma  ! surface coverage of the chains

read(8, *), basura
read(8, *), csalt  ! salt concentration in bulk (Molar)

read(8, *), basura
read(8, *), pKa    ! pKa of weak polyacid segments

read(8, *), basura
read(8, *), pHbulk ! bulk pH

read(8, *), basura
read(8, *), st     ! polymer-polymer attraction strenght in kBT
 
read(8, *), basura
read(8, *), Xulimit  ! cutoff for porr sv interaction in lattice sites

read(8, *), basura
read(8, *), infile ! read input from file?

end subroutine
