!###############################################################################
!     
!     Simple brush: Standard Molecular Theory Program 
!    
!     Calculates a weak polyelectrolyte brush in poor sv conditions 
!     Calculates free-energy
!
!     MT Jan2019
!         
!###############################################################################

implicit none
integer cc, ccc

print*, 'Program Simple Brush'
print*, 'GIT Version: ', _VERSION

call readinput  ! reads input variables from file
call init       ! initialize system dependent variables
call allocation ! allocates memory
call kais       ! generates coefficients for poor-sv interactions
call creador    ! create chains conformations

cc = 1
ccc = 1

! START HERE LOOPING OVER VARAIBLES TO SCAN, USE CC AND CCC as COUNTERS

call solve               ! solves the molecular theory
call fe(cc, ccc)         ! calculates and saves free energy to disk
call saveall(cc, ccc)    ! save results to disk

! END HERE LOOP

call endall     ! clean up and terminate

end

subroutine endall
stop
end subroutine
