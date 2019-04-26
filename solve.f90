subroutine solve

use system
use solver
use bulk
implicit none

integer n
integer i
real*8 x1(3*dimz),xg1(3*dimz)
integer ier

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Set ups the initial guess and calls kinsol
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

n = dimz ! size of system

! Set ups initial guess 

if(infile.eq.2) then ! aleready solve one case, recover guess from xflag
  do i = 1, 3*n 					!!G_
    xg1(i) = xflag(i)     
    x1(i) = xflag(i)
  enddo
endif

if(infile.eq.0) then ! no initial guess provided, use bulk
  do i=1,n
    xg1(i)=xsolbulk
    x1(i)=xsolbulk
  enddo

  do i=n+1, n*2
     xg1(i)=0.0d0
     x1(i)=0.0d0
  enddo

  do i=n+n+1, n*3				!!G_
     xg1(i)=0.0
     x1(i)=0.0
  enddo

endif
 



if(infile.eq.1) then ! initial guess provided in file in.txt
  open(unit=45, file='in.txt')
  do i=1, 3*n					!!G_
    read(45,*)xg1(i)
    x1(i) = xg1(i)
  enddo
  close(45)
endif 

! Solve the molecular theory      

call call_kinsol(x1, xg1, ier)

! Checks for convergence 

if((ier.lt.0).or.(norma.gt.error)) then ! failed...

         print*, 'Error in solver: ', ier
         print*, 'norm ', norma
         print*, 'st', st
         print*, 'pH', pHbulk
         call endall
endif    

! Converged ok
! stores xflag

xflag = x1 ! xflag serves as a initial guess for next iteratio
infile = 2 ! avoids reading infile again

return
end subroutine


