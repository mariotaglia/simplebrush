subroutine saveall(cc,ccc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Saves results to disk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use chainsdat
use system
use results
use solver
use molecules
implicit none


character*5  title
integer cc, ccc
character*24 filename
integer i, n

n = dimz

! Polymer
title = 'avpol'
call savetodisk(avpol, title, cc ,ccc)

! Solvent
title = 'avsol'
call savetodisk(xh, title, cc, ccc)

! Cationes
! Cations
title = 'avpos'
call savetodisk(xpos, title, cc, ccc)

! Anions
title = 'avneg'
call savetodisk(xneg, title, cc, ccc)

! H+
title = 'avHpl'
call savetodisk(xHplus, title, cc, ccc)

! OH-
title = 'avOHm'
call savetodisk(xOHmin, title, cc,ccc)

! fdis
title = 'frdis'
call savetodisk(fdis, title, cc, ccc)

! Potencial electrostatico
title = 'poten'
call savetodisk(psi, title, cc, ccc)

! System information

 write(filename,'(A7, I3.3, A1, I3.3, A4)')'system.',cc,'.',ccc,'.dat'
 open (unit=510, file=filename)

         write(510,*)'st          = ',st ! residual size of iteration vector
         write(510,*)'fnorm       = ',norma ! residual size of iteration vector
         write(510,*)'length seg  = ',0.5 ! value see subroutine cadenas
         write(510,*)'delta       = ',delta
         write(510,*)'vsol        = ',vsol
         write(510,*)'vsol        = ',vpol
         write(510,*)'vsalt       = ',vsalt*vsol
         write(510,*)'csalt       = ',csalt
         write(510,*)'pHbulk      = ',pHbulk
         write(510,*)'pKa         = ',pKa
         write(510,*)'pK0         = ',-dlog(K0)/dlog(10.0D0)
         write(510,*)'K0          = ',K0
         write(510,*)'zpos        = ',zpos
         write(510,*)'zneg        = ',zneg
         write(510,*)'cuantas     = ',cuantas
         write(510,*)'iterations  = ',iter
         write(510,*)'sigma       = ',sigma

close(510)

! Saves solver vector
           
write(filename,'(A6, I3.3, A4)')'out.', ccc, '.dat'
open(unit=45, file=filename)
do i = 1, 2*n
  write(45, *)xflag(i)
enddo
close(45)

end subroutine

subroutine savetodisk(array, title, counter, counter2)

use system

integer scx, scy
integer dimzview

integer ix, iy, iz, i, jx, jy, jz
real*8 array(dimz)
real*8 arrayz(dimz)
integer counter, counter2
character*5 title
character*6 titlez
character*30 filename, tempc

write(filename,'(A5 , A1, I3.3, A1, I3.3, A4)')   & 
  title,'.', counter,'.', counter2, '.dat'
open(unit=45, file=filename)
do iz=1,dimz
  write(45,*)(iz-0.5)*delta,array(iz)
enddo
close(45)

return
end subroutine


