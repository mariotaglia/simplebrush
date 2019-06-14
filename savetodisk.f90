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

! Polymer A
title = 'avpoA'
call savetodisk(avpolA, title, cc ,ccc)

! Polymer B
title = 'avpoB'
call savetodisk(avpolB, title, cc ,ccc)

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

! fdisA
title = 'fdANC'
call savetodisk(fdisANC, title, cc, ccc)

! fdisA
title = 'fdANa'
call savetodisk(fdisANa, title, cc, ccc)

! fdisB
title = 'fdBNC'
call savetodisk(fdisBNC, title, cc, ccc)

! fdisAas
title = 'fdAas'
call savetodisk(fdisAas, title, cc, ccc)

! fdisBas
title = 'fdBas'
call savetodisk(fdisBas, title, cc, ccc)

! fdisA
title = 'fdBCl'
call savetodisk(fdisBCl, title, cc, ccc)

! mfdisAas
title = 'mfdis'
call savetodisk(M, title, cc, ccc)

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
         write(510,*)'vpol        = ',vpol
         write(510,*)'vsalt       = ',vsalt*vsol
         write(510,*)'csalt       = ',csalt
         write(510,*)'pHbulk      = ',pHbulk
         write(510,*)'pKaA        = ',pKaA
         write(510,*)'pK0A         = ',-dlog(K0A)/dlog(10.0D0)
         write(510,*)'K0A          = ',K0A
			write(510,*)'pKEo         = ',pKEo	
         write(510,*)'K0Eo          = ',K0Eo
			write(510,*)'pKaB         = ',pKaB
         write(510,*)'pK0B         = ',-dlog(K0B)/dlog(10.0D0)
         write(510,*)'K0B          = ',K0B	 
			write(510,*)'zpos        = ',zpos
         write(510,*)'zneg        = ',zneg
         write(510,*)'cuantas     = ',cuantas
         write(510,*)'newcuantas     = ',newcuantas
         write(510,*)'iterations  = ',iter
         write(510,*)'sigmaA       = ',sigmaA
         write(510,*)'sigmaB      = ',sigmaB
			write(510,*)'Version GIT	  = ',_VERSION
close(510)

! Saves solver vector
           
write(filename,'(A6, I3.3, A4)')'outin.', ccc, '.dat'
open(unit=45, file=filename)
do i = 1, 3*n
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


