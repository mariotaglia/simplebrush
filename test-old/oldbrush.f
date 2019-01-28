!###############################################################################
!     
!     Simple brush: Standard Molecular Theory Program 
!    
!     Calculates a weak polyelectrolyte brush in poor sv conditions 
!     Calculates free-energy
!
!
!    
!###############################################################################

      implicit none

      real*8 pi
 
      real*8 Na
      parameter (Na=6.02d23)

      real*8 delta

      real*8 fmedio, sumpol

      real*8 pend, ord

      integer dimzview

      integer cuantas, dimz ,long

      parameter (cuantas=200000)
      parameter (long = 28)
      parameter (dimz = 40)

      real*8 in1(cuantas,long,3)  ! Posicion de cada segmento

      integer i,j, k, ix, iy, iz, ii
      integer infile, flag 

      character basura
      character*24 filename 
      character*5  title


      real*8 xtotal(0:dimz+1) ! xtotal para poor solvent
      real*8 norma, error

      real*8 looped

     

C-----  varables de la resolucion -----------

      real*8 x1(2*dimz),xg1(2*dimz)
      real*8 xflag(2*dimz), xflag2(2*dimz), algo
  
       
      real*8 errel, fnorm
      integer n, itmax, iter

      integer cc,ccc,cccc

      integer seed

C--------------------------------------------

! Volumen
      
      real*8 vsol
      real*8 vpol 
      real*8 vsalt, vsalt2
 

      real*8 eps(dimz), eps1

! Volumen fraction

      real*8 avpol(dimz), xh(dimz)
      real*8 qtot(dimz) ! Carga total
      real*8 psi(dimz) ! potencial
      real*8 xpos(dimz) ! pos ion
      real*8 xpos2(dimz) ! pos ion
      real*8 xneg(dimz) ! neg ioni
      real*8 xHplus(dimz) ! H+
      real*8 xOHmin(dimz) ! OH-
      
! Bulk

      real*8 xsolbulk ! volume fraction solvent in bulk
      real*8 xposbulk, xposbulk2, xnegbulk, xsalt,xsalt2,csalt,csalt2
      real*8 stok
      real*8 expmupos, expmuneg, expmupos2
         
      real*8 sigma, sigmadelta,  mupol1, mupol2, mupold
      real*8 sigma1, sigma2
      integer*1 sigmaflag

! Charge

      real*8 zpos, zpos2, zneg, zpol
      real*8 sigmaq
      real*8 lb
      real*8 constq
      real*8 betae

! Poor solvent

      real*8 st
      real*8 Xu, Au, Bu, Cu

! Weak pol

      real*8 fdis(dimz) 
      real*8 Kw, pKw, pKa, Ka, pHbulk, expmuHplus, expmuOHmin
      real*8 Kb, pKb
      real*8 xHplusbulk, cHplus, pOHbulk, xOHminbulk, cOHmin
      real*8 K0

! Free energy

      real*8 pro(cuantas), q

      real*8 Free_energy, F_Mix_s, F_Mix_pos, F_Mix_pos2
      real*8 F_Mix_neg, F_Mix_Hplus
      real*8 Free_energy2, sumpi, sumrho, sumel, sum, mupol, pilat
      real*8 F_Mix_OHmin, F_Conf, F_Eq, F_vdW, F_eps, F_electro
      real*8 medio, abajo

C--------------------------------------------

! Kinsol

      integer *4 ier ! Kinsol error flag
      integer *8 neq ! Kinsol number of equations
      real*8 f(2*dimz)

! IMSL

      external fcnelect

C---------------------------------------------

! Barre

      real*8 sts(100)
      integer nst


      integer nspl
      parameter (nspl=1000)     ! numero de puntos para spline cubico
    
      real*8 mudlist_spl(2,nspl)
      real*8 temp



      REAL*8 mudlist1(2,100) 
      REAL*8, ALLOCATABLE :: mudlist(:,:) ! Cargas capacitica y faradaica
      REAL*8 yspl(nspl), xspl(nspl)

      real*8 shift

      common /segme/ in1
      common /volum/ vpol, vsol, vsalt, vsalt2
      common /potent/ eps, constq, sigmaq, sigma
      common /resul1/ avpol,xpos,xpos2, xneg, qtot, xHplus, xOHmin, fdis
      common /charge/ zpol, zpos,zpos2, zneg

      common /delta/ delta
      common /bulk/ expmupos,expmupos2,expmuneg,expmuHplus,expmuOHmin,K0
      common /iter/ iter 
 
      common /ps/ st, abajo, medio      

      common /psize/ neq ! Kinsol

      external fcn

      common /norma/ norma

      common /seed1/ seed

      common /pro/ pro, q

      common /xtotal/ xtotal

      common /shift/ shift

c--------------------------------------------------


! Variables

      shift = 1.0

      pi=dacos(-1.0d0)          ! pi = arccos(-1) 

      delta = 0.5
      lb = 0.714 ! bjerrum lenght in nm

      zpos = 1.0
      zpos2 = 3.0
      zneg = -1.0
      
      vsol = 0.030     
      vsalt=((4.0/3.0)*pi*(0.27)**3)/vsol  ! volume salt in units of vsol 0.2=radius salt  
      vsalt2=((4.0/3.0)*pi*(0.27)**3)/vsol  ! volume salt in units of vsol 0.2=radius salt  
      vpol= 0.095/vsol! ((4.0/3.0)*pi*(0.2)**3)/vsol  ! volume polymer segment in units of vsol 

      constq=delta*delta*4.0*pi*lb/vsol   ! multiplicative factor in poisson eq  

      pKw = 14

      error = 1e-6 ! para comparar con la norma...

      betae = 38.94             ! beta * e

c----------------------------------------------

! Coeficientes poor-solvent

      call kai


c----------------------------------------------

      n=dimz
      
      errel=1d-6
      itmax=200

! Lee variables de resolucion

      read(8, *), basura
      read(8, *), sigma

      read(8, *), basura
      read(8, *),eps1
 
      read(8,*), basura
      read(8,*), sigmaq

      sigmaq = sigmaq *(4.0*pi*lb*delta)

      read(8,*), basura
      read(8, *), zpol

      read(8, *), basura
      read(8, *), csalt

      csalt2 = 0.0

      read(8, *), basura
      read(8, *), pKa

      read(8, *), basura
      read(8, *), pHbulk

      read(8, *), basura
      read(8, *), infile

      read(8, *), basura
      read(8, *), nst
      read(8, *), (sts(i), i=1, nst)

      read(8, *), basura
      read(8, *), sigmadelta

! Array eps

      eps(1) = 0.0

      do iz = 2, dimz

      eps(iz) = 0.0

      enddo

C---------------------------------------------------

! Crea cadenas.......

      call creador ! Genera cadenas
      call pxs ! Genera cadenas

c----------------------------------------------------------
! Archivos comunes a todos los casos...

c      real*8 Free_energy, F_Mix_s, F_Mix_pos, F_Mix_neg, F_Mix_Hplus
c      real*8 F_Mix_OHmin, F_Conf, F_Eq, F_vdW, F_eps, F_electro

c       open(unit=301, file='F_tot.dat')
c       open(unit=302, file='F_mixs.dat')
c       open(unit=303, file='F_mixpos.dat')
c       open(unit=304, file='F_mixneg.dat')
c       open(unit=305, file='F_mixH.dat')
c       open(unit=306, file='F_mixOH.dat')
c       open(unit=307, file='F_conf.dat')
c       open(unit=308, file='F_eq.dat')
c       open(unit=309, file='F_vdW.dat')
c       open(unit=310, file='F_eps.dat')
c       open(unit=311, file='F_electro.dat')
c       open(unit=312, file='F_tot2.dat')
c       open(unit=314, file='pilat.dat')
c       open(unit=315, file='F_mixpos2.dat')
            
C----------------------------------------------------------

! LOOP

C----------------------------------------------------------

      sigmaflag = 0
      sigma1 = sigma * vsol / (delta) ! dimensionless sigma
      sigma2 = (sigma+sigmadelta) * vsol / (delta) ! dimensionless sigma



     
      open(unit=301, file='mupol.dat')
      open(unit=302, file='mupold.dat')
      open(unit=303, file='mupol_spl.dat')
      open(unit=304, file='zeros')
      open(unit=305, file='fmedio')

      ccc = 1
      do while(ccc.le.nst) ! loop en st

 555  print*, 'ccc', ccc, 'de', nst
      print*, 'q',q

      print*, 'sigmaflag', sigmaflag

      st = sts(ccc)
      if(sigmaflag.eq.0) then
      sigma = sigma1
      else
      sigma = sigma2
      endif
        
   
  

! Variables dependientes del input...


 257  Ka=10**(-pKa)
      pKb = 14 - pKa
      Kb = 10**(-pKb)
      cHplus = 10**(-pHbulk)    ! concentration H+ in bulk
      xHplusbulk = (cHplus*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol
      pOHbulk= pKw -pHbulk
      cOHmin = 10**(-pOHbulk)   ! concentration OH- in bulk
      xOHminbulk = (cOHmin*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol  

      xsalt=(csalt*Na/(1.0d24))*(vsalt*vsol)   ! volume fraction salt,csalt in mol/l 
      xsalt2=(csalt2*Na/(1.0d24))*(vsalt2*vsol)   ! volume fraction salt,csalt in mol/l 

      if(pHbulk.le.7) then  ! pH<= 7
            xposbulk=xsalt/zpos
            xnegbulk=
     &  -xsalt/zneg - xsalt2/zneg +(xHplusbulk -xOHminbulk) *vsalt ! NaCl+ HCl  
      else                  ! pH >7 
            xposbulk=xsalt/zpos +(xOHminbulk -xHplusbulk) *vsalt ! NaCl+ NaOH   
            xnegbulk=-xsalt/zneg - xsalt2/zneg
      endif

      xposbulk2 = xsalt2/zpos2            ! Lantano o calcio
      

         xsolbulk=1.0 -xHplusbulk -xOHminbulk -
     &        xnegbulk -xposbulk - xposbulk2

         K0 = (Kb*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant 

         expmupos = xposbulk /xsolbulk**vsalt
         expmupos2 = xposbulk2 /xsolbulk**vsalt2
         expmuneg = xnegbulk /xsolbulk**vsalt
         expmuHplus = xHplusbulk /xsolbulk   ! vsol = vHplus 
         expmuOHmin = xOHminbulk /xsolbulk   ! vsol = vOHmin 

! Guess inicial

      if(infile.eq.2) then

      if((ccc.eq.1).and.(sigmaflag.eq.0)) then ! primer caso para cada pH = lee xfile
      do i = 1, 2*n ! dejar aca, asegura guess correcto para
      xg1(i) = xflag2(i)
      x1(i) = xflag2(i)
      enddo
      endif
      if(ccc.ne.1) then
      do i = 1, 2*n ! dejar aca, asegura guess correcto para infile = 2 (cc > 1)
      xg1(i) = xflag(i)     
      x1(i) = xflag(i)
      enddo
      endif
      endif

      if(infile.eq.0) then
      do i=1,n
         xg1(i)=xsolbulk
	 x1(i)=xsolbulk
      enddo

      do i=n+1, n*2
         xg1(i)=0.0d0
         x1(i)=0.0d0
      enddo

      endif
 
      if(infile.eq.1) then

      open(unit=45, file='in.txt')
      do i=1, 2*n
      read(45,*)xg1(i)
      x1(i) = xg1(i)
      enddo
      close(45)
      endif 
      

C--------------------------------------------------------------
C               +++ RESOLUCION +++ 
C--------------------------------------------------------------

          call call_kinsol(x1, xg1, ier)

! Recupera xh y psi (NO SON COMMON!)

            do iz=1,dimz

            xh(iz)=x1(iz)

            psi(iz)=x1(iz+n)

            enddo

! Fin de rutina del solver

! Chequea si exploto... => Sistema anti-crash


         if((ier.lt.0).or.
     &   (norma.gt.error)) then ! exploto...

         print*, 'Error en solver: ', ier
         print*, 'norma ', norma
         print*, 'q ', q
         print*, 'st de error', st
         print*, 'pH de error', pHbulk
       
         algo = (st + stok)/2

c         nsigma = nsigma + 1


          write(1010,*)'Fallo st ', st, ' Paso a ', algo

         st = algo
         flag = 1

         goto 257
         endif    

! No exploto, guardo xflag

         do i = 1, 2*n
         xflag(i) = x1(i) ! xflag sirve como input para la proxima iteracion
         enddo

         if((ccc.eq.1).and.(sigmaflag.eq.0)) then ! guarda para loop en pH

         do i = 1, 2*n
         xflag2(i) = x1(i) ! xflag sirve como input para la proxima iteracion
         enddo

         endif

         infile = 2 ! no vuelve a leer infile

         stOk = st

         if(flag.eq.1) then  ! habia un error...
         print*, 'Recupero del error'
         print*, 'OK', stOK
         flag = 0
         goto 555
         endif


C----------------------------------------------------------
C  OUTPUT!
C----------------------------------------------------------

      if(sigmaflag.eq.0) then

! Guarda infile


           
      write(filename,'(A6, I3.3, A4)')'out.', ccc, '.dat'
      open(unit=45, file=filename)
      do i = 1, 2*n
      write(45, *)x1(i)
      enddo
      close(45)

!!!!!!!!!!!!!!!!!!! Guarda archivos !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      cc = 1

! Polimero

      title = 'avpol'
      call savetodisk(avpol, title, cc ,ccc)

! Solvente

      title = 'avsol'
      call savetodisk(xh, title, cc, ccc)

! Cationes

      title = 'avpos'
      call savetodisk(xpos, title, cc, ccc)

! Cationes2 

      title = 'avpo2'
      call savetodisk(xpos2, title, cc, ccc)


! Aniones

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


!!!!!!!!!!!!!!!!! Informacion del sistema !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         write(filename,'(A7, I3.3, A1, I3.3, A4)')'system.',cc,'.',ccc,
     &     '.dat'
         open (unit=510, file=filename)

         write(510,*)'st          = ',st ! residual size of iteration vector
         write(510,*)'fnorm       = ',norma ! residual size of iteration vector
         write(510,*)'length seg  = ',0.5 ! value see subroutine cadenas
         write(510,*)'delta       = ',delta
         write(510,*)'vsol        = ',vsol
         write(510,*)'vsol        = ',vpol
         write(510,*)'vsalt       = ',vsalt*vsol
         write(510,*)'vsalt2       = ',vsalt2*vsol
         write(510,*)'csalt       = ',csalt
         write(510,*)'pHbulk      = ',pHbulk
         write(510,*)'pKw         = ',pKw
         write(510,*)'pKa         = ',pKa
         write(510,*)'pK0         = ',-dlog(K0)/dlog(10.0D0)
         write(510,*)'K0          = ',K0
         write(510,*)'sigmaq      = ',sigmaq /(4.0*pi*lb*delta),sigmaq
         write(510,*)'zpos        = ',zpos
         write(510,*)'zpos2        = ',zpos2
         write(510,*)'zneg        = ',zneg
         write(510,*)'cuantas     = ',cuantas
         write(510,*)'iterations  = ',iter
         write(510,*)'sigma       = ',sigma

         close(510)
 
         endif ! sigmaflag  

 


c----------------------------------------------------------------------------------------------
c        Calcula Energia libre
C----------------------------------------------------------------------------------------------

      Free_Energy = 0.0
      Free_Energy2 = 0.0

! 1. Mezcla solvente

      F_Mix_s = 0.0 

      do iz = 1, dimz
      F_Mix_s = F_Mix_s + xh(iz)*(dlog(xh(iz))-1.0)
      F_Mix_s = F_Mix_s - xsolbulk*(dlog(xsolbulk)-1.0)
      enddo      

      F_Mix_s = F_Mix_s * delta/vsol
      Free_Energy = Free_Energy + F_Mix_s

! 2. Mezcla ion positivo

      F_Mix_pos = 0.0 

      do iz = 1, dimz

      F_Mix_pos = F_Mix_pos + xpos(iz)*(dlog(xpos(iz)/vsalt)-1.0  
     &  - dlog(expmupos) + dlog(vsalt))

      F_Mix_pos = F_Mix_pos - xposbulk*(dlog(xposbulk/vsalt)-1.0
     &  - dlog(expmupos) + dlog(vsalt))

      enddo
      F_Mix_pos = F_Mix_pos * delta/vsol/vsalt
      Free_Energy = Free_Energy + F_Mix_pos

! 2. Mezcla ion positivo 2

C      F_Mix_pos2 = 0.0
C
c      do iz = 1, dimz
c      F_Mix_pos2 = F_Mix_pos2 + xpos2(iz)*(dlog(xpos2(iz)/vsalt2)-1.0
c     &  - dlog(expmupos2) + dlog(vsalt2))
c
c      F_Mix_pos2 = F_Mix_pos2 - xposbulk2*(dlog(xposbulk2/vsalt2)-1.0
c     &  - dlog(expmupos2) + dlog(vsalt2))
c
c      enddo
c
c      F_Mix_pos2 = F_Mix_pos2 * delta/vsol/vsalt2
c      Free_Energy = Free_Energy + F_Mix_pos2

! 3. Mezcla ion negativo

      F_Mix_neg = 0.0

      do iz = 1, dimz
      F_Mix_neg = F_Mix_neg + xneg(iz)*(dlog(xneg(iz)/vsalt)-1.0
     & - dlog(expmuneg) + dlog(vsalt))

      F_Mix_neg = F_Mix_neg - xnegbulk*(dlog(xnegbulk/vsalt)-1.0
     & - dlog(expmuneg) + dlog(vsalt))

      enddo 
      F_Mix_neg = F_Mix_neg * delta/vsol/vsalt
      Free_Energy = Free_Energy + F_Mix_neg

! 4. Mezcla protones

      F_Mix_Hplus = 0.0

      do iz = 1, dimz
      F_Mix_Hplus = F_Mix_Hplus + xHplus(iz)*(dlog(xHplus(iz))-1.0 
     & -dlog(expmuHplus))

      F_Mix_Hplus = F_Mix_Hplus - xHplusbulk*(dlog(xHplusbulk)-1.0
     & -dlog(expmuHplus))

      enddo
      F_Mix_Hplus = F_Mix_Hplus * delta/vsol
      Free_Energy = Free_Energy + F_Mix_Hplus

! 5. Mezcla hidroxilos

      F_Mix_OHmin = 0.0

      do iz = 1, dimz
      F_Mix_OHmin = F_Mix_OHmin + xOHmin(iz)*(dlog(xOHmin(iz))-1.0
     & -dlog(expmuOHmin))

      F_Mix_OHmin = F_Mix_OHmin - xOHminbulk*(dlog(xOHminbulk)-1.0
     & -dlog(expmuOHmin))

      enddo
      F_Mix_OHmin = F_Mix_OHmin * delta/vsol
      Free_Energy = Free_Energy + F_Mix_OHmin

! 6. Entropia interna polimero

      F_Conf = 0.0

      do i = 1, cuantas

         F_Conf = F_Conf + (pro(i)/q)*dlog((pro(i))/q)
     &  /vsol*delta*sigma

      enddo

      Free_Energy = Free_Energy + F_Conf

! 7. Chemical Equilibrium

      F_Eq = 0.0 
            

      do iz  = 1, dimz

      F_Eq = F_Eq + fdis(iz)*dlog(fdis(iz))*avpol(iz)/vpol
      F_Eq = F_Eq + (1.0-fdis(iz))*dlog(1.0-fdis(iz))*avpol(iz)/vpol
      F_Eq = F_Eq + (1.0-fdis(iz))*dlog(K0)*avpol(iz)/vpol
      F_Eq = F_Eq + (1.0-fdis(iz))*(-dlog(expmuOHmin))*avpol(iz)/vpol

      enddo

      F_eq = F_eq *delta/vsol

      Free_Energy = Free_Energy + F_Eq

! 8.vdW ! Ojo, los kai son negativos => atraccion

       F_vdW = 0.0

      do iz = 1, dimz
      F_vdW = F_vdW - 0.5000*delta*xtotal(iz)*(xtotal(iz)*medio 
     &  + (xtotal(iz+1)+xtotal(iz-1))*abajo)*st

      enddo


      Free_Energy = Free_Energy + F_vdW

! 9. Electrostatic ! VER ESTO...

      F_electro = 0.0    

c      do iz  = 1, dimz
c
c      F_electro = F_electro + delta*psi(iz)*qtot(iz)/2.0/vsol
c
c      enddo

c      F_electro = F_electro + sigmaq*psi(0)/2.0

       do iz = 2, dimz-1
  
      F_electro = F_electro + delta*psi(iz)*qtot(iz)
      F_electro = F_electro -0.5/constq
     &   *(psi(iz+1)-psi(iz))*(psi(iz)-psi(iz-1))*delta

      enddo


      F_electro = F_electro + delta*psi(1)*qtot(1)
      F_electro = F_electro + delta*psi(dimz)*qtot(dimz)

      F_electro = F_electro -0.5/constq
     &   *(psi(2)-psi(1))*(sigmaq)*delta


      Free_Energy = Free_Energy + F_electro

! 10. Pol-Sup

      F_eps = 0.0 
      do iz = 1, dimz
      F_eps = F_eps - eps(iz)
      enddo
      Free_Energy = Free_Energy + F_eps

      print*, 'A', F_eq
      print*, '1', Free_energy
 

! Segun Peng-Gong




      Free_Energy2 = 0.0


         sumpi = 0.0
         sumrho=0.0
         sumel=0.0



         do i=1,dimz
            
            sumpi = sumpi+dlog(xh(i))     
            sumpi = sumpi-dlog(xsolbulk)     

            sumrho = sumrho + ( - xh(i) -xHplus(i) -xOHmin(i) 
     &               -(xpos(i)+xneg(i))/vsalt - xpos2(i)/vsalt2 )! sum over  rho_i i=+,-,si
            sumrho = sumrho - ( - xsolbulk -xHplusbulk -xOHminbulk
     &               -(xposbulk+xnegbulk)/vsalt - xposbulk2/vsalt2)! sum over  rho_i i=+,-,si

            sumel = sumel - qtot(i)*psi(i)/2.0       
                                ! electrostatic part free energy
         
         enddo
         
         sumpi = (delta/vsol)*sumpi
         sumrho = (delta/vsol)*sumrho
         sumel = (delta/vsol)*sumel

         sum = sumpi + sumrho + sumel

         Free_Energy2 = -(delta/(vsol))*sigma*dlog(q) + sum -F_vdW 

         Print*, '2',Free_Energy2

! Calcula mupol

         mupol = -dlog(q/shift)
         print*, 'mupol', mupol 

         if(sigmaflag.eq.0) then
         mupol1 = mupol
         else
         mupol2 = mupol
 
         mupold = (mupol2-mupol1)/sigma 


         write(301, *)st, mupol
         write(302, *)st, mupold

         mudlist1(2,ccc) = mupold
         mudlist1(1,ccc) = st

         endif

! Calcula pilat

         pilat = sigma*delta/vsol*mupol - Free_energy


! Guarda energia libre


c      real*8 Free_energy, F_Mix_s, F_Mix_pos, F_Mix_neg, F_Mix_Hplus
c      real*8 F_Mix_OHmin, F_Conf, F_Eq, F_vdW, F_eps, F_electro

 
c         write(301,*)csalt, pH, Free_energy
c         write(302,*)csalt, pH, F_Mix_s 
c         write(303,*)looped, F_Mix_pos
c         write(303,*)looped, F_Mix_pos2
c         write(304,*)looped, F_Mix_neg
c         write(305,*)looped, F_Mix_Hplus
c         write(306,*)looped, F_Mix_OHmin
c         write(307,*)looped, F_Conf
c         write(308,*)looped, F_Eq
c         write(309,*)looped, F_vdW
c         write(310,*)looped, F_eps
c         write(311,*)looped, F_electro
c         write(312,*)looped, Free_energy2
c         write(313,*)looped, mupol
c         write(314,*)looped, pilat


c--------------------------- FIN DE ENERGIA LIBRE -----------------
   
! fmedio
  
      sumpol = 0.0
      fmedio = 0.0

      do iz = 1, dimz
      fmedio = fmedio + avpol(iz)*fdis(iz)
      sumpol = sumpol + avpol(iz)
      enddo
 
      fmedio = fmedio/sumpol

      write(305,*)st, fmedio

      if(sigmaflag.eq.0) then
      sigmaflag = 1
      goto 555
      endif

      sigmaflag = 0 

      ccc = ccc + 1

      enddo ! ccc

      ALLOCATE (mudlist(2,nst))       ! dimensiona arrays de voltametria (Pot, Charge)

      do i = 1, nst
      mudlist(1, i)=mudlist1(1,i)
      mudlist(2, i)=mudlist1(2,i)
      enddo

      call fitspl (yspl, xspl, mudlist, nst)
      do i = 1, nspl
      mudlist_spl(1,i) = xspl(i)
      mudlist_spl(2,i) = yspl(i)
      WRITE (303,*)xspl(i), yspl(i)
      end do

c Busca ceros

      do i = 1, nspl-1
      temp = mudlist_spl(2,i)/mudlist_spl(2,i+1)
      if(temp.lt.0.0) then                      ! uno positivo y el otro negativo
      WRITE (304,*)pHbulk,(mudlist_spl(1, i)+mudlist_spl(1, i+1))/2
      endif
      enddo

! close files

c      close(301)
c      close(302)
c      close(303)
c      close(304)
c      close(305)
c     close(306)
c      close(307)
c      close(308)
c      close(309)
c      close(310)
c      close(311)
c      close(312)
c      close(313)
c      close(314)
c      close(315)
c
 
!      call average(xh) ! Calcula el volume fraction de las cadenas sin usar CPC

!!!!!!!!!!!!!!!!!!! Matriz en formato vtk !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      stop
      end
      
C*************************************************************

      subroutine fkfun(x,f,ier)

      implicit none
 

      integer cuantas, dimz ,long
      parameter (cuantas=200000)
      parameter (long = 28)
      parameter (dimz = 40)

      integer ntot
      real*8 x(2*dimz),f(2*dimz)

      real*8 in1(cuantas,long,3)  ! Posicion de cada segmento

      integer*1 pz(cuantas,long)

      real*8 q, pro(cuantas), protemp, protemp1

      integer i,j, k, ix, iy, iz, ii, ax, ay, az, temp

      real*8 temp2

      integer iter
      real*8 norma
      real*8 sigma
      real*8 delta


C-------------------------------------

! Volumen

      real*8 vsol
      real*8 vpol
      real*8 vsalt
      real*8 vsalt2

      real*8 eps(dimz)

! Volume fraction

      real*8 xh(dimz)
      real*8 avpol(dimz)
      real*8 xpot(dimz)

      real*8 qtot(dimz)
      real*8 xpos(dimz) ! pos ion
      real*8 xpos2(dimz) ! pos ion
      real*8 xneg(dimz) ! neg ion
      real*8 xHplus(dimz) ! H+
      real*8 xOHmin(dimz) ! OH-

      real*8 psi(0:dimz+1) ! psi se define asi para evitar problemas al construir las fs
       
! Bulk 

      real*8 expmupos,expmupos2,  expmuneg

! Charge

      real*8 zpos,zpos2, zneg, zpol
      real*8 sigmaq
      real*8 constq

! poor solvent 

      real*8 st
      real*8 xtotal(0:dimz+1) ! xtotal para poor solvent
      real*8 Xu, Au, Bu, Cu, medio, abajo

 
!  weakpol

      real*8 fdis(dimz) 
      real*8 expmuHplus,expmuOHmin, K0
      
! Kinsol

      integer*8 neq
      integer*4 ier

      real*8 shift

C-----------------------------------------------------

! Common variables

      common /segme/ in1
      common /volum/ vpol, vsol, vsalt, vsalt2
      common /potent/ eps, constq, sigmaq, sigma
      common /resul1/ avpol, xpos, xpos2,xneg, qtot,xHplus, xOHmin, fdis
      common /bulk/ expmupos,expmupos2,expmuneg,expmuHplus,expmuOHmin,K0
      common /charge/ zpol, zpos,zpos2, zneg

      common /chainz/ pz
      common /iter/ iter

      common /ps/ st, abajo, medio      

      common /norma/ norma
      common /delta/ delta

      common /pro/ pro, q
      common /xtotal/ xtotal

      common /shift/ shift

! Recupera xh y psi desde x()

      ntot = dimz 

            do iz=1,dimz

            xh(iz)=x(iz)
            psi(iz)=x(iz+ntot)
      
            enddo
      
! Condicion de borde psi en z


      psi(dimz+1) = 0.0 ! psibulk = 0.0
      psi(0) = sigmaq + psi(1) !psimetal ! potencial constante

      sigmaq = psi(0) - psi(1)

! Fracciones de volumen inicial	y fdis

            do iz=1,dimz

           avpol(iz)=0

           xpos(iz) = expmupos*(xh(iz)**vsalt)
     &     *dexp(-psi(iz)*zpos) ! ion plus volume fraction

           xpos2(iz) = expmupos2*(xh(iz)**vsalt2)
     &     *dexp(-psi(iz)*zpos2) ! ion plus volume fraction

           xneg(iz) = expmuneg*(xh(iz)**vsalt)
     &     *dexp(-psi(iz)*zneg) ! ion neg volume fraction
           xHplus(iz) = expmuHplus*(xh(iz))
     &     *dexp(-psi(iz))           ! H+ volume fraction
           xOHmin(iz) = expmuOHmin*(xh(iz))
     &     *dexp(+psi(iz))           ! OH-  volume fraction
           fdis(iz) =
     &     1.0 /(1.0 + xOHmin(iz)/(K0*xh(iz)) )
   
            enddo

! Calculo de xtotal para poor solvent

! en el lattice

            do iz=1,dimz

          xtotal(iz) = 1.0 - xpos(iz) - xneg(iz) - xpos2(iz) 
     &    - xh(iz) - xHplus(iz) - xOHmin(iz) ! xtotal es todo menos solvente e iones

            enddo

      xtotal(dimz+1) = 0.0 ! xtotal en bulk = 0.0
      xtotal(0) = 0.0 ! xtotal en la superficie = 0.0


      q = 0.0

      do iz = 1, dimz
c      xpot(iz) = xh(iz)**vpol* 
c     & dexp(eps(iz)-psi(iz)*zpol)


      xpot(iz) = xh(iz)**vpol*
     & dexp(eps(iz)-psi(iz)*zpol+st*vsol*vpol*
     & (medio*xtotal(iz)+abajo*(xtotal(iz+1)+xtotal(iz-1))))
     & /fdis(iz)

      enddo

      do i=1,cuantas

         pro(i)=1.0*shift

         do j=1,long

         az = pz(i, j)         
         pro(i) = pro(i) * xpot(az)
  
         enddo
            q=q+pro(i)

c            print*,i, pro

            do j = 1,long
            az = pz(i, j)
            avpol(az) = avpol(az) + pro(i)*sigma*vpol
            end do

         enddo


         do iz=1, dimz            ! norma avpol
               avpol(iz)=avpol(iz)/q
         enddo

         temp2 = 0.0

         do iz = 1, dimz
         temp2 = temp2 + avpol(iz)
         enddo
c         print*, temp2, sigma, vpol


C----------------------------------------------------------------------------------------------
C   Construye Ecuaciones a resolver 
C----------------------------------------------------------------------------------------------

! Qtot

      do iz=1,dimz

         qtot(iz) = 
     &   (zpos*xpos(iz)+zneg*xneg(iz))/vsalt + zpos2*xpos2(iz)/vsalt2+
     &   avpol(iz)*zpol/vpol*fdis(iz) + xHplus(iz)-xOHmin(iz) 

      enddo

! Volume fraction

      do iz=1,dimz
               
               f(iz)=
     &  avpol(iz) + xh(iz) + xneg(iz) 
     &  + xpos(iz) + xpos2(iz) + xHplus(iz) + xOHmin(iz) - 1.000000d0
   
      enddo

! Poisson eq.

            do iz=1,dimz

               f(iz+ntot)=
     & psi(iz+1) -2*psi(iz) + psi(iz-1) +
     & qtot(iz)*constq

               f(iz+ntot)=f(iz+ntot)/(-2.0)
      enddo
 
      iter = iter + 1

      norma = 0.0

      do i = 1, 2*ntot
      norma = norma +(f(i))**2    
      enddo

      print*, iter, norma, q
      ier = 0.0
   
c      do i = 1, dimz
c      print*,i, xpos2(i)
c      enddo


 
      return

 

      end

C*************************************************************



C*************************************************************

C**************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!  Esta subrutina se encarga de poner a todas los segmentos dentro del slab
!
! 


      subroutine pxs

      implicit none
    
      integer cuantas, dimz ,long
      real*8 delta 

      parameter (cuantas=200000)
      parameter (long = 28)
      parameter (dimz = 40)

      real*8 in1(cuantas,long,3)  ! Posicion de cada segmento

      integer i,j, k, ix, iy, iz, ii

      integer*1 pz(cuantas,long)

      common /segme/ in1
      common /chainz/ pz
      common /delta/delta



      do i=1,cuantas

         do j=1,long

            pz(i,j)=int((in1(i,j,1))/delta)+1

            if(pz(i, j).gt.dimz) then
                Print*,'Incrementar dimz', pz(i, j)
                stop
            endif         
        
            
 
         enddo
         
c         print*, i, pz(i, 10)

      enddo

      return
      end
      

      subroutine creador

      implicit none

      real*8 lseg,delta

      integer long,k,vx(4),vy(4)

      integer cuantas,seed,total,ix(3)

      parameter (lseg=0.5,
     &     cuantas=200000)

      parameter (long = 28)


      integer i,il,u1,u2,iii,ii,ll
      integer j,ncha

      real*8 indax, inday, indaz
      
      real*8 xend(3,200),rands,chains(3,200,100)
      real*8 altx,alty,altz,x(200),y(200),xp(200),yp(200)
      real*8 rij,theta,theta1,pi, rn1, rn2
      
      real*8 in1(cuantas,long,3)

      integer total1,iglobal

      common /segme/ in1
      common /seed1/ seed

      pi=dacos(-1.0d0)

      il=0
      iglobal=1

      theta1= 0 !30.0*pi/180.0
     
      do while (il.lt.cuantas)
         
         call cadenas72mr(chains,ncha)
         
         do i=1,ncha

            do ll=1,1

               il=il+1

               if(il.gt.cuantas) goto 100

               theta=ll*theta1

c               rn1=(rands(seed)-0.5) !OJO
c               rn2=(rands(seed)-0.5) 


               do j=1,long

                  xp(j)=chains(2,j,i)
                  yp(j)=chains(3,j,i)
               
                  ! Genero las cadenas en 0,0

                  x(j)=xp(j)*cos(theta)+yp(j)*sin(theta) !+ rn1
                  y(j)=-xp(j)*sin(theta)+yp(j)*cos(theta) !+ rn2

                  altx=x(j)
                  alty=y(j)
                  altz=chains(1,j,i)
                  
                  indax=altx
                  inday=alty
                  indaz=altz
                  
                
         
                  in1(il,j,2)=indax
                  in1(il,j,3)=inday
                  in1(il,j,1)=indaz
               
               enddo
               
            enddo

         enddo

      enddo
      

 100  return

!  if(rank.eq.0) then
! print*,'# de cadenas (pretencion)',cuantas
! print*,'total',il-1
! print*,'SEED de salida-->',seed
! print*,'numero fort.??-->',iii
! print*,seed,iii+1,total1
! endif
 
      !return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC





C****************************************************************

      subroutine cadenas72mr(chains,ncha)
      
      implicit none
	
      integer long,seed,ncha
      real*8 chains(3,200,100)
     
      parameter (long = 28)

 
      integer i,state,ii,j,ive,jve
      real*8 rn,state1,pi,sitheta,cotheta,dista,lseg

      real*8 siphip,cophip,rands
      character*1 test
      real*8 m(3,3),mm(3,3),tt(3,3),tp(3,3),tm(3,3)
      
      parameter (lseg=0.5)
      real*8 x(3),xend(3,200),xendr(3,200)

      common /seed1/ seed
      
      pi=acos(-1.0000000e0)
      sitheta=sin(68.0*pi/180.0)
      cotheta=cos(68.0*pi/180.0)
      siphip=sin(120.0*pi/180.0)
      cophip=cos(120.0*pi/180.0)
      
 223  x(1)=lseg
      x(2)=0.0
      x(3)=0.0
      
      xend(1,1)=lseg
      xend(2,1)=0.0
      xend(3,1)=0.0
      
      tt(1,1)=cotheta
      tt(1,2)=sitheta
      tt(1,3)=0.0
      tt(2,1)=sitheta
      tt(2,2)=-cotheta
      tt(2,3)=0.0
      tt(3,1)=0.0
      tt(3,2)=0.0
      tt(3,3)=-1.0
      
      tp(1,1)=cotheta
      tp(1,2)=sitheta
      tp(1,3)=0.0
      tp(2,1)=sitheta*cophip
      tp(2,2)=-cotheta*cophip
      tp(2,3)=siphip
      tp(3,1)=sitheta*siphip
      tp(3,2)=-cotheta*siphip
      tp(3,3)=-cophip
      
      tm(1,1)=cotheta
      tm(1,2)=sitheta
      tm(1,3)=0.0
      tm(2,1)=sitheta*cophip
      tm(2,2)=-cotheta*cophip
      tm(2,3)=-siphip
      tm(3,1)=-sitheta*siphip
      tm(3,2)=cotheta*siphip
      tm(3,3)=-cophip
      
 222  rn=rands(seed)
      
      state1=0.0
c     state1=rn*pi*2.0
c     print *,state1,'state1'
      
      m(1,1)=cotheta
      m(1,2)=sitheta
      m(1,3)=0.0
      
      m(2,1)=cos(state1)*sitheta
      m(2,2)=-cos(state1)*cotheta
      m(2,3)=sin(state1)
      m(3,1)=sin(state1)*sitheta
      m(3,2)=-sin(state1)*cotheta
      m(3,3)=-cos(state1)
      
      x(1)=m(1,1)*lseg
      x(2)=m(2,1)*lseg
      x(3)=m(3,1)*lseg
      
      xend(1,2)=lseg+x(1)
      xend(2,2)=x(2)
      xend(3,2)=x(3)
      
      do 10 i=3,long
         
 123     rn=rands(seed)
         state=int(rn*3)
c     print*,'state',state

         if (state.eq.3) then 
            state=2
         endif
         if (state.eq.0) then
c     statelog='t' 
            
            call mrrrr(m,tt,mm)
            do 30 ii=1,3
               do 40 j=1,3
                  m(ii,j)=mm(ii,j)
 40            continue
 30         continue
            
            
         elseif (state.eq.1) then
c     rn=rands(seed)
c     if (rn.gt.0.25) goto 123
c     statelog='+'
            
            call mrrrr(m,tp,mm)
            do 31 ii=1,3
               do 41 j=1,3
                  m(ii,j)=mm(ii,j)
 41            continue
 31         continue

         elseif (state.eq.2) then
c     rn=rands(seed)
c     if (rn.gt.0.25) goto 123
c     statelog='-'
            
            call mrrrr(m,tm,mm)
            do 32 ii=1,3
               do 42 j=1,3
                  m(ii,j)=mm(ii,j)
 42            continue
 32         continue
            
         endif
         
         x(1)=m(1,1)*lseg
         x(2)=m(2,1)*lseg
         x(3)=m(3,1)*lseg
         
         xend(1,i)=xend(1,i-1)+x(1)
         xend(2,i)=xend(2,i-1)+x(2)
         xend(3,i)=xend(3,i-1)+x(3)
         
c     if (xend(1,i).lt.0.0) goto 222
         
 10   continue
      
      dista=0.0
      do 300 ive=4,long
         do 310 jve=1,ive-3
            dista=(xend(1,jve)-xend(1,ive))**(2.0)
            dista=dista+(xend(2,jve)-xend(2,ive))**(2.0)
            dista=dista+(xend(3,jve)-xend(3,ive))**(2.0)
            dista=sqrt(dista)
            if (dista.lt.lseg) then
cccccccccwrite(2,*)'noself-a.',k
               goto 222
            endif
 310     continue
 300  continue

      ncha=0
      do 400 i=1,300

         test='S'
         call rota36(xend,xendr,long,test)
         if (test.eq.'N') goto 400
         ncha=ncha+1
         do 401 j=1,long
            chains(1,j,ncha)=xendr(1,j)
            chains(2,j,ncha)=xendr(2,j)
            chains(3,j,ncha)=xendr(3,j)
 401     continue
         if (ncha.eq.25) goto 402
         
 400  continue
 402  if (ncha.eq.0) goto 223
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine rota36(xend,xendr,n,test)
      
      real*8 xend(3,200),rands,xendr(3,200)
      character*1 test
      integer n,seed
      real*8 fac,fac1,fac2,sbe,cbe,sal,cal,sga
      
      common /seed1/ seed
      
      pi=acos(-1.0)

      fac=rands(seed)
      fac1=rands(seed)
      fac2=rands(seed)
      alfa=fac*2*pi
      cbe=2.0d0*fac1-1.0d0
      gama=fac2*2*pi

      sbe=(1-cbe**2)**0.5
      cal=cos(alfa)
      sal=sin(alfa)
      cga=cos(gama)
      sga=sin(gama)

      do 1 i=1,n

         a=xend(1,i)
         b=xend(2,i)
         c=xend(3,i)

         xendr(1,i)=a*(-cbe*sal*sga+cal*cga)
     &        -b*(cbe*sal*cga+cal*sga)+c*sbe*sal
         xendr(2,i)=a*(cbe*cal*sga+sal*cga)+
     &        b*(cbe*cal*cga-sal*sga)-c*sbe*cal
         xendr(3,i)=a*sbe*sga+b*sbe*cga+c*cbe
 
 1    continue

      do 2 i=1,n
         if (xendr(1,i).lt.0.0) test='N'  
 2    continue

      return
      end
      
C****************************************************************

      subroutine mrrrr(a,b,c)

      real*8 a(3,3),b(3,3),c(3,3)

      do 1 i=1,3
         do 1 j=1,3
            c(i,j)=0
 1    continue

      do 2 i=1,3
         do 2 j=1,3
            do 2 k=1,3
               c(i,j)=c(i,j)+a(i,k)*b(k,j)
 2    continue

      return
      end 
        
C****************************************************************
C **********************************************************************
        double precision FUNCTION RANDS (SEED)
C **********************************************************************

C-----  this is a special function for random number generation
C        on 32-bit machines that do not support long integer
C        multiplication and truncation.  the technique used is to do
C        the multiplication and addition in parts, by splitting all
C       integers in a 'high' and a 'low' part.  the algorithm is
C-----        exact, and should give machine-independent results.

C-----        the algorithm implemented is (following d.e. knuth):
C        seed = seed*1592653589 + 453816691
C        if (seed.lt.0) seed = seed + 1 + 2147483647
C-----        note that 1592653589 = 48603*2**15 + 30485

C 32768 = 2^15, 65536 = 2^16, 4.65661287308E-10 = 2^(-31)

        INTEGER SEED, I1, I2

        I2 = MOD (SEED, 32768) * 30485
        I1 = MOD (SEED / 32768 * 30485, 65536) + MOD (MOD (SEED, 32768)
     X    * 48603, 65536) + I2 / 32768 + MOD (SEED / 32768, 2) * 32768 +
     X     13849
        I2 = MOD (I2, 32768) + 12659
        I1 = I1 + I2 / 32768
        I2 = MOD (I2, 32768)
        I1 = MOD (I1, 65536)
        SEED = I1 * 32768 + I2
        RANDS = REAL(I1 * 32768 + I2) * 4.65661287308E-10

        RETURN
        END


c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* * *
c     The routine kpsol is the preconditioner solve routine. It must have
c     that specific name be used in order that the c code can find and link
c     to it.  The argument list must also be as illustrated below:

      subroutine fkpsol(udata, uscale, fdata, fscale, vv, ftem, ier)

      implicit none

      integer  dimz
      parameter (dimz = 40)

      integer ier
      integer*8 neq, i
      double precision pp
      double precision udata(*), uscale(*), fdata(*), fscale(*)
      double precision vv(*), ftem(*)

      common /pcom/ pp(2*dimz)
      common /psize/ neq

      do  i = 1, neq
         vv(i) = vv(i) * pp(i)
      enddo
      ier = 0

      return
      end

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* * *
c     The routine kpreco is the preconditioner setup routine. It must have
c     that specific name be used in order that the c code can find and link
c     to it.  The argument list must also be as illustrated below:

      subroutine fkpset(udata, uscale, fdata, fscale,vtemp1,vtemp2, ier)

      implicit none

      integer dimz 
      parameter (dimz = 40)


      integer ier
      integer*8 neq, i
      double precision pp
      double precision udata(*), uscale(*), fdata(*), fscale(*)
      double precision vtemp1(*), vtemp2(*)


      common /pcom/ pp(2*dimz)
      common /psize/ neq

      do i = 1, neq/2
         pp(i) = 0.1 / (1.0+(1.0-udata(i))*exp(1.0-udata(i)))
      enddo

      do i = neq/2+1, neq
         pp(i) = 1.0
      enddo

      ier = 0


      return
      end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subrutina que llama a kinsol
!

      subroutine call_kinsol(x1, xg1, ier)

      integer dimz 
      parameter (dimz = 40)

      real*8 x1(2*dimz), xg1(2*dimz)

      integer*8 iout(15) ! Kinsol additional output information
      real*8 rout(2) ! Kinsol additional out information

      integer*8 msbpre
      real*8 fnormtol, scsteptol
      real*8 scale(2*dimz)
      real*8 constr(2*dimz)

      integer*4  globalstrat, maxl, maxlrst

      integer *4 ier ! Kinsol error flag
      integer *8 neq ! Kinsol number of equations

      common /psize/ neq ! Kinsol


! INICIA KINSOL



      neq = 2*dimz
      msbpre  = 5 ! maximum number of iterations without prec. setup (?)
      fnormtol = 1.0d-8 ! Function-norm stopping tolerance
      scsteptol = 1.0d-8 ! Function-norm stopping tolerance

      maxl = 10 ! maximum Krylov subspace dimesion (?!?!?!) ! Esto se usa para el preconditioner
      maxlrst = 2 ! maximum number of restarts

      globalstrat = 0

      call fnvinits(3, neq, ier) ! fnvinits inits NVECTOR module
      if (ier .ne. 0) then       ! 3 for Kinsol, neq ecuantion number, ier error flag (0 is OK)
      print*, 'SUNDIALS_ERROR: FNVINITS returned IER = ', ier
         stop

      endif

      call fkinmalloc(iout, rout, ier)    ! Allocates memory and output additional information
      if (ier .ne. 0) then
         print*, 'SUNDIALS_ERROR: FKINMALLOC returned IER = ', ier
         stop
      endif


      call fkinsetiin('MAX_SETUPS', msbpre, ier)  ! Additional input information
      call fkinsetrin('FNORM_TOL', fnormtol, ier)
      call fkinsetrin('SSTEP_TOL', scsteptol, ier)

         do i = 1, neq  !constraint vector
         constr(i) = 0.0
         enddo
         do i = 1, neq/2  !constraint vector
         constr(i) = 2.0
         enddo



         call fkinsetvin('CONSTR_VEC', constr, ier) ! constraint vector

c       CALL FKINSPTFQMR (MAXL, IER)
c      call fkinspgmr(maxl, maxlrst, ier) !  Scale Preconditioned GMRES solution of linear system (???)
      call fkinspbcg(maxl, ier) !  Scale Preconditioned BCG


      if (ier .ne. 0) then
      print*, 'SUNDIALS_ERROR: FKINSPGMR returned IER = ', ier
         call fkinfree ! libera memoria
         stop
      endif
      call fkinspilssetprec(1, ier) ! preconditiones

      do i = 1, neq ! scaling vector
      scale(i) = 1.0
      enddo

      do i = 1, neq ! Initial guess
      x1(i) = xg1(i)
      enddo


      call fkinsol(x1, globalstrat, scale, scale, ier)         ! Llama a kinsol
      if (ier .lt. 0) then
      print*, 'SUNDIALS_ERROR: FKINSOL returned IER = ', ier
      print*, 'Linear Solver returned IER = ', iout(9)
c         call fkinfree

      endif

      return
      end


      subroutine savetodisk(array, title, counter, counter2)

      integer dimz
      parameter (dimz = 40)

      integer scx, scy
      integer dimzview

      integer ix, iy, iz, i, jx, jy, jz
      real*8 array(dimz)
      real*8 arrayz(dimz)
      integer counter, counter2
      character*5 title
      character*6 titlez
      character*30 filename, tempc
      real*8 delta

      common /delta/ delta


      ! Variables

      dimzview = 20 ! maximo en z

      ! Graba material en crudo
c
c      write(filename,'(A5, A1, I3.3, A1, I3.3, A4)') 
c     &   title,'.', counter,'.', counter2, '.raw' 
c      open(unit=45, file=filename)
c
c          do iz = 1, dimz
c            write(45,*)array(iz)
c          enddo
c      close (45)

      ! Integra y graba promedios en z

      titlez = title // 'z'
      write(filename,'(A6 , A1, I3.3, A1, I3.3, A4)')
     &  titlez,'.', counter,'.', counter2, '.dat' 
      open(unit=45, file=filename)
      do iz=1,dimz
         write(45,*)(iz-0.5)*delta,array(iz)
      enddo
      close(45)

      return
      end

      subroutine fitspl (yspl, xspl, mulist, points)

      integer points, nspl, i, n
      parameter (nspl=1000)
      REAL*8 yspl(nspl), xspl(nspl), yp1, ypn
      REAL*8 ireal, nsplreal
      REAL*8 mulist(2,nspl), x(nspl), y(nspl), y2(nspl)  ! Dimensiona por exceso

      n = points - 1

      do i = 1, n
      x(i) = mulist(1, i)
      y(i) = mulist(2, i)
      end do

      yp1 = 1E31
      ypn = 1E31
      call splines(x, y, n, yp1, ypn, y2)

      nsplreal = nspl
      do  i = 1, nspl
      ireal = i
      xspl(i) = ((ireal-1)/(nsplreal-1))*(x(n)-x(1))+x(1)
      call splint(x, y, y2, n, xspl(i), yspl(i))
      end do
      end


      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL*8 x,y,xa(n),y2a(n),ya(n)
c Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
c xai~Rs in order), and given the array y2a(1:n), which is the output from spline above,
cand given a value of x, this routine returns a cubic-spline interpolated value y.
      INTEGER k,khi,klo
      REAL*8 a,b,h
      klo=1
c We will find the right place in the table by means of bisection.
cThis is optimal if sequential calls to this routine are at random values of x. If sequential calls are in order, and closely
c spaced, one would do better to store previous values of
c klo and khi and test if they remain appropriate on the
c next call.
      khi=n
    1 if (khi-klo.gt.1) then
      k=(khi+klo)/2
      if(xa(k).gt.x)then
      khi=k
      else
      klo=k
      endif
      goto 1
      ENDIF ! klo and khi now bracket the input value of x.
      h=xa(khi)-xa(klo)
!      if (h.eq.0.) pause !~Rbad xa input in splint~R The xa~Rs must be distinct.
      a=(xa(khi)-x)/h !Cubic spline polynomial is now evaluated.
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     & ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      END



      SUBROUTINE splines(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL*8 yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)

      INTEGER i,k
      REAL*8 p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then !The lower boundary condition is set either to be natural
      y2(1)=0.
      u(1)=0.
      else !or else to have a specified first derivative.
      y2(1)=-0.5
      u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1 !This is the decomposition loop of the tridiagonal algorithm. y2 and u are used for temporary storage of the decomposed factors.  
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*y2(i-1)+2.
      y2(i)=(sig-1.)/p
      u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     & /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn.gt..99e30) then !The upper boundary condition is set either to be ~Snatural~T
      qn=0.
      un=0.
      else !or else to have a specified first derivative.
      qn=0.5
      un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1 ! This is the backsubstitution loop of the tridiagonal algorithm
      y2(k)=y2(k)*y2(k+1)+u(k)
      enddo
      return
      END



      subroutine kai

      implicit none


      real*8 suma
      real*8 lseg ! largo del segmento

      parameter (lseg=0.5)
      real*8 delta ! delta


      integer dimR
      integer seed

      integer MCsteps ! numero de steps de MC

      real*8 R,theta,z, x, y

      real*8 vsol, vpol

      real*8 rn

      integer i, ii

      real*8 rands

      real*8 pi

      real*8 medio, abajo

      real*8 st
      real*8 radio ! radio del poro en nm

      real*8 x1,x2,y1, y2, z1, z2, vect

       integer j

C---------------------------------------------

      common /delta/ delta
      common /ps/ st, abajo, medio
      common /nanoporo/ radio


      print*,'Kai calculation'
      open(unit=111, file='kais.dat')

      pi=dacos(-1.0d0)          ! pi = arccos(-1) 

      suma = 0.0

      abajo = 0.0
      medio = 0.0

      seed = 1010
c      MCsteps = 1000000000
      MCsteps = 60**3

      do i = 1, MCsteps

      z = 3.0*(rands(seed)-0.5)*delta ! numero al azar entre -1.5*delta y 1.5*delta
      x = 3.0*(rands(seed)-0.5)*delta ! numero al azar entre -1.5*delta y 1.5*delta
      y = 3.0*(rands(seed)-0.5)*delta ! numero al azar entre -1.5*delta y 1.5*delta

      vect = sqrt(x**2 + y**2 + z**2) ! vector 

! coordenadas del segmento (x1,y1,z1) y del punto a integrar (x2,y2,z2)

      if(vect.gt.(1.5*delta))goto 15 ! No esta dentro de la esfera del cut-off   
      if(vect.lt.lseg) goto 15 ! esta dentro de la esfera del segmento

      if(z.gt.(delta/2.0))abajo = abajo+ ((lseg/vect)**6) 
      if((z.lt.(delta/2.0)).and.(z.gt.(-delta/2.0)))
     & medio = medio + ((lseg/vect)**6) 


 15   enddo

      medio = medio/MCsteps*(3.0*delta)**3
      abajo = abajo/MCsteps*(3.0*delta)**3

       print*, 'abajo', abajo
       print*, 'medio', medio


      close(111)

      end


