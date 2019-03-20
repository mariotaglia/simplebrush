subroutine creador

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calls cadenas to generate the set of random conformations and then calls pxs
! to put them in the lattice
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use system
use chainsdat
implicit none

integer i, il
real*8 chains(3,long,100)
integer ncha
real*8 in1(long,3)

integer vx(4),vy(4)
integer total,ix(3)
integer u1,u2,iii,ii,ll
real*8 indax, inday, indaz
real*8 xend(3,long)
real*8 altx,alty,altz,x(long),y(long),xp(long),yp(long)
real*8 rij,theta,theta1,pi, rn1, rn2
integer total1,iglobal
integer j, k
il=0
newcuantas=0
do while (il.lt.cuantas)
         
call cadenas(chains,ncha)
         
 do i=1,ncha
  il=il+1
  if(il.gt.cuantas) exit
  do j = 1,long
   do k = 1,3
     in1(j,k)=chains(k,j,ncha)
   enddo
  enddo
  call pxs(in1,il)
 enddo ! ncha
enddo ! while
return
end subroutine


subroutine cadenas(chains,ncha)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This routine generates the polymer conformations
!  
!  Outputs: chains(3,:,:) coordinates of a subset of conformations in real space
!           ncha          number of conformations in the subset
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


use chainsdat
use rand
use const
implicit none

integer ncha
real*8 chains(3,long,100)
integer i,state,ii,j,ive,jve
real*8 rn, state1,sitheta,cotheta,dista
real*8 siphip,cophip,rands
logical test

logical testsa ! self avoiding?

real*8 m(3,3),mm(3,3),tt(3,3),tp(3,3),tm(3,3)
real*8 x(3),xend(3,long),xendr(3,long)

! Initialize angles of rotation matrix
      
sitheta=sin(68.0*pi/180.0)
cotheta=cos(68.0*pi/180.0)
siphip=sin(120.0*pi/180.0)
cophip=cos(120.0*pi/180.0)
 
ncha = 0
     
do while (ncha.eq.0)

 x(1)=lseg
 x(2)=0.0
 x(3)=0.0
      
! first segment, note that in this code the first segment is not on the surface but at a distance lseg from it
 xend(1,1)=lseg
 xend(2,1)=0.0
 xend(3,1)=0.0

! rotation matrixes
     
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
    

 testsa = .false. ! not self-avoiding
  
 do while (testsa.eqv..false.) ! loop until a self-avoiding conformation is found

! second segment

  rn=rands(seed)
  state1=0.0
      
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
      
! other segments

  do i=3,long
      rn=rands(seed)
      state=int(rn*3)

      select case (state)
         case (0)
            call mrrrr(m,tt,mm)
         case (1)
            call mrrrr(m,tp,mm)
         case (2, 3) ! state shouldn't be 3, but just in case rounding errors give 3
            call mrrrr(m,tm,mm)
      end select

      m = mm
         
      x(1)=m(1,1)*lseg
      x(2)=m(2,1)*lseg
      x(3)=m(3,1)*lseg
         
      xend(1,i)=xend(1,i-1)+x(1)
      xend(2,i)=xend(2,i-1)+x(2)
      xend(3,i)=xend(3,i-1)+x(3)
  enddo ! i

!!!! Check for self-avoidance       
         
      testsa=.true.
      dista=0.0
      do ive=4,long
         do jve=1,ive-3
            dista=(xend(1,jve)-xend(1,ive))**(2.0)
            dista=dista+(xend(2,jve)-xend(2,ive))**(2.0)
            dista=dista+(xend(3,jve)-xend(3,ive))**(2.0)
            dista=sqrt(dista)
            if (dista.lt.lseg)testsa=.false. ! not self-avoding
         enddo
       enddo

 enddo ! self-avoding?

 ncha=0 ! number of chains

 do i=1,300
   test = .true. ! avoidance with the wall
   call rota36(xend,xendr,long,test) ! rotate conformation
   if (test.eqv..false.)cycle
         ncha=ncha+1
         chains(1,:,ncha)=xendr(1,:)
         if (ncha.eq.25)exit
 enddo
enddo ! outer loop ncha  

return
end

subroutine rota36(xend,xendr,n,test)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generates a ramdom rotations of conformation xend
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

use rand
use chainsdat
use const
implicit none

real*8 xend(3,long),rands,xendr(3,long)
logical test
integer n
real*8 fac,fac1,fac2,sbe,cbe,sal,cal,sga
real*8 a,b,c   
real*8 alfa, cga, gama
integer i
 
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

do i=1,n

   a=xend(1,i)
   b=xend(2,i)
   c=xend(3,i)

   xendr(1,i)=a*(-cbe*sal*sga+cal*cga)-b*(cbe*sal*cga+cal*sga)+c*sbe*sal
   xendr(2,i)=a*(cbe*cal*sga+sal*cga)+b*(cbe*cal*cga-sal*sga)-c*sbe*cal
   xendr(3,i)=a*sbe*sga+b*sbe*cga+c*cbe

enddo
do i=1,n
   if (xendr(1,i).lt.0.0) test=.false.
enddo

return
end
      
!****************************************************************

subroutine mrrrr(a,b,c)

implicit none
real*8 a(3,3),b(3,3),c(3,3)
integer i,j,k

do i=1,3
   do j=1,3
      c(i,j)=0
   enddo
enddo

do i=1,3
 do j=1,3
   do k=1,3
      c(i,j)=c(i,j)+a(i,k)*b(k,j)
   enddo
 enddo
enddo

return
end subroutine
