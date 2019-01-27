!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     The routine kpsol is the preconditioner solve routine. It must have
!     that specific name be used in order that the c code can find and link
!     to it.  The argument list must also be as illustrated below:

      subroutine fkpsol(udata, uscale, fdata, fscale, vv, ftem, ier)

      use chainsdat
      implicit none

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

!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     The routine kpreco is the preconditioner setup routine. It must have
!     that specific name be used in order that the c code can find and link
!     to it.  The argument list must also be as illustrated below:

      subroutine fkpset(udata, uscale, fdata, fscale,vtemp1,vtemp2, ier)

      use chainsdat
      implicit none

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
! This routine calls kinsol solver
!

      subroutine call_kinsol(x1, xg1, ier)

      
      use chainsdat
      implict none

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

! INIT KINSOL

      neq = 2*dimz
      msbpre  = 5 ! maximum number of iterations without prec. setup (?)
      fnormtol = 1.0d-8 ! Function-norm stopping tolerance
      scsteptol = 1.0d-8 ! Function-norm stopping tolerance

      maxl = 500 ! maximum Krylov subspace dimesion
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

!       CALL FKINSPTFQMR (MAXL, IER)
!       call fkinspgmr(maxl, maxlrst, ier) !  Scale Preconditioned GMRES solution of linear system 
        call fkinspbcg(maxl, ier) !  Scale Preconditioned BCG


      if (ier .ne. 0) then
      print*, 'SUNDIALS_ERROR: FKINSPGMR returned IER = ', ier
         call fkinfree ! free memory
         stop
      endif
      call fkinspilssetprec(1, ier) ! preconditioner

      do i = 1, neq ! scaling vector
      scale(i) = 1.0
      enddo

      do i = 1, neq ! Initial guess
      x1(i) = xg1(i)
      enddo


      call fkinsol(x1, globalstrat, scale, scale, ier)         ! call kinsol
      if (ier .lt. 0) then
      print*, 'SUNDIALS_ERROR: FKINSOL returned IER = ', ier
      print*, 'Linear Solver returned IER = ', iout(9) 

      call fkinfree ! free memory

      endif

      return
      end subroutine



