      SUBROUTINE PARDISM(DISIT,REFOR,LNODS,
     .                  NDOFN,NELEM,NEVAB,NNODE,NPOIN,NTOTV,
     .                  IDSK2,INRHS,IPASS,KPASS,KPORE,KRESL,NEQNS,NLAST,
     .                  IFFIX,
     .                  GSTDI,GSTLO,GSTUP,CSTIF,ELOAD,LNUEQ,LPONT,
     .                  DISIM,FOREL,ALOAD,DELTA,DISIC,LOCAL,
     .                  iiter,miter,
     .                  NPREL,NGRUP,NPROP,NMATS,
     .                  NDATA,NPREV,NSTAT,NMATX,NDIME,
     .                  NDISR,NDISO,
     .                  MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .                  ELVAR,ELMAT,WORK1,TEMPN,DISPR,DTEMP,
     .                  VNORM,DISTO,COORD,NPOIC,IAUGM,INFRI,COFRI,
     .                  NACTI,LACTI)
C***********************************************************************
C
C**** THIS ROUTINE SOLVES THE SET OF LINEAR EQUATIONS USING THE
C         Direct SOLVER PARDISO 
C
C     SKYDES
C     PARDISM
C            - PARASSM
C            - SKYRHS 
C         !  - SKYBAK              ! Direct Skyline Solver
C         !  - SKYITE  - SKYBAK    ! Iterative PCG Solver
C                      - SKYPRO
C            - SKYRHS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
c      parameter (nprdssize=20000000)
      INCLUDE 'mkl_pardiso.f77'
C
      COMMON/CONTNCM/NOCOI,NNODN,LARGC,LICOI,NSKIC !Agregado para evaluar contacto no coincidente
      COMMON/SOLVERA/KRENU,KSOLV,KSYMM,NWIDT,MITCG,NBUFA,NPRIR,
     .               MITGM,MKRYL,IPGMR,KPARAM, NTHMCM, NTHSOM
      COMMON/SOLVERB/TOLCG,TOLC1,TOLGM
C
      COMMON/LOGUN/LUDTS,LUSOL,LUFRO,LUFRH,LUDAT,LUPRI,LURES,
     .             LUSO2,LUFR2,LUPOS,LURST,LUBFG,LUPIP,LUPAN,LUINF,
     .             LUGEO,LUSET,LUMAT,LUINI,LULOA,LUFIX,LUIN1,
     .             LUTUN,LUCON,LUACT,LUFAN,
     .             LUCU1,LUCU2,LUCU3,LUCU4,LUCU5,LUCU6,LUCU7,
     .             LUCU8,LUCU9,LUC10
C
      DIMENSION DISIT(NTOTV),   REFOR(NTOTV),  LNODS(NNODE,NELEM)
      DIMENSION LNUEQ(NTOTV),   LPONT(NTOTV)      
      DIMENSION GSTDI(*),       GSTLO(*),      GSTUP(*), 
     .          CSTIF(NEVAB,*), ELOAD(*),      DISIM(*),
     .          FOREL(*),       ALOAD(*),      DELTA(*),
     .          DISIC(*),       LOCAL(*)
      DIMENSION IFFIX(*)
C
      DIMENSION MATNO(NELEM),
     .          PROEL(NPREL,NGRUP), PROPS(NPROP,NMATS),
     .          ELDAT(NDATA),       ELPRE(NPREV),
     .          ELVAR(NSTAT),       ELMAT(NMATX),
     .          WORK1(*)
      DIMENSION TEMPN(NPOIN,2),     DISPR(NTOTV,NDISR),
     .          DTEMP(NPOIN),       VNORM(NTOTV)
      DIMENSION DISTO(NTOTV,NDISO), COORD(NDIME,NPOIN)
      DIMENSION INFRI(*),           COFRI(*),
     .          LACTI(NELEM)
     
      INTERFACE
       SUBROUTINE  PARASSM(CSTIF,VALUES,LCOLUMNS,IRWNDX,
     .             LNUEQ,KSYMM,
     .             LNODS,NDOFN,NELEM,NEQNS,NEVAB,NLAST,NNODE,
     .             NWIDT,NTOTV,NPOIN,
     .             NPREL,NGRUP,NPROP,NMATS,
     .             NDATA,NPREV,NSTAT,NMATX,NDIME,
     .             NDISR,NDISO,
     .             MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .             ELVAR,ELMAT,WORK1,TEMPN,DISPR,DTEMP,
     .             VNORM,DISTO,COORD,INFRI,COFRI,LACTI,
     .             NPRDSSIZE,IPRDSINIT,CHANGEPHASE)
       IMPLICIT REAL*8(A-H,O-Z)
       REAL*8, ALLOCATABLE, INTENT(INOUT) :: values(:)
       INTEGER*4, ALLOCATABLE, INTENT(INOUT) :: IRWNDX(:)
       INTEGER*4, ALLOCATABLE, INTENT(INOUT) :: LCOLUMNS(:)
       INTEGER*4 CHANGEPHASE
      INCLUDE 'addi_om.f'
      INCLUDE 'auxl_om.f'
      DIMENSION CSTIF(NEVAB,NEVAB), LNODS(NNODE,*), LNUEQ(*)
C
      DIMENSION MATNO(NELEM),
     .          PROEL(NPREL,NGRUP), PROPS(NPROP,NMATS),
     .          ELDAT(NDATA),       ELPRE(NPREV),
     .          ELVAR(NSTAT),       ELMAT(NMATX),
     .          WORK1(*)
      DIMENSION TEMPN(NPOIN,2),     DISPR(NTOTV,NDISR),
     .          DTEMP(NPOIN),       VNORM(NTOTV)
      DIMENSION DISTO(NTOTV,NDISO), COORD(NDIME,NPOIN)
      DIMENSION INFRI(*),           COFRI(*),
     .          LACTI(NELEM)
       END SUBROUTINE
      END INTERFACE
      
      real*8 x(NEQNS)
      real*8, ALLOCATABLE :: values(:)
      INTEGER*4, ALLOCATABLE :: IRWNDX(:)
      INTEGER*4, ALLOCATABLE :: LCOLUMNS(:)

      INTEGER pt(64), iprdsinit,nprdssize,CHANGEPHASE
      INTEGER iparm(64)
      save pt, iprdsinit, IRWNDX, lcolumns,iprdsinitflag
      data iprdsinitflag /0/
     
      NPRDSSIZE=NEQNS**1.3
      IERROR=0
      IF (.NOT. ALLOCATED(LCOLUMNS)) THEN 
       ALLOCATE(LCOLUMNS(NPRDSSIZE),STAT=IERROR)
      ELSE
       NPRDSSIZE=SIZE(LCOLUMNS)
      END IF
      IF (IERROR.NE.0) 
     . CALL RUNEND("I CAN'T ALLOCATE LCOLUMNS IN PARDISM ")
      
      IERROR=0
      IF (.NOT. ALLOCATED(values)) 
     . ALLOCATE(values(nprdssize),STAT=IERROR)
      IF (IERROR.NE.0) 
     . CALL RUNEND("I CAN'T ALLOCATE values IN PARDISM ")
     
      IF (.NOT. ALLOCATED(IRWNDX)) 
     . ALLOCATE(IRWNDX(NEQNS+1),STAT=IERROR)
      IF (IERROR.NE.0) 
     . CALL RUNEND("I CAN'T ALLOCATE IRWNDX IN PARDISM ") 

      IF(IPRDSINITFLAG.EQ.0) THEN
        iprdsinit=0
       IPRDSINITFLAG=1
      end if     
     
     
C
C**** IDENTIFY THE DESTINATION OF THE EQUATIONS INTO PROFILE
C     NEQNS,NLAST,LNUEQ,LPONT
C
      IF((KPASS.EQ.0).OR.((KPORE.EQ.2).AND.(NEQNS.GT.0))) THEN
       REWIND IDSK2
       READ(IDSK2) NEQNS,LNUEQ
       KPASS=1
      ENDIF
C-------------------------------------------      
C      write(6,*)'NTHSOM=',NTHSOM
C-------------------------------------------
C
C**** ASSEMBLE & FACTORISE MATRIX IF NECESSARY
C
      IF(KRESL.EQ.1) THEN
       CALL PARASSM(CSTIF,VALUES,LCOLUMNS,IRWNDX,
     .             LNUEQ,KSYMM,
     .             LNODS,NDOFN,NELEM,NEQNS,NEVAB,NLAST,NNODE,
     .             NWIDT,NTOTV,NPOIN,
     .             NPREL,NGRUP,NPROP,NMATS,
     .             NDATA,NPREV,NSTAT,NMATX,NDIME,
     .             NDISR,NDISO,
     .             MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .             ELVAR,ELMAT,WORK1,TEMPN,DISPR,DTEMP,
     .             VNORM,DISTO,COORD,INFRI,COFRI,LACTI,
     .             NPRDSSIZE,IPRDSINIT,CHANGEPHASE)
C
      ENDIF
C
C**** STORE OR READ BACK THE FACTORISED GLOBAL STIFFNESS MATRICES
C

C****************************************************************
C      IF(KPORE.EQ.2) THEN
C       IF(NEQNS.GT.0)
C     .  CALL SKYTAP(KRESL,KSYMM,GSTDI,GSTLO,GSTUP,IDSK2,NLAST,
C     .              NEQNS)
C      ENDIF
C****************************************************************

C
C**** IF NECESSARY MODIFY RHS FOR NON-ZERO BOUNDARY CONDITIONS
C
      IF(INRHS.EQ.1)
     . CALL SKYRHS(CSTIF,DISIT,REFOR,LNODS,NDOFN,NELEM,NEVAB,
     .             NNODE,NPOIN,DISIM,FOREL,    1,
     .             NPREL,NGRUP,NPROP,NMATS,
     .             NDATA,NPREV,NSTAT,NMATX,NTOTV,NDIME,
     .             NDISR,NDISO,
     .             MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .             ELVAR,ELMAT,WORK1,TEMPN,DISPR,DTEMP,
     .             VNORM,DISTO,COORD,INFRI,COFRI,LACTI)
C
C**** COMPRESS RHS
C
      DO ITOTV=1,NTOTV
       IEQNS=LNUEQ(ITOTV)
       IF(IEQNS.GT.0) ELOAD(IEQNS)=REFOR(ITOTV)
      ENDDO
C
C**** SOLVE THE EQUATIONS 
C
! PARDISO
C------------------------------------------------------------------
         maxfct = 1
         mnum = 1
         nrhs = 1
      do i = 1, 64
         iparm(i) = 0
      end do
      iparm(1) = 1 ! no solver default
      iparm(2) = 2 ! fill-in reordering from METIS
      iparm(3) = NTHSOM ! numbers of processors - to be replaced by 4/NTHSO
      iparm(4) = 0 ! no iterative-direct algorithm
      iparm(5) = 0 ! no user fill-in reducing permutation
      iparm(6) = 1 ! overwrite b with x
      iparm(7) = 0 ! not in use
      iparm(8) = 1 ! numbers of iterative refinement steps
      iparm(9) = 0 ! not in use
      iparm(10) = 13 ! perturbe the pivot elements with 1E-13
      iparm(11) = 0 ! don't use nonsymmetric permutation and scaling MPS
      if(KSYMM.EQ.0) iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
      iparm(12) = 0 ! not in use
      iparm(13) = 0
      if(KSYMM.EQ.0) iparm(13) = 1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric)
      iparm(14) = 0 ! Output: number of perturbed pivots
      iparm(15) = 0 ! not in use
      iparm(16) = 0 ! not in use
      iparm(17) = 0 ! not in use
      iparm(18) = -1 ! Output: number of nonzeros in the factor LU
      iparm(19) = -1 ! Output: Mflops for LU factorization
      iparm(20) = 0 ! Output: Numbers of CG Iterations
      ierror = 0 ! initialize error flag
      msglvl = 0 ! print statistical information
      mtype = -2 ! real symmetric indefinite
      if(KSYMM.EQ.0) mtype = 11 ! real unsymmetric
      call mkl_set_num_threads(NTHSOM)     !to be replaced by 4/NTHSO
      if (iprdsinit .eq. 0) then
C         write (*,*) 'Initializing PARDISO solver'
         iprdsinit = 1
C.. Initiliaze the internal solver memory pointer. This is only
C necessary for the FIRST call of the PARDISO solver.
         do i = 1, 64
          pt(i) = 0
         end do
C.. Reordering and Symbolic Factorization, This step also allocates
C all memory that is necessary for the factorization
         iphase = 11 ! analysis only
         CALL pardiso (pt(1), maxfct, mnum, mtype, iphase, 
     .        neqns, VALUES, irwndx(1), lcolumns(1),
     .        idum, nrhs, iparm, msglvl, eload, x, ierror)
!$    timend=omp_get_wtime()
C      write (6,*)'Pardiso 11 takes',timend-timstart
      endif
!$    timstart=omp_get_wtime()

c     CALL CHECKSM (ELOAD,NEQNSF,COEFSUM)
c     write (6,*) 'kprob, eload antes de pardiso 23: ',
c    .     kprob,coefsum

      ncoef=IRWNDX(NEQNS+1)-1
C      CALL CHECKSM (VALUES,NCOEF,COEFSUM)
C      write (6,*) 'ncoef, loc-values, antes 23,coefsum: ',
C     .      ncoef,loc(values),coefsum,values(1),values(ncoef)

      IF(CHANGEPHASE.EQ.1) THEN
      iphase = 13 ! Analiza para la factorización de una nueva matriz, factoriza y resueleve.
      ELSE
      iphase = 23 ! numerical factorization, solve & iterative refinement 
      ENDIF
      
      CALL pardiso (pt(1), maxfct, mnum, mtype, iphase, 
     .     neqns, VALUES, irwndx(1), lcolumns(1),
     .     idum, nrhs, iparm, msglvl, eload, x, ierror)
     
!$    timend=omp_get_wtime()
C    write (6,*)'Pardiso 23 takes',timend-timstart

c     isumirow=0
c     ncoef=IRWNDX(NEQNSF+1,kprob)-1
c     do k=1,neqnsf+1
c       isumirow = isumirow + irwndx(k,kprob)
c     end do
c     write (6,*) 'kprob,ncoef,IRWNDX: ', kprob,ncoef,isumirow

C      CALL CHECKSM (VALUES,NCOEF,COEFSUM)
C      write (6,*) 'ncoef, loc-values, despues 23 , coefsum: ',
C     .      ncoef,loc(values),coefsum,values(1),values(ncoef)

c     CALL CHECKSM (ELOAD,NEQNSF,COEFSUM)
c     write (6,*) 'kprob, eload despues de pardiso 23: ',
c    .     kprob,coefsum

      if (ierror.ne.0) then
       write(LURES,*) 'pardiso error=', ierror
       stop
      end if

      if (iprdsinit .eq. 0) then   !comentar iprdsinit=1
        ddum = 0.0d0
        iphase=-1
        CALL pardiso (pt(1), maxfct, mnum, mtype, iphase, 
     .       neqns, ddum, idum, idum,
     .       idum, nrhs, iparm, msglvl, ddum, ddum, ierror)
      end if

C
      IF (ierror .NE. 0) THEN
         WRITE(LURES,*) 'The following ERROR was detected: ', ierror
         STOP
      END IF
C--------------------------------------------------------------------
C
C**** REALLOCATE SOLUTION ARRAY INTO TOTAL NTOTV
C
      DO ITOTV=1,NTOTV
       IEQNS=LNUEQ(ITOTV)
       IF(IEQNS.GT.0) DISIT(ITOTV)=ELOAD(IEQNS)
      ENDDO
C
C**** COMPUTE REACTIONS (only necessary for coupled prob. with kpore=2)
C
      IF(IPASS.EQ.2) THEN
       CALL SKYRHS(CSTIF,DISIT,REFOR,LNODS,NDOFN,NELEM,NEVAB,
     .             NNODE,NPOIN,DISIM,FOREL,    0,
     .             NPREL,NGRUP,NPROP,NMATS,
     .             NDATA,NPREV,NSTAT,NMATX,NTOTV,NDIME,
     .             NDISR,NDISO,
     .             MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .             ELVAR,ELMAT,WORK1,TEMPN,DISPR,DTEMP,
     .             VNORM,DISTO,COORD,INFRI,COFRI,LACTI)
      ENDIF
C
      RETURN
      END
