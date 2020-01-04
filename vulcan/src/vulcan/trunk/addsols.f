      SUBROUTINE ADDSOLS(NEQNSS,NFRONS,NLASTS,NSTIFS,WORK1S)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE MEMORY REQUIREMENTS & ADDRESS FOR 
C     SOLVER
C
C***********************************************************************
C
C**** GENERAL DIMENSIONS
C
      INCLUDE 'para_oms.f'
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_oms.f'
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
    
      REAL*8, ALLOCATABLE, INTENT(INOUT) :: WORK1S(:)
      REAL*8, ALLOCATABLE                :: TEMP(:)
      INTEGER                            :: NSIZEWORK1S,ERROR

      NSIZEWORK1S = SIZE(WORK1S)
C
      ISOLAS=1
      IF(NMEMO7S.EQ.1) ISOLAS=1+LSTIFS
C
C**** SKYLINE SOLVER
C
      IF(KSOLVS.EQ.0) THEN
       IUNSYS=0
       IF(KSYMMS.EQ.0) IUNSYS=1
       IFITES=0
       IF(MITCGS.GT.0) IFITES=1
       ISOLVS( 1)=ISOLAS                           ! GSTDI(NEQNS)
       ISOLVS( 2)=ISOLVS( 1)+NEQNSS                ! GSTLO(NLAST*IUNSY)
       ISOLVS( 3)=ISOLVS( 2)+NLASTS*IUNSYS         ! GSTUP(NLAST)
       ISOLVS( 4)=ISOLVS( 3)+NLASTS                ! CSTIF(NEVAC,NEVAC)
       ISOLVS( 5)=ISOLVS( 4)+NEVACS*NEVACS         ! ELOAD(NEQNS)
       ISOLVS( 6)=ISOLVS( 5)+NEQNSS                ! LNUEQ(NTOTV)
       ISOLVS( 7)=ISOLVS( 6)+(NTOTVS*ICHALS+4)/8   ! LPONT(NTOTV)
       ISOLVS( 8)=ISOLVS( 7)+(NTOTVS*ICHALS+4)/8   ! DISIM(NEVAC)
       ISOLVS( 9)=ISOLVS( 8)+NEVACS                ! FOREL(NEVAC)
       ISOLVS(10)=ISOLVS( 9)+NEVACS                ! ALOAD(NEQNS*IFITE)
       ISOLVS(11)=ISOLVS(10)+NEQNSS*IFITES         ! DELTA(NEQNS*IFITE)
       ISOLVS(12)=ISOLVS(11)+NEQNSS*IFITES         ! DISIC(NEQNS*IFITE)
       ISOLVS(13)=ISOLVS(12)+NEQNSS*IFITES         ! LOCAL(NEVAC*IFITE)
       LSOLVS    =ISOLVS(13)+(NEVACS*IFITES*ICHALS+4)/8
      ENDIF
C
C**** FRONTAL SOLVER
C
      IF(KSOLVS.EQ.1) THEN
       IUNSYS=0
       IF(KSYMMS.EQ.0) IUNSYS=1
       ISOLVS(1) =ISOLAS                           ! EQRHS(NBUFA)
       ISOLVS(2) =ISOLVS(1)+NBUFAS                 ! EQUAT(NFRON,NBUFA)
       ISOLVS(3) =ISOLVS(2)+NFRONS*NBUFAS          ! GLOAD(NFRON)
       ISOLVS(4) =ISOLVS(3)+NFRONS                 ! GSTIF(NSTIF)
       ISOLVS(5) =ISOLVS(4)+NSTIFS                 ! LOCEL(NEVAC)
       ISOLVS(6) =ISOLVS(5)+(NEVACS*ICHALS+4)/8    ! NACVA(NFRON)
       ISOLVS(7) =ISOLVS(6)+(NFRONS*ICHALS+4)/8    ! NAMEV(NBUFA)
       ISOLVS(8) =ISOLVS(7)+(NBUFAS*ICHALS+4)/8    ! NDEST(NEVAC)
       ISOLVS(9) =ISOLVS(8)+(NEVACS*ICHALS+4)/8    ! NPIVO(NBUFA)
       ISOLVS(10)=ISOLVS(9)+(NBUFAS*ICHALS+4)/8    ! VECRV(NFRON)
       ISOLVS(11)=ISOLVS(10)+NFRONS                ! CSTIF(NEVAC,NEVAC)
       ISOLVS(12)=ISOLVS(11)+NEVACS*NEVACS   ! EQCOL(IUNSY*NFRON*NBUFA)
       LSOLVS    =ISOLVS(12)+IUNSYS*NFRONS*NBUFAS
      ENDIF
C
C**** PCG SOLVER
C
      IF(KSOLVS.EQ.2) THEN
       NSIZES=NDOFCS
       IF(NWIDTS.LT.0) NSIZES=NDOFCS*NDOFCS
       ISOLVS( 1)=ISOLAS                           ! GSTDI(NPOIN*NSIZE)
       ISOLVS( 2)=ISOLVS( 1)+NPOINS*NSIZES         ! CSTIF(NEVAC,NEVAC)
       ISOLVS( 3)=ISOLVS( 2)+NEVACS*NEVACS         ! DISIM(NEVAC)
       ISOLVS( 4)=ISOLVS( 3)+NEVACS                ! FOREL(NEVAC)
       ISOLVS( 5)=ISOLVS( 4)+NEVACS                ! ALOAD(NTOTV)
       ISOLVS( 6)=ISOLVS( 5)+NTOTVS                ! DELTA(NTOTV)
       LSOLVS    =ISOLVS( 6)+NTOTVS
      ENDIF
C
C**** GMRES SOLVER (not implemented yet)
C

C
C**** EXPLICIT SOLUTION (NO SOLVER)
C
      IF(KSOLVS.EQ.4) THEN
       ISOLVS( 1)=ISOLAS                           ! CSTIF(NEVAC,NEVAC)
       ISOLVS( 2)=ISOLVS( 1)+NEVACS*NEVACS         ! DISIM(NEVAC)
       ISOLVS( 3)=ISOLVS( 2)+NEVACS                ! FOREL(NEVAC)
       ISOLVS( 4)=ISOLVS( 3)+NEVACS                ! CSTII(NEVAC)
       ISOLVS( 5)=ISOLVS( 4)+NEVACS                ! CSTIT(NTOTV)
       LSOLVS    =ISOLVS( 5)+NTOTVS
      ENDIF
C
C**** MEMORY CONTROL
C
      IF(LSOLVS.GT.NSIZEWORK1S) THEN
        WRITE(LURESS,*) "IN ADDSOLS REALLOCATE WORK1S FROM: ", 
     .                   NSIZEWORK1S, "TO: ", LSOLVS
        ERROR=0
        ALLOCATE(TEMP(LSOLVS),STAT=ERROR)
        IF (ERROR.NE.0)CALL RUNEND("CANT REALLOCATE WORK1S IN ADDSOLS")
        TEMP(1:NSIZEWORK1S)=WORK1S(:)
        CALL MOVE_ALLOC(FROM=TEMP,TO=WORK1S)
      END IF

C
C**** DEFINE STARTING POSITION FOR SOLVER MEMORY
C
      LBYTSS=LSOLVS*8
      IRELES=1
C
      RETURN
      END