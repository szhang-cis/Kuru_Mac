      SUBROUTINE STIFMXT(ELDATT,ELPRET,ELVART,ELMATT,LNODST,MATNOT,
     .                   PROELT,PROPST,WORK1T,VELCTT,DISITT,COORDT,
     .                   ADVELT,TEMPIT,PREAST,TGAPST,DISTOT,BOUCHT,
     .                   DISPLT,FPCHAT,LACTIT)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE STIFFNESS & COUPLING MATRICES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'          ! thermal-mechanical
      INCLUDE 'nued_om.f'          ! thermal-microstructural
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION MATNOT(NELEMT),        LNODST(NNODET,NELEMT),
     .          PROELT(NPRELT,NGRUPT), PROPST(NPROPT,NMATST),
     .          ELDATT(NDATAT),        ELPRET(NPREVT),
     .          ELVART(NSTATT),        ELMATT(NMATXT),
     .          WORK1T(*),             VELCTT(*),
     .          DISITT(*),             COORDT(NDIMET,NPOINT)
      DIMENSION ADVELT(NTOTVT*NDIMET), TEMPIT(NPOINT)
      DIMENSION PREAST(NPOINT),        TGAPST(NPOINT)
      DIMENSION DISTOT(*)
      DIMENSION BOUCHT(NPOINT),        DISPLT(NTOTVM),
     .          FPCHAT(NFPCH,NPOINT),  LACTIT(NELEMT)
C
      CALL CPUTIMT(TIME1T)
C
C**** LOOP ON ELEMENTS
C
      DO 1000 IELEMT=1,NELEMT
      IF(NACTIT.EQ.1) THEN
       IF(LACTIT(IELEMT).EQ.0) THEN     ! deals with non active elements
        IF(NMEMO6.EQ.0) THEN
         N=IMATXT(3)-IMATXT(2)
         CALL VECZER(N,ELMATT(IMATXT(2)))              ! K
         IF(KDYNAT.EQ.1) THEN
          N=IMATXT(4)-IMATXT(3)
          CALL VECZER(N,ELMATT(IMATXT(3)))             ! M
         ENDIF
        ELSE
         N=IMATXT(2)-IMATXT(1)
         CALL VECZER(N,ELMATT(IMATXT(1)))              ! C
        ENDIF
        GO TO 1001
       ENDIF
      ENDIF
      LGRUPT=MATNOT(IELEMT)
      LMATST=INT(PROELT(1,LGRUPT))
      NNODLT=INT(PROELT(2,LGRUPT))
      LTYPET=INT(PROELT(5,LGRUPT))
C
      NNODXT=NNODLT
      IF(NOCOIT.EQ.1) THEN                    ! non-coincident mesh
       IF(LTYPET.EQ.104) THEN
        NOCOLT=INT(PROELT(11,LGRUPT))
        IF(NOCOLT.EQ.1) THEN
         NNODNT=INT(PROELT(10,LGRUPT))
         NNODXT=NNODNT
        ENDIF
       ENDIF
      ENDIF
C
C**** READ ELDAT & ELVAR FROM DATA BASE
C
      IF(NMEMO1.EQ.0.OR.NMEMO2.EQ.0)
     . CALL DATBAST(ELDATT,    1,    2)
C
      IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
     . CALL DATBAST(ELVART,    3,    2)    ! current
C
      IF(NMEMO5.EQ.1) THEN
C
C**** SCALAR GATHER OPERATIONS ( DISTOT ---> ELDIST )
C
       CALL GATHER(DISTOT,NDOFCT,NPOINT,WORK1T(ISTIFT(39)),NDOFCT,
     .             NNODXT,
     .             LNODST(1,IELEMT))
      ENDIF
C
C**** SCALAR GATHER OPERATIONS
C
      IF(KDYNAT.EQ.1) THEN
       CALL GATHER(VELCTT,NDOFCT,NPOINT,WORK1T(ISTIFT(13)),NDOFCT,
     .             NNODXT,
     .             LNODST(1,IELEMT))
      ENDIF
C
C**** SCALAR GATHER OPERATIONS
C
      CALL GATHER(DISITT,NDOFCT,NPOINT,WORK1T(ISTIFT(17)),NDOFCT,NNODLT,
     .            LNODST(1,IELEMT))
C
C**** SCALAR GATHER OPERATIONS ( ADVELT: advective velocity)
C
      IF(ICONVT.EQ.1)
     . CALL GATHER(ADVELT,NDIMET,NTOTVT,WORK1T(ISTIFT(25)),NDIMET,
     .             NNODLT,
     .             LNODST(1,IELEMT))
C
C**** GATHER NODAL COORDINATES
C
      IF(NMEMO1.EQ.1) THEN    ! coordinates in a global array
       CALL GATHER(COORDT,NDIMET,NPOINT,WORK1T(ISTIFT(18)),NDIMET,
     .             NNODLT,LNODST(1,IELEMT))
      ENDIF
C
C**** SCALAR GATHER OPERATIONS ( TEMPIT ---> WORK1T(ISTIFT(40)) )
C
      CALL GATHER(TEMPIT,NDOFCT,NPOINT,WORK1T(ISTIFT(40)),NDOFCT,NNODLT,
     .            LNODST(1,IELEMT))     ! initial temperature
C
C**** SCALAR GATHER OPERATIONS ( BOUCHT ---> WORK1T(ISTIFT(42)) )
C
      IF(NMEMO10.EQ.1)
     . CALL GATHER(BOUCHT,NDOFCT,NPOINT,WORK1T(ISTIFT(42)),NDOFCT,
     .             NNODLT,
     .             LNODST(1,IELEMT))       ! boundary changes due to rho
C
      IF(ITERME.GT.0) THEN                       ! bidirectional coupled
       IAUXX=0
       IF(ITERMG.GT.0) THEN                      ! gap dependency
        IF(NMEMO11.EQ.0) THEN
C
C**** SCALAR GATHER OPERATIONS ( PREAST ---> WORK1T(ISTIFT(35)) )
C
         CALL GATHER(PREAST,NDOFCT,NPOINT,WORK1T(ISTIFT(35)),NDOFCT,
     .               NNODLT,LNODST(1,IELEMT))
C
C**** SCALAR GATHER OPERATIONS ( TGAPST ---> WORK1T(ISTIFT(36)) )
C
         CALL GATHER(TGAPST,NDOFCT,NPOINT,WORK1T(ISTIFT(36)),NDOFCT,
     .               NNODLT,LNODST(1,IELEMT))
        ELSE
C
C**** SCALAR GATHER OPERATIONS ( DISPLT ---> WORK1T(ISTIFT(44)) )
C
         IAUXX=1
         CALL GATHER(DISPLT,NDOFCM,NPOINT,WORK1T(ISTIFT(44)),NDOFCM,
     .               NNODLT,
     .               LNODST(1,IELEMT))
        ENDIF                        ! nmemo11.eq.0
       ENDIF                         ! itermg.gt.0
C
       IF((ITERMD.EQ.1.OR.ITERMF.EQ.1).AND.
     .                             IAUXX.EQ.0) THEN     ! deformed shape
C
C**** SCALAR GATHER OPERATIONS ( DISPLT ---> WORK1T(ISTIFT(44)) )
C
        CALL GATHER(DISPLT,NDOFCM,NPOINT,WORK1T(ISTIFT(44)),NDOFCM,
     .              NNODLT,
     .              LNODST(1,IELEMT))
       ENDIF                         ! itermd.eq.1. ... .and.iauxx.eq.0
      ENDIF                          ! iterme.gt.0
C
C**** GATHER PHASE-CHANGE VECTOR
C
      IF(NFILL.EQ.1)
     . CALL GATHER(FPCHAT,NFPCH,NPOINT,WORK1T(ISTIFT(45)),NFPCH,
     .             NNODLT,LNODST(1,IELEMT))
C
C**** CALL ELEMENT PROCESSOR TO PERFORM REQUIRED OPERATIONS
C
      CALL ELMLIBT(LNODST(1,IELEMT),PROELT(1,LGRUPT),PROPST(1,LMATST),
     .             ELDATT,ELPRET,ELVART,ELMATT,DUMMYT,WORK1T,     3)
C
C**** WRITE ELMAT TO DATA BASE
C
 1001 IF(NMEMO6.EQ.0) THEN
       CALL DATBAST(ELMATT(IMATXT(2)),    7,    1)   ! K
       IF(KDYNAT.EQ.1)
     .  CALL DATBAST(ELMATT(IMATXT(3)),    8,    1)  ! M
      ELSE
       CALL DATBAST(ELMATT(IMATXT(1)),    6,    1)   ! C
      ENDIF
C
 1000 CONTINUE
C
      CALL CPUTIMT(TIME2T)
      CPUSFT=CPUSFT+(TIME2T-TIME1T)
C
      RETURN
      END
