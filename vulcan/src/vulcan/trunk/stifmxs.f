      SUBROUTINE STIFMXS(ELDATS,ELPRES,ELVARS,ELMATS,LNODSS,MATNOS,
     .                   PROELS,PROPSS,WORK1S,VELCTS,DISITS,COORDS,
     .                   ADVELS,TEMPIS,PREASS,TGAPSS,DISTOS,BOUCHS,
     .                   DISPLS,FPCHAS)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE STIFFNESS & COUPLING MATRICES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_oms.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'          ! thermal-mechanical
      INCLUDE 'nued_om.f'          ! thermal-microstructural
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      DIMENSION MATNOS(NELEMS),        LNODSS(NNODES,NELEMS),
     .          PROELS(NPRELS,NGRUPS), PROPSS(NPROPS,NMATSS),
     .          ELDATS(NDATAS),        ELPRES(NPREVS),
     .          ELVARS(NSTATS),        ELMATS(NMATXS),
     .          WORK1S(*),             VELCTS(*),
     .          DISITS(*),             COORDS(NDIMES,NPOINS)
      DIMENSION ADVELS(NTOTVS*NDIMES), TEMPIS(NPOINS)
      DIMENSION PREASS(NPOINS),        TGAPSS(NPOINS)
      DIMENSION DISTOS(*)
      DIMENSION BOUCHS(NPOINS),        DISPLS(NTOTVMS),
     .          FPCHAS(NFPCH,NPOINS)
C
      CALL CPUTIMS(TIME1S)
C
C**** LOOP ON ELEMENTS
C
      DO 10 IELEMS=1,NELEMS
      LGRUPS=MATNOS(IELEMS)
      LMATSS=INT(PROELS(1,LGRUPS))
      NNODLS=INT(PROELS(2,LGRUPS))
C
C**** READ ELDAT & ELVAR FROM DATA BASE
C
      IF(NMEMO1S.EQ.0.OR.NMEMO2S.EQ.0)
     . CALL DATBASS(ELDATS,    1,    2)
C
      IF(NMEMO3S.EQ.1.OR.NMEMO4S.EQ.1.OR.NMEMO5S.EQ.0)
     . CALL DATBASS(ELVARS,    3,    2)    ! current
C
      IF(NMEMO5S.EQ.1) THEN
C
C**** SCALAR GATHER OPERATIONS ( DISTOS ---> ELDISS )
C
       CALL GATHER(DISTOS,NDOFCS,NPOINS,WORK1S(ISTIFS(39)),NDOFCS,
     .             NNODLS,
     .             LNODSS(1,IELEMS))
      ENDIF
C
C**** SCALAR GATHER OPERATIONS
C
      IF(KDYNAS.EQ.1) THEN
       CALL GATHER(VELCTS,NDOFCS,NPOINS,WORK1S(ISTIFS(13)),NDOFCS,
     .             NNODLS,
     .             LNODSS(1,IELEMS))
      ENDIF
C
C**** SCALAR GATHER OPERATIONS
C
      CALL GATHER(DISITS,NDOFCS,NPOINS,WORK1S(ISTIFS(17)),NDOFCS,NNODLS,
     .            LNODSS(1,IELEMS))
C
C**** SCALAR GATHER OPERATIONS ( ADVELS: advective velocity)
C
      IF(ICONVS.EQ.1)
     . CALL GATHER(ADVELS,NDIMES,NTOTVS,WORK1S(ISTIFS(25)),NDIMES,
     .             NNODLS,
     .             LNODSS(1,IELEMS))
C
C**** GATHER NODAL COORDINATES
C
      IF(NMEMO1S.EQ.1) THEN    ! coordinates in a global array
       CALL GATHER(COORDS,NDIMES,NPOINS,WORK1S(ISTIFS(18)),NDIMES,
     .             NNODLS,LNODSS(1,IELEMS))
      ENDIF
C
C**** SCALAR GATHER OPERATIONS ( TEMPIS ---> WORK1S(ISTIFS(40)) )
C
      CALL GATHER(TEMPIS,NDOFCS,NPOINS,WORK1S(ISTIFS(40)),NDOFCS,NNODLS,
     .            LNODSS(1,IELEMS))     ! initial temperature
C
C**** SCALAR GATHER OPERATIONS ( BOUCHS ---> WORK1S(ISTIFS(42)) )
C
      IF(NMEMO10S.EQ.1)
     . CALL GATHER(BOUCHS,NDOFCS,NPOINS,WORK1S(ISTIFS(42)),NDOFCS,
     .             NNODLS,
     .             LNODSS(1,IELEMS))       ! boundary changes due to rho
C
      IF(ITERME.GT.0) THEN                       ! bidirectional coupled
       IAUXX=0
       IF(ITERMG.GT.0) THEN                      ! gap dependency
        IF(NMEMO11S.EQ.0) THEN
C
C**** SCALAR GATHER OPERATIONS ( PREASS ---> WORK1S(ISTIFS(35)) )
C
         CALL GATHER(PREASS,NDOFCS,NPOINS,WORK1S(ISTIFS(35)),NDOFCS,
     .               NNODLS,LNODSS(1,IELEMS))
C
C**** SCALAR GATHER OPERATIONS ( TGAPSS ---> WORK1S(ISTIFS(36)) )
C
         CALL GATHER(TGAPSS,NDOFCS,NPOINS,WORK1S(ISTIFS(36)),NDOFCS,
     .               NNODLS,LNODSS(1,IELEMS))
        ELSE
C
C**** SCALAR GATHER OPERATIONS ( DISPLS ---> WORK1S(ISTIFS(44)) )
C
         IAUXX=1
         CALL GATHER(DISPLS,NDOFCMS,NPOINS,WORK1S(ISTIFS(44)),NDOFCMS,
     .               NNODLS,
     .               LNODSS(1,IELEMS))
        ENDIF                        ! nmemo11s.eq.0
       ENDIF                         ! itermg.gt.0
C
       IF(ITERMD.EQ.1.AND.IAUXX.EQ.0) THEN     ! deformed shape
C
C**** SCALAR GATHER OPERATIONS ( DISPLS ---> WORK1S(ISTIFS(44)) )
C
        CALL GATHER(DISPLS,NDOFCMS,NPOINS,WORK1S(ISTIFS(44)),NDOFCMS,
     .              NNODLS,
     .              LNODSS(1,IELEMS))
       ENDIF                        ! itermd.eq.1.and.iauxx.eq.0
      ENDIF                          ! iterme.gt.0
C
C**** GATHER PHASE-CHANGE VECTOR
C
      IF(NFILLS.EQ.1)
     . CALL GATHER(FPCHAS,NFPCH,NPOINS,WORK1S(ISTIFS(45)),NFPCH,
     .             NNODLS,LNODSS(1,IELEMS))
C
C**** CALL ELEMENT PROCESSOR TO PERFORM REQUIRED OPERATIONS
C
      CALL ELMLIBS(LNODSS(1,IELEMS),PROELS(1,LGRUPS),PROPSS(1,LMATSS),
     .             ELDATS,ELPRES,ELVARS,ELMATS,WORK1S,       3)
C
C**** WRITE ELMAT TO DATA BASE
C
      IF(NMEMO6S.EQ.0) THEN
       CALL DATBASS(ELMATS(IMATXS(2)),    7,    1)   ! K
       IF(KDYNAS.EQ.1)
     .  CALL DATBASS(ELMATS(IMATXS(3)),    8,    1)  ! M
      ELSE
       CALL DATBASS(ELMATS(IMATXS(1)),    6,    1)   ! C
      ENDIF
C
   10 CONTINUE
C
      CALL CPUTIMS(TIME2S)
      CPUSFS=CPUSFS+(TIME2S-TIME1S)
C
      RETURN
      END
