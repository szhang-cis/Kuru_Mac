      SUBROUTINE FORCINS(DISTOS,ELDATS,ELPRES,ELVARS,ELMATS,ELOADS,
     .                   HEADSS,LNODSS,MATNOS,PROELS,PROPSS,WORK1S,
     .                   ELELTS,VELCTS,VEL1CS,DISPLS,PREASS,TGAPSS,
     .                   COORDS,TEMPIS,ADVELS,BOUCHS,FPCHAS)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE RESISTING FORCES
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
      COMMON/JACOBSSA/IERORS,KERORS
C
      DIMENSION MATNOS(NELEMS),        LNODSS(NNODES,NELEMS),
     .          PROELS(NPRELS,NGRUPS), PROPSS(NPROPS,NMATSS),
     .          ELDATS(NDATAS),        ELPRES(NPREVS),
     .          ELVARS(NSTATS),        ELMATS(NMATXS),
     .          WORK1S(*)
      DIMENSION DISTOS(*),             ELOADS(*),
     .          ELELTS(*),             HEADSS(*),
     .          VELCTS(*),             VEL1CS(*)
      DIMENSION DISPLS(NTOTVMS)
      DIMENSION PREASS(*),             TGAPSS(*)
      DIMENSION COORDS(NDIMES,NPOINS)
      DIMENSION TEMPIS(NPOINS),        ADVELS(NTOTVS*NDIMES),
     .          BOUCHS(NPOINS),        FPCHAS(NFPCH,NPOINS)
C
      KUNLDS=0
      IF(LARGES.NE.0) KERORS=0
C
C**** LOOP OVER ELEMENTS
C
      DO 1000 IELEMS=1,NELEMS
      LGRUPS=MATNOS(IELEMS)
      LMATSS=INT(PROELS(1,LGRUPS))
      NNODLS=INT(PROELS(2,LGRUPS))
C
      IF(LARGES.NE.0) IERORS=0
C
C**** READ ELEM. DATA FROM DATA BASE
C
      IF(KPROBS.EQ.3) THEN
       CALL DATBASS(ELMATS(IMATXS(2)),    7,    2) ! ESTIF
      ELSE
       IF(NMEMO1S.EQ.0.OR.NMEMO2S.EQ.0)
     .  CALL DATBASS(ELDATS,    1,    2)   ! geometrical data
C
       IF(NMEMOS.EQ.1)
     .  CALL DATBASS(ELPRES,    2,    2)   ! current values
C
       IF(NMEMO3S.EQ.1.OR.NMEMO4S.EQ.1.OR.NMEMO5S.EQ.0)
     .  CALL DATBASS(ELVARS,    5,    2)   ! last converged
C
       IF(NMEMO5S.EQ.0) THEN
C
C**** SCALAR GATHER OPERATIONS ( DISTOS ---> ELDISS(ELVARS) )
C
        CALL GATHER(DISTOS,NDOFCS,NPOINS,ELVARS,         NDOFCS,NNODLS,
     .              LNODSS(1,IELEMS))
       ELSE
C
C**** SCALAR GATHER OPERATIONS ( DISTOS ---> ELDISS )
C
        CALL GATHER(DISTOS,NDOFCS,NPOINS,WORK1S(IFORCS(44)),NDOFCS,
     .              NNODLS,
     .              LNODSS(1,IELEMS))
       ENDIF                   ! nmemo5s.eq.0
      END IF                   ! kprobs.eq.3
C
C**** SCALAR GATHER OPERATIONS ( VELCTS ---> WORK1S(IFORCS(12)) )
C
      IF(KDYNAS.EQ.1) THEN
       CALL GATHER(VELCTS,NDOFCS,NPOINS,WORK1S(IFORCS(12)),NDOFCS,
     .             NNODLS,
     .             LNODSS(1,IELEMS))     ! current temp. rate
      ENDIF
C
C**** SCALAR GATHER OPERATIONS ( VEL1CS ---> WORK1S(IFORCS(33)) )
C
c     IF(KDYNAS.EQ.1) THEN
c      IF(IMICR.EQ.1) THEN
c       IF(IMICO.EQ.0)                   ! weak microst. coupling
c    .   CALL GATHER(VEL1CS,NDOFCS,NPOINS,WORK1S(IFORCS(33)),NDOFCS,
c    .               NNODLS,
c    .               LNODSS(1,IELEMS))     ! predicted temp. rate
c      ENDIF
c     ENDIF
C
C**** SCALAR GATHER OPERATIONS ( TEMPIS ---> WORK1S(IFORCS(45)) )
C
      CALL GATHER(TEMPIS,NDOFCS,NPOINS,WORK1S(IFORCS(45)),NDOFCS,NNODLS,
     .            LNODSS(1,IELEMS))     ! initial temperature
C
C**** SCALAR GATHER OPERATIONS ( BOUCHT ---> WORK1T(IFORCT(46)) )
C
c     IF(NMEMO10S.EQ.1)
c    .  CALL GATHER(BOUCHT,NDOFCT,NPOINT,WORK1T(IFORCT(46)),NDOFCT,
c    .              NNODLT,
c    .              LNODST(1,IELEMT))   ! boundary changes due to rho
C
c     IF(ITERME.GT.0) THEN              ! bidirectional coupled
c      IAUXX=0
c      IF(ITERMG.GT.0) THEN             ! gap dependency
c       IF(NMEMO11.EQ.0) THEN
C
C**** SCALAR GATHER OPERATIONS ( PREAST ---> WORK1T(IFORCT(14)) )
C
c        CALL GATHER(PREAST,NDOFCT,NPOINT,WORK1T(IFORCT(14)),NDOFCT,
c    .               NNODLT,
c    .               LNODST(1,IELEMT))
C
C**** SCALAR GATHER OPERATIONS ( TGAPST ---> WORK1T(IFORCT(15)) )
C
c        CALL GATHER(TGAPST,NDOFCT,NPOINT,WORK1T(IFORCT(15)),NDOFCT,
c    .               NNODLT,
c    .               LNODST(1,IELEMT))
c       ELSE
C
C**** SCALAR GATHER OPERATIONS ( DISPLT ---> WORK1T(IFORCT(13)) )
C
c        IAUXX=1
c        CALL GATHER(DISPLT,NDOFCM,NPOINT,WORK1T(IFORCT(13)),NDOFCM,
c    .               NNODLT,
c    .               LNODST(1,IELEMT))
c       ENDIF                       ! nmemo11.eq.0
c      ENDIF                        ! itermg.gt.0
C
c      IF(ITERMD.EQ.1.AND.IAUXX.EQ.0) THEN     ! deformed shape
C
C**** SCALAR GATHER OPERATIONS ( DISPLT ---> WORK1T(IFORCT(13)) )
C
c       CALL GATHER(DISPLT,NDOFCM,NPOINT,WORK1T(IFORCT(13)),NDOFCM,
c    .              NNODLT,
c    .              LNODST(1,IELEMT))
c      ENDIF                        ! itermd.eq.1.and.iauxx.eq.0
c     ENDIF                         ! iterme.gt.0
C
C**** SCALAR GATHER OPERATIONS ( ADVELT ---> WORK1T(IFORCT(34)) )
C
      IF(ICONVS.EQ.1)
     . CALL GATHER(ADVELS,NDIMES,NTOTVS,WORK1S(IFORCS(34)),NDIMES,
     .             NNODLS,
     .             LNODSS(1,IELEMS))
C
C**** GATHER NODAL COORDINATES
C
      IF(NMEMO1S.EQ.1) THEN    ! coordinates in a global array
       CALL GATHER(COORDS,NDIMES,NPOINS,WORK1S(IFORCS(16)),NDIMES,
     .             NNODLS,LNODSS(1,IELEMS))
      ENDIF
C
C**** GATHER PHASE-CHANGE VECTOR
C
      IF(NFILLS.EQ.1)
     . CALL GATHER(FPCHAS,NFPCH,NPOINS,WORK1S(IFORCS(48)),NFPCH,
     .             NNODLS,LNODSS(1,IELEMS))
C
C**** CALL ELEMENT PROCESSOR TO PERFORM REQUIRED OPERATIONS
C
      CALL ELMLIBS(LNODSS(1,IELEMS),PROELS(1,LGRUPS),PROPSS(1,LMATSS),
     .             ELDATS,ELPRES,ELVARS,ELMATS,WORK1S,       4)
C
      IF(LARGES.NE.0) THEN
       IF(IERORS.NE.0) GO TO 1000
      ENDIF
C
C**** WRITE CURRENT STATE VARIABLES TO DATA BASE
C 
      IF(NMEMO3S.EQ.1.OR.NMEMO4S.EQ.1.OR.NMEMO5S.EQ.0)
     . CALL DATBASS(ELVARS,    3,    1)
C
C**** SCALAR SCATTER OPERATION ( WORK1S --> ELOADS )
C
      CALL SCATER(WORK1S(IFORCS(1)), NDOFNS,NNODLS,ELOADS,NDOFCS,NPOINS,
     .            LNODSS(1,IELEMS))
C
C**** SCALAR SCATTER OPERATION ( WORK1S --> ELELTS )
C
      CALL SCATER(WORK1S(IFORCS(11)),NDOFNS,NNODLS,ELELTS,NDOFCS,NPOINS,
     .            LNODSS(1,IELEMS))
C
C**** SCALAR SCATTER OPERATION ( WORK1S --> BOUCHS )
C
      IF(NMEMO10S.EQ.1)
     . CALL SCATERA(WORK1S(IFORCS(46)),NDOFCS,NNODLS,BOUCHS,NDOFCS,
     .              NPOINS,LNODSS(1,IELEMS))
C
 1000 CONTINUE
C
      IF(LARGES.NE.0) THEN
       IF(KERORS.NE.0)
     .  CALL RUNENDS('ERROR IN JACOBIAN MATRIX           ')
      ENDIF
C
      RETURN
      END
