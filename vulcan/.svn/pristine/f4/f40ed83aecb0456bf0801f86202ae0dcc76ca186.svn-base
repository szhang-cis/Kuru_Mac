      SUBROUTINE SETMTXS(ELDATS,ELPRES,ELVARS,ELMATS,LNODSS,MATNOS,
     .                   PROELS,PROPSS,COORDS,TEMPIS,FPCHAS,DISPLS,
     .                   WORK1S)
C***********************************************************************
C
C**** THIS ROUTINE SETS UP SOME CONSTANT MATRICES FOR FUTURE USE:
C
C         - ROTATION MATRIX             (RMAT1)
C         - ELASTIC CONSTANTS MATRIX    (EPMTX)
C         - NODAL SMOOTHING MATRICES    (SMATR;RMATR,FMATR)
C         - GAUSSIAN SMOOTHING MATRICES (EMASS)
C
C     IT ALSO EVALUATES ONCE AND FOR ALL:
C
C         - SHAPE FUNCTIONS          (SHAPE)
C         - CARTESIAN DERIVATIVES    (CARTD)
C         - DIFFERENTIAL VOLUMEN     (DVOLU)
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
      INCLUDE 'nuec_om.f'
      INCLUDE 'nued_om.f'
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
     .          COORDS(NDIMES,NPOINS), TEMPIS(NPOINS),
     .          FPCHAS(NFPCH,NPOINS),  DISPLS(NTOTVMS),
     .          WORK1S(*)
C
      CALL CPUTIMS(TIME1S)
C
C**** LOOP ON ELEMENTS
C
      KERORS=0
      DO 1000 IELEMS=1,NELEMS
      LGRUPS=MATNOS(IELEMS)
      LMATSS=INT(PROELS(1,LGRUPS))
      NNODLS=INT(PROELS(2,LGRUPS))
      LTYPES=INT(PROELS(5,LGRUPS))
C
      IERORS=0
C
      IF(NMEMO1S.EQ.0) THEN     ! coordinates in an elemental array
C
C**** READ ELDAT ( coordinates ) FROM DATA-BASE
C
       CALL DATBASS(ELDATS,    1,    2)
C
      ELSE                     ! coordinates in a global array
       CALL GATHER(COORDS,NDIMES,NPOINS,WORK1S(ISETMS(6)),NDIMES,
     .             NNODLS,LNODSS(1,IELEMS))
      ENDIF
C
C**** READ ELVAR FOR MICROSTRUCTURAL PROBLEMS
C
c     IF(IMICR.EQ.1) THEN             ! to be improved
c      CALL DATBASS(ELVARS,    3,    2)
c     ENDIF
C
C**** SCALAR GATHER OPERATIONS ( TEMPIT ---> WORK1T(ISETMT(13)) )
C
      CALL GATHER(TEMPIS,NDOFCS,NPOINS,WORK1S(ISETMS(13)),NDOFCS,NNODLS,
     .            LNODSS(1,IELEMS))     ! initial temperature
C
      IF(ITERME.GT.0) THEN          ! bidirectional coupled
       IF(ITERMD.EQ.1) THEN         ! deformed shape
C
C**** SCALAR GATHER OPERATIONS ( DISPLT ---> WORK1T(ISETMT(18)) )
C
        CALL GATHER(DISPLS,NDOFCMS,NPOINS,WORK1S(ISETMS(18)),NDOFCMS,
     .              NNODLS,
     .              LNODSS(1,IELEMS))
       ENDIF                        ! itermd.eq.1
      ENDIF                         ! iterme.gt.0
C
C**** SET UP ARRAY ELDAT & WRITE IT TO DATA BASE
C
      CALL ELMLIBS(LNODSS(1,IELEMS),PROELS(1,LGRUPS),PROPSS(1,LMATSS),
     .             ELDATS,ELPRES,ELVARS,ELMATS,WORK1S,       1)
C
      IF(IERORS.NE.0) GO TO 1000
C
C**** WRITE ELDAT TO DATA-BASE
C
      IF(NMEMO1S.EQ.0.OR.NMEMO2S.EQ.0)
     . CALL DATBASS(ELDATS,    1,    1)
C
C**** WRITE ELVAR FOR MICROSTRUCTURAL PROBLEMS (initial conditions)
C
C     ELVAR will be considered in adddatt.f as last converged values
C
C
C     ELVAR IS NOT WRITTEN FOR RESTART PROBLEMS (the initial conditions
C     are in *.pan file)
C
c     IF(INITIS.EQ.0) THEN
c      IF(IMICR.EQ.1)
c    .  CALL DATBASS(ELVARS,    3,    1)
c     ENDIF
C
      IF(LTYPES.EQ.5)
     . CALL SCATERA(WORK1S(ISETMS(14)),NFPCH,NNODLS,FPCHAS,NFPCH,
     .              NPOINS,LNODSS(1,IELEMS))
C
 1000 CONTINUE
C
      IF(KERORS.NE.0)
     . CALL RUNENDS('ERROR IN JACOBIAN MATRIX           ')
C
      CALL CPUTIMS(TIME2S)
      CPUSTS=CPUSTS+(TIME2S-TIME1S)
C
      RETURN
      END
