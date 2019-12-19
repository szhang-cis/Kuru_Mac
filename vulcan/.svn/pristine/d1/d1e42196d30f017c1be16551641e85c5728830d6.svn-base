      SUBROUTINE STARTSS(DISTOS,DISPRS,ELDATS,ELPRES,ELVARS,ELMATS,
     .                   HEADSS,HTLODS,IFFIXS,LNODSS,LPNTNS,MATNOS,
     .                   PROELS,PROPSS,REFORS,RLOADS,TLOADS,COORDS,
     .                   DISPLS,PWORKS,PREASS,TGAPSS,TEMPIS,ADVELS,
     .                   FPCHAS,LACTIS,WORK1S,INDEXS)
C***********************************************************************
C
C**** THIS ROUTINE INITIALISES NODAL AND GAUSSIAN VARIABLES
C
C     Note: NFURES=2 (future restart) & IRESTT=1 (restart) have to be
C           improved
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
      DIMENSION MATNOS(NELEMS),         LNODSS(NNODES,NELEMS),
     .          PROELS(NPRELS,NGRUPS),  PROPSS(NPROPS,NMATSS),
     .          ELDATS(NDATAS),         ELPRES(NPREVS),
     .          ELVARS(NSTATS),         ELMATS(NMATXS),
     .          COORDS(NDIMES,NPOINS),  WORK1S(*)
      DIMENSION DISTOS(NTOTVS,*),       DISPRS(NTOTVS,*), 
     .          HEADSS(NPOINS,*),       HTLODS(NHLODS,NSUBFS,*), 
     .          IFFIXS(NTOTVS,*),       LPNTNS(*),
     .          REFORS(*),              RLOADS(*), 
     .          TLOADS(NTOTVS,*)
      DIMENSION DISPLS(NTOTVMS),        PWORKS(NPOINS,3),
     .          PREASS(NPREAS,NPOINS),  TGAPSS(NPOINS),
     .          TEMPIS(NPOINS,2),       ADVELS(NTOTVS*NDIMES),
     .          FPCHAS(NFPCH,NPOINS),   LACTIS(NELEMS)
C
      IF(IRESTS.EQ.0) THEN
C***********************************************************************
C
C**** THIS IS A FRESH RUN
C
C***********************************************************************
       IF(INDEXS.EQ.1) THEN
C
C**** READ IN MOST OF THE PROBLEM DATA
C
        CALL INPDATS(ELDATS,ELPRES,ELVARS,ELMATS,LNODSS,MATNOS,PROELS,
     .               PROPSS,COORDS,WORK1S)
C
        IF(IEVFI.EQ.0) RETURN
C
C**** RENUMBER NODES
C
        CALL RENUMNS(KRENUS,LNODSS,LPNTNS,NELEMS,NNODES,NPOINS,WORK1S)  
C
C**** WRITE TO DATA BASE
C
        IF(NFURESS.EQ.2)
     .   CALL RSTAR2S(LPNTNS,LNODSS,MATNOS,PROELS,PROPSS,    1)
C
C**** INITIALIZE GLOBAL AND ELEMENTAL VARIABLES
C
        CALL ZEROTES(ELPRES,ELVARS,DISTOS,HEADSS,TLOADS,RLOADS,
     .               REFORS,DISPLS,PWORKS,PREASS,TGAPSS,TEMPIS,
     .               ADVELS,FPCHAS)
C
C**** READ IN INITIAL GLOBAL AND ELEMENTAL VARIABLES
C
        IF(INITIS.EQ.1)
     .   CALL PREVOSS(DISTOS,ELDATS,ELPRES,ELVARS,ELMATS,HEADSS,LNODSS,
     .                MATNOS,PROELS,PROPSS,TEMPIS(1,1),WORK1S)
        IF(INITIS.EQ.2)
     .   CALL RSTAR5S(TEMPIS,DISTOS,ELPRES,ELVARS)
       ELSE
C
C**** CALCULATE THE CONSTANT MATRICES AND ARRAYS OF THE PROBLEM
C
        CALL SETMTXS(ELDATS,ELPRES,ELVARS,ELMATS,LNODSS,MATNOS,PROELS,
     .               PROPSS,COORDS,TEMPIS(1,1),FPCHAS,DISPLS,
     .               WORK1S)
       ENDIF                  ! indexs.eq.1
C
      ELSE
C
C***********************************************************************
C
C**** THIS IS A RESTART RUN
C
C***********************************************************************
       IF(INDEXS.EQ.1) THEN
C
        IF(IEVFI.EQ.0) RETURN
C
C**** READ GENERAL DATA FROM DATA BASE
C
        CALL RSTAR2S(LPNTNS,LNODSS,MATNOS,PROELS,PROPSS,    2)
C
C**** READ LAST CONVERGED VALUES FROM DATA BASE
C
        CALL RSTAR3S(DISPRS,DISTOS,ELPRES,ELVARS,HEADSS,HTLODS,IFFIXS,
     .               REFORS,RLOADS,TLOADS,    2)
       ENDIF                      ! indexs.eq.1
C
      ENDIF                       ! irests.eq.0
C
C***********************************************************************
C
C**** DUMP GPCOD TO TAPE FOR POSTPROCESS IF NECESSARY
C
      IF(INDEXS.EQ.1) RETURN
C
      IF(NWPOSS.EQ.1) THEN
       IF(NMEMO1S.EQ.0) THEN
        CALL OUTPOSS(ELDATS,LNODSS,MATNOS,PROELS,PROPSS,
     .               WORK1S(ISTARS(1)),WORK1S(ISTARS(2)),
     .               ELPRES,ELVARS,ELMATS,DISPLS,WORK1S)
       ELSE
        CALL OUTPOSS(ELDATS,LNODSS,MATNOS,PROELS,PROPSS,
     .               COORDS,WORK1S(ISTARS(1)),
     .               ELPRES,ELVARS,ELMATS,DISPLS,WORK1S)
       ENDIF
      ENDIF
C
      RETURN
      END
