      SUBROUTINE TMSTEPS(DISITS,DISPRS,DISTOS,ELDATS,ELPRES,ELVARS,
     .                   ELMATS,HEADSS,HTLODS,IFFIXS,PRESCS,LNODSS,
     .                   MATNOS,PROELS,PROPSS,REFORS,RLOADS,RLOAHS,
     .                   FICTOS,TFICTS,TLOADS,DISPLS,PWORKS,PREASS,
     .                   TGAPSS,COORDS,TEMPIS,ADVELS,FPCHAS,
     .                   LPNTNS,WORK1S)
C***********************************************************************
C
C**** THIS ROUTINE ADVANCES IN TIME, UPDATE TIME AND PRESCR. VARIABLES,
C     AND PERFORMS THE PREDICTOR STAGE
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
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
      DIMENSION DISITS(NTOTVS),               DISPRS(NTOTVS,3), 
     .          DISTOS(NTOTVS,3),             ELDATS(NDATAS),      
     .          ELPRES(NPREVS),               ELVARS(NSTATS),          
     .          ELMATS(NMATXS),               HEADSS(NPOINS,4),    
     .          HTLODS(NHLODS,NSUBFS,NFUNCS), IFFIXS(NTOTVS,2),    
     .          LNODSS(NNODES,NELEMS),        MATNOS(NELEMS),
     .          PROELS(NPRELS,NGRUPS),        PROPSS(NPROPS,NMATSS),
     .          REFORS(NTOTVS),               RLOADS(NTOTVS),      
     .          TLOADS(NTOTVS,2),             WORK1S(*),
     .          DISPLS(NTOTVMS),              PWORKS(NPOINS,3)
      DIMENSION PREASS(NPREAS,NPOINS),        TGAPSS(NPOINS),
     .          COORDS(NDIMES,NPOINS)
      DIMENSION PRESCS(NTOTVS,2),             RLOAHS(NTOTVS,NFUNCS),
     .          FICTOS(NFUNCS),               TFICTS(NFUNCS)
      DIMENSION TEMPIS(NPOINS,2),             ADVELS(NTOTVS*NDIMES),
     .          FPCHAS(NFPCH,NPOINS),         LPNTNS(*)
C
      IF(KARCLS.EQ.0.OR.ISTEPS.LE.2) THEN
C
C**** UPDATE DTIME AND TTIME
C
       CALL TMARCHS
C
C**** INCREMENT LOAD FACTORS, APPLIED LOADS AND NON-TENSIONAL STRAINS
C
       CALL INCREMS(ELDATS,ELPRES,ELVARS,ELMATS,HEADSS,HTLODS,
     .              IFFIXS(1,1),  PRESCS,
     .              LNODSS,MATNOS,PROELS,PROPSS,RLOADS,RLOAHS,
     .              FICTOS,TFICTS,TLOADS)
C
C**** PREDICTOR PHASE
C
       CALL PREDICS(DISITS,DISPRS,DISTOS)
C
C**** CALCULATE RESIDUAL FORCES
C
       CALL RESIDUS(DISITS,DISPRS,DISTOS,ELDATS,ELPRES,ELVARS,ELMATS,
     .              HEADSS,IFFIXS,LNODSS,MATNOS,PROELS,PROPSS,REFORS,
     .            TLOADS,DISPLS,PWORKS,PREASS(1,1),TGAPSS,COORDS,TEMPIS,
     .              ADVELS,FPCHAS,WORK1S)
      ENDIF
C
      RETURN
      END 