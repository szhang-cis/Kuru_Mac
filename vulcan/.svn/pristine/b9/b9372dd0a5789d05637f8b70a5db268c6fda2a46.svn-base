      SUBROUTINE ESEPASS(DISITS,DISPRS,IFFIXS,LNODSS,REFORS,RLOADS,
     .                   TLOADS,WORK1S,
     .                   MATNOS,PROELS,PROPSS,ELDATS,ELPRES,
     .                   ELVARS,ELMATS,
     .                   DISTOS,COORDS,
     .                   ADVELS,TEMPIS,PREASS,TGAPSS,DISPLS,FPCHAS)
C***********************************************************************
C
C**** THIS ROUTINE PERFOMS THE SOLUTION FOR UNCOUPLED PROBLEM
C     (CONSTRAINED METHODS)
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
      DIMENSION DISITS(NTOTVS,*), DISPRS(NTOTVS,*), IFFIXS(NTOTVS,*),
     .          LNODSS(NNODES,*), REFORS(NTOTVS,*), RLOADS(*),
     .          TLOADS(NTOTVS,*)
      DIMENSION WORK1S(*)
C
      DIMENSION MATNOS(NELEMS),
     .          PROELS(NPRELS,NGRUPS), PROPSS(NPROPS,NMATSS),
     .          ELDATS(NDATAS),        ELPRES(NPREVS),
     .          ELVARS(NSTATS),        ELMATS(NMATXS)
      DIMENSION DISTOS(NTOTVS,3),      COORDS(NDIMES,NPOINS)
      DIMENSION ADVELS(NTOTVS*NDIMES), TEMPIS(NPOINS,2)
      DIMENSION PREASS(NPOINS),        TGAPSS(NPOINS)
      DIMENSION DISPLS(NTOTVMS),       FPCHAS(NFPCH,NPOINS)
C
C**** SOLVE EQUATIONS  ( FIRST PASS )
C
      CALL ANTESOS(REFORS(1,1),REFORS(1,2),DISITS(1,2),
     .             WORK1S(IQUASS(1)),WORK1S(IQUASS(2)))
C
      CALL SOLVERS(DISITS(1,1),REFORS(1,1),IFFIXS(1,1),
     .             IMPOSS,KPORES,KRESLS,
     .             LNODSS,NDOFNS,NELEMS,NEVABS,NNODES,NPOINS,
     .             NTOTVS,     1,IITERS,WORK1S,miters,
     .             NPRELS,NGRUPS,NPROPS,NMATSS,
     .             NDATAS,NPREVS,NSTATS,NMATXS,NDIMES,NTOTVMS,NFPCH,
     .             MATNOS,PROELS,PROPSS,ELDATS,ELPRES,
     .             ELVARS,ELMATS,
     .             DISTOS,COORDS,
     .             ADVELS,TEMPIS,PREASS,TGAPSS,DISPLS,FPCHAS)
C
      CALL POSTSOS(DISITS(1,1),DISITS(1,2),REFORS(1,1),REFORS(1,2),
     .             WORK1S(IQUASS(1)),WORK1S(IQUASS(2)))
C
      RETURN
      END
