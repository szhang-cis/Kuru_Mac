      SUBROUTINE MODUSS0(LNODSS,MATNOS,PROELS,PROPSS,COORDS,HTLODS,
     .                   IFFIXS,PRESCS,RLOADS,RLOAHS,FICTOS,TFICTS,
     .                   DISITS,DISPRS,DISTOS,HEADSS,REFORS,TLOADS,
     .                   LPNTNS,ELDATS,ELPRES,ELVARS,ELMATS,DISPLS,
     .                   PWORKS,PREASS,TGAPSS,TEMPIS,ADVELS,FPCHAS,
     .                   LACTIS,
     .                   WORK1S)
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8, ALLOCATABLE, INTENT(INOUT) :: WORK1S(:)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'   ! thermal-mechanical
      INCLUDE 'nued_om.f'   ! thermal-microstructural
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
      INCLUDE 'inpo_oms.f'
      INTERFACE
       INCLUDE 'intervs.inc'
      END INTERFACE
C
      DIMENSION LNODSS(NNODES,NELEMS),  MATNOS(NELEMS),
     .          PROELS(NPRELS,NGRUPS),  PROPSS(NPROPS,NMATSS),
     .          COORDS(NDIMES,NPOINS),  HTLODS(NHLODS,NSUBFS,NFUNCS)
      DIMENSION IFFIXS(NTOTVS,2),       PRESCS(NTOTVS,2),
     .          RLOADS(NTOTVS),         RLOAHS(NTOTVS,NFUNCS),
     .          FICTOS(NFUNCS),         TFICTS(NFUNCS)
      DIMENSION DISITS(NTOTVS,2),       DISPRS(NTOTVS,3),
     .          DISTOS(NTOTVS,3),       HEADSS(NPOINS,4),
     .          REFORS(NTOTVS,2),       TLOADS(NTOTVS,2)
      DIMENSION LPNTNS(NPOINS),         ELDATS(NDATAS),
     .          ELPRES(NPREVS),         ELVARS(NSTATS),
     .          ELMATS(NMATXS),         DISPLS(NTOTVMS),
     .          PWORKS(NPOINS,3)
      DIMENSION PREASS(NPREAS,NPOINS),  TGAPSS(NPOINS),
     .          TEMPIS(NPOINS,2),       ADVELS(NTOTVS*NDIMES),
     .          FPCHAS(NFPCH,NPOINS),   LACTIS(NELEMS)
C
C**** GLOBAL CONVERGENCE ON OUTPUT
C
      NCKGLO=0
C
C***********************************************************************
C
C**** START
C
C***********************************************************************
C
C**** START THE COMPUTATIONS FOR THE THERMAL PROBLEM
C
      CALL STARTSS(DISTOS,DISPRS,ELDATS,ELPRES,ELVARS,ELMATS,
     .             HEADSS,HTLODS,IFFIXS,LNODSS,LPNTNS,MATNOS,
     .             PROELS,PROPSS,REFORS,RLOADS,TLOADS,COORDS,
     .             DISPLS,PWORKS,PREASS,TGAPSS,TEMPIS,ADVELS,
     .             FPCHAS,LACTIS,WORK1S,     1)
C
C**** CONTINUES THE COMPUTATIONS FOR THE THERMAL PROBLEM
C
      CALL STARTSS(DISTOS,DISPRS,ELDATS,ELPRES,ELVARS,ELMATS,
     .             HEADSS,HTLODS,IFFIXS,LNODSS,LPNTNS,MATNOS,
     .             PROELS,PROPSS,REFORS,RLOADS,TLOADS,COORDS,
     .             DISPLS,PWORKS,PREASS,TGAPSS,TEMPIS,ADVELS,
     .             FPCHAS,LACTIS,WORK1S,     2)
C
C***********************************************************************
C
C**** SOLUTION OF THE THERMAL PROBLEM
C
C***********************************************************************
C
C**** LOOP ON TIME INTERVALS (FOR THERMAL PROBLEM)
C
      DO WHILE (.TRUE.)                                
C
C**** START TIME INTERVAL FOR THE THERMAL PROBLEM
C
       CALL INTERVS(DISTOS,ELDATS,ELPRES,ELVARS,ELMATS,
     .              HTLODS,IFFIXS,PRESCS,LNODSS,LPNTNS,
     .              MATNOS,PROELS,PROPSS,RLOADS,RLOAHS,
     .              FICTOS,TFICTS,ADVELS,
     .              HEADSS,TLOADS,COORDS,TEMPIS,FPCHAS,PREASS,
     .              DISPLS,LACTIS,WORK1S)
C
C**** INDEX TO PRINT INITIAL CONDITIONS
C
       IPRCOS=IPRCOS+1
C
C**** LOOP OVER TIME STEPS (FOR THERMAL PROBLEMS)
C
       DO ISTEPS = KSTEPS,NSTEPS
C
C**** ADVANCE ON TIME AND PREDICT RESPONSE FOR THE THERMAL PROBLEM
C
        DO ISSTEPS=1,NSSTEPS             !  SUBINCREMENTS
C
         CALL TMSTEPS(DISITS,DISPRS,DISTOS,ELDATS,ELPRES,
     .               ELVARS,ELMATS,HEADSS,HTLODS,IFFIXS,
     .               PRESCS,LNODSS,MATNOS,PROELS,PROPSS,
     .               REFORS,RLOADS,RLOAHS,FICTOS,TFICTS,
     .               TLOADS,DISPLS,PWORKS,PREASS,TGAPSS,
     .               COORDS,TEMPIS,ADVELS,FPCHAS,
     .               LPNTNS,WORK1S)
C
C**** LOOP ON ITERATIONS FOR THE THERMAL PROBLEM
C
         DO WHILE
     .    ((NCHEKS.NE.0).AND.(IITERS.LE.MITERS))
C
C**** ITERATIVE CORRECTION FOR THE THERMAL PROBLEM
C
          CALL ITERATS(DISITS,DISPRS,DISTOS,ELDATS,ELPRES,
     .                 ELVARS,ELMATS,HEADSS,IFFIXS,PRESCS,
     .                 LNODSS,MATNOS,PROELS,PROPSS,REFORS,
     .                 RLOADS,FICTOS,TLOADS,DISPLS,PWORKS,
     .                 PREASS,TGAPSS,COORDS,TEMPIS,ADVELS,
     .                 FPCHAS,WORK1S)
C
         ENDDO
C
         IF(NCHEKS.NE.0)
     .    CALL RUNENDS(' NOT CONVERGED                    ') 
C
C**** CHECK FOR ANALYSIS TERMINATION FOR THE THERMAL PROBLEM
C
         IF(NSSTEPS.GT.1)
     .    CALL FINISHS(DISPRS,DISTOS,ELPRES,ELVARS,HEADSS,HTLODS,
     .                 IFFIXS,REFORS,RLOADS,TLOADS,PWORKS,TEMPIS)
C
        ENDDO             ! LOOP OF SUBINCREMENTS
C
C***********************************************************************
C
C**** FINISH
C
C***********************************************************************
C
C**** CHECK FOR ANALYSIS TERMINATION FOR THE THERMAL PROBLEM
C
        IF(NSSTEPS.EQ.1)
     .   CALL FINISHS(DISPRS,DISTOS,ELPRES,ELVARS,HEADSS,HTLODS,
     .                IFFIXS,REFORS,RLOADS,TLOADS,PWORKS,TEMPIS)
C
       ENDDO                                 ! END OF TIME STEPS LOOP
C
      ENDDO                                  ! END OF TIME INTERVAL 
C
      END
