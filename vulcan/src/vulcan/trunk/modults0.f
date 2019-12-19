      SUBROUTINE MODULTS0(LNODST,MATNOT,PROELT,PROPST,COORDT,HTLODT,
     .                    IFFIXT,PRESCT,RLOADT,RLOAHT,FICTOT,TFICTT,
     .                    DISITT,DISPRT,DISTOT,HEADST,REFORT,TLOADT,
     .                    LPNTNT,ELDATT,ELPRET,ELVART,ELMATT,DISPLT,
     .                    PWORKT,PREAST,TGAPST,TEMPIT,ADVELT,FPCHAT,
     .                    LACTIT,
     .                    WORK1T,
     .                    LNODSS,MATNOS,PROELS,PROPSS,COORDS,HTLODS,
     .                    IFFIXS,PRESCS,RLOADS,RLOAHS,FICTOS,TFICTS,
     .                    DISITS,DISPRS,DISTOS,HEADSS,REFORS,TLOADS,
     .                    LPNTNS,ELDATS,ELPRES,ELVARS,ELMATS,DISPLS,
     .                    PWORKS,PREASS,TGAPSS,TEMPIS,ADVELS,FPCHAS,
     .                    LACTIS,
     .                    WORK1S)
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8, ALLOCATABLE, INTENT(INOUT) :: WORK1T(:),WORK1S(:)

      INTERFACE
       INCLUDE 'intervt.inc'
       INCLUDE 'intervt0.inc'
       INCLUDE 'tmstept.inc'
       INCLUDE 'intervs.inc'
       INCLUDE 'couple0s.inc'
      END INTERFACE
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'    ! thermal-mechanical
      INCLUDE 'nued_om.f'    ! thermal-microstructural
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
      INCLUDE 'inpo_oms.f'
C
C**** THERMAL VARIABLES
C
      DIMENSION LNODST(NNODET,NELEMT), MATNOT(NELEMT),
     .          PROELT(NPRELT,NGRUPT),    
     .          PROPST(NPROPT,NMATST), COORDT(NDIMET,NPOINT),
     .          HTLODT(NHLODT,NFUNCT)
      DIMENSION IFFIXT(NTOTVT,2),      PRESCT(NTOTVT,2),
     .          RLOADT(NTOTVT),        RLOAHT(NTOTVT,NFUNCT),
     .          FICTOT(NFUNCT),        TFICTT(NFUNCT)
      DIMENSION DISITT(NTOTVT,2),  
     .          DISPRT(NTOTVT,3),      DISTOT(NTOTVT,3),
     .          HEADST(NPOINT,4),
     .          REFORT(NTOTVT,2),      TLOADT(NTOTVT,2),
     .          LPNTNT(NPOINT)
      DIMENSION ELDATT(NDATAT),        ELPRET(NPREVT),
     .          ELVART(NSTATT),
     .          ELMATT(NMATXT),        DISPLT(NTOTVM),
     .          PWORKT(NPOINT,3)
      DIMENSION PREAST(NPREAT,NPOINT), TGAPST(NPOINT),
     .          TEMPIT(NPOINT,2),      ADVELT(NTOTVT*NDIMET),
     .          FPCHAT(NFPCH,NPOINT),  LACTIT(NELEMT)
C
C**** MICROSTRUCTURAL VARIABLES
C
      DIMENSION LNODSS(NNODES,NELEMS), MATNOS(NELEMS),
     .          PROELS(NPRELS,NGRUPS),
     .          PROPSS(NPROPS,NMATSS), COORDS(NDIMES,NPOINS),
     .          HTLODS(NHLODS,NFUNCS)
      DIMENSION IFFIXS(NTOTVS,2),      PRESCS(NTOTVS,2),
     .          RLOADS(NTOTVS),        RLOAHS(NTOTVS,NFUNCS),
     .          FICTOS(NFUNCS),        TFICTS(NFUNCS)
      DIMENSION DISITS(NTOTVS,2),
     .          DISPRS(NTOTVS,3),      DISTOS(NTOTVS,3),
     .          HEADSS(NPOINS,4),
     .          REFORS(NTOTVS,2),      TLOADS(NTOTVS,2),
     .          LPNTNS(NPOINS)
      DIMENSION ELDATS(NDATAS),        ELPRES(NPREVS),
     .          ELVARS(NSTATS),
     .          ELMATS(NMATXS),        DISPLS(NTOTVMS),
     .          PWORKS(NPOINS,3)
      DIMENSION PREASS(NPOINS),        TGAPSS(NPOINS),
     .          TEMPIS(NPOINS),        ADVELS(NTOTVS*NDIMES),
     .          FPCHAS(NFPCH,NPOINS),  LACTIS(NELEMS)
C
C**** GLOBAL CONVERGENCE ON OUTPUT
C
      NCKGLO=0
C
      IF(IEVFI.EQ.0) THEN
C
C***********************************************************************
C
C**** START
C
C***********************************************************************
C
C**** STARTS THE COMPUTATIONS FOR THE THERMAL PROBLEM
C
       CALL STARTST(DISTOT,DISPRT,ELDATT,ELPRET,ELVART,ELMATT,
     .              HEADST,HTLODT,IFFIXT,LNODST,LPNTNT,MATNOT,
     .              PROELT,PROPST,REFORT,RLOADT,TLOADT,COORDT,
     .              DISPLT,PWORKT,PREAST,TGAPST,TEMPIT,ADVELT,
     .              FPCHAT,LACTIT,WORK1T,     1)
C
C**** ASSIGNS MACROSCOPICAL PROPERTIES TO MICRO PROPERTIES
C
       CALL ASSPROT(PROPST,PROPSS,     0)
C
C**** STARTS THE COMPUTATIONS FOR THE MICROSTRUCTURAL PROBLEM
C
       CALL STARTSS(DISTOS,DISPRS,ELDATS,ELPRES,ELVARS,ELMATS,
     .              HEADSS,HTLODS,IFFIXS,LNODSS,LPNTNS,MATNOS,
     .              PROELS,PROPSS,REFORS,RLOADS,TLOADS,COORDS,
     .              DISPLS,PWORKS,PREASS,TGAPSS,TEMPIS,ADVELS,
     .              FPCHAS,LACTIS,WORK1S,     1)
C
C**** ASSIGNS MICROSCOPICAL PROPERTIES TO THERMAL PROPERTIES
C
       CALL ASSPROT(PROPST,PROPSS,     1)
C
C**** CONTINUES THE COMPUTATIONS FOR THE THERMAL PROBLEM
C
       CALL STARTST(DISTOT,DISPRT,ELDATT,ELPRET,ELVART,ELMATT,
     .              HEADST,HTLODT,IFFIXT,LNODST,LPNTNT,MATNOT,
     .              PROELT,PROPST,REFORT,RLOADT,TLOADT,COORDT,
     .              DISPLT,PWORKT,PREAST,TGAPST,TEMPIT,ADVELT,
     .              FPCHAT,LACTIT,WORK1T,     2)
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
        CALL INTERVT0(DISTOT,ELDATT,ELPRET,ELVART,ELMATT,
     .                HTLODT,IFFIXT,PRESCT,LNODST,LPNTNT,
     .                MATNOT,PROELT,PROPST,RLOADT,RLOAHT,
     .                FICTOT,TFICTT,ADVELT,
     .                HEADST,TLOADT,COORDT,TEMPIT,FPCHAT,PREAST,
     .                DISPLT,LACTIT,WORK1T)
C
C**** INDEX TO PRINT INITIAL CONDITIONS
C
        IPRCOT=IPRCOT+1
C
C**** LOOP OVER TIME STEPS (FOR THERMAL PROBLEMS)
C
        DO ISTEPT = KSTEPT,NSTEPT
C
C**** ADVANCE ON TIME AND PREDICT RESPONSE FOR THE THERMAL PROBLEM
C
         DO ISSTEPT=1,NSSTEPT             !  SUBINCREMENTS
C
          CALL TMSTEPT(DISITT,DISPRT,DISTOT,ELDATT,ELPRET,
     .                 ELVART,ELMATT,HEADST,HTLODT,IFFIXT,
     .                 PRESCT,LNODST,MATNOT,PROELT,PROPST,
     .                 REFORT,RLOADT,RLOAHT,FICTOT,TFICTT,
     .                 TLOADT,DISPLT,PWORKT,PREAST,TGAPST,
     .                 COORDT,TEMPIT,ADVELT,FPCHAT,LACTIT,
     .                 LPNTNT,WORK1T)
C
C**** LOOP ON ITERATIONS FOR THE THERMAL PROBLEM
C
          DO WHILE
     .     ((NCHEKT.NE.0).AND.(IITERT.LE.MITERT))
C
C**** ITERATIVE CORRECTION FOR THE THERMAL PROBLEM
C
           CALL ITERATT(DISITT,DISPRT,DISTOT,ELDATT,ELPRET,
     .                  ELVART,ELMATT,HEADST,IFFIXT,PRESCT,
     .                  LNODST,MATNOT,PROELT,PROPST,REFORT,
     .                  RLOADT,FICTOT,TLOADT,DISPLT,PWORKT,
     .                  PREAST,TGAPST,COORDT,TEMPIT,ADVELT,
     .                  FPCHAT,LACTIT,HTLODT,WORK1T)
C
          ENDDO
C
          IF(NCHEKT.NE.0)
     .     CALL RUNENDT(' NOT CONVERGED                    ') 
C
C**** CHECK FOR ANALYSIS TERMINATION FOR THE THERMAL PROBLEM
C
          IF(NSSTEPT.GT.1)
     .     CALL FINISHT(DISPRT,DISTOT,ELPRET,ELVART,HEADST,HTLODT,
     .                  IFFIXT,REFORT,RLOADT,TLOADT,PWORKT,TEMPIT)
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
         IF(NSSTEPT.EQ.1)
     .    CALL FINISHT(DISPRT,DISTOT,ELPRET,ELVART,HEADST,HTLODT,
     .                 IFFIXT,REFORT,RLOADT,TLOADT,PWORKT,TEMPIT)
C
        ENDDO                                 ! END OF TIME STEPS LOOP
C
       ENDDO                                  ! END OF TIME INTERVAL 
C
      ELSE                                    ! ievfi=1
C
C***********************************************************************
C
C**** START
C
C***********************************************************************
C
C**** START THE COMPUTATIONS FOR THE THERMAL PROBLEM
C
       CALL STARTST(DISTOT,DISPRT,ELDATT,ELPRET,ELVART,ELMATT,
     .              HEADST,HTLODT,IFFIXT,LNODST,LPNTNT,MATNOT,
     .              PROELT,PROPST,REFORT,RLOADT,TLOADT,COORDT,
     .              DISPLT,PWORKT,PREAST,TGAPST,TEMPIT,ADVELT,
     .              FPCHAT,LACTIT,WORK1T,     1)
C
C**** ASSIGNS MACROSCOPICAL PROPERTIES TO MICRO PROPERTIES
C
       CALL ASSPROT(PROPST,PROPSS,     0)
C
C**** START THE COMPUTATIONS FOR THE MICROSTRUCTURAL PROBLEM
C
       CALL STARTSS(DISTOS,DISPRS,ELDATS,ELPRES,ELVARS,ELMATS,
     .              HEADSS,HTLODS,IFFIXS,LNODSS,LPNTNS,MATNOS,
     .              PROELS,PROPSS,REFORS,RLOADS,TLOADS,COORDS,
     .              DISPLS,PWORKS,PREASS,TGAPSS,TEMPIS,ADVELS,
     .              FPCHAS,LACTIS,WORK1S,     1)
C
C**** ASSIGNS MICROSCOPICAL PROPERTIES TO THERMAL PROPERTIES
C
       CALL ASSPROT(PROPST,PROPSS,     1)
C
C**** CHECKS MECHANICAL & THERMAL DATA
C
c      CALL CHECKPS(MATNOS,PROELS,PROPSS,
c    .              MATNOT,PROELT,PROPST)
C
C**** CONTINUES THE COMPUTATIONS FOR THE THERMAL PROBLEM
C
       CALL STARTST(DISTOT,DISPRT,ELDATT,ELPRET,ELVART,ELMATT,
     .              HEADST,HTLODT,IFFIXT,LNODST,LPNTNT,MATNOT,
     .              PROELT,PROPST,REFORT,RLOADT,TLOADT,COORDT,
     .              DISPLT,PWORKT,PREAST,TGAPST,TEMPIT,ADVELT,
     .              FPCHAT,LACTIT,WORK1T,     2)
C
C**** CONTINUES THE COMPUTATIONS FOR THE MICROSTRUCTURAL PROBLEM
C
       CALL STARTSS(DISTOS,DISPRS,ELDATS,ELPRES,ELVARS,ELMATS,
     .              HEADSS,HTLODS,IFFIXS,LNODSS,LPNTNS,MATNOS,
     .              PROELS,PROPSS,REFORS,RLOADS,TLOADS,COORDS,
     .              DISPLS,PWORKS,PREASS,TGAPSS,TEMPIS,ADVELS,
     .              FPCHAS,LACTIS,WORK1S,     2)
C
C***********************************************************************
C
C**** ESTABLISH THE INITIAL CONDITIONS FOR THE MECHANICAL PROBLEM
C
C***********************************************************************
C
c      IF(INICOU.EQ.0)
c    .  CALL COUINI(LNODS,MATNO,PROEL,PROPS,COORD,HTLOD,
c    .              IFFIX,PRESC,RLOAD,RLOAH,FICTO,TFICT,
c    .              DISIT,DISPR,DISTO,HEADS,REFOR,TLOAD,
c    .              LPNTN,ELDAT,ELPRE,ELVAR,ELMAT,TEMPN,
c    .              DTEMP,INFRI,COFRI,PWORK,PREAS,TGAPS,
c    .              VNORM,FPCHA,LACTI,DACTI,
c    .              WORK1,
c    .              LNODST,MATNOT,PROELT,PROPST,COORDT,HTLODT,
c    .              IFFIXT,PRESCT,RLOADT,RLOAHT,FICTOT,TFICTT,
c    .              DISITT,DISPRT,DISTOT,HEADST,REFORT,TLOADT,
c    .              LPNTNT,ELDATT,ELPRET,ELVART,ELMATT,DISPLT,
c    .              PWORKT,PREAST,TGAPST,TEMPIT,ADVELT,FPCHAT,
c    .              LACTIT,
c    .              WORK1T)
C
C***********************************************************************
C
C**** SOLUTION OF THE COUPLED PROBLEM
C
C***********************************************************************
C
C**** INITIALISE SOME MECHANICAL PARAMETERS
C
c      IF(IREST.EQ.0.AND.INITI.NE.2) THEN
c       ITIME=0
c       TTIME=0.0D+00
c      ENDIF
C
C**** LOOP ON TIME INTERVALS (FOR BOTH PROBLEMS)
C
       DO WHILE (.TRUE.)                                
C
C**** START TIME INTERVAL FOR THE THERMAL PROBLEM
C
        CALL INTERVT(DISTOT,ELDATT,ELPRET,ELVART,ELMATT,
     .               HTLODT,IFFIXT,PRESCT,LNODST,LPNTNT,
     .               MATNOT,PROELT,PROPST,RLOADT,RLOAHT,
     .               FICTOT,TFICTT,ADVELT,
     .               HEADST,TLOADT,COORDT,TEMPIT,FPCHAT,PREAST,
     .               DISPLT,LACTIT,WORK1T)
C
C**** INDEX TO PRINT INITIAL CONDITIONS
C
        IPRCOT=IPRCOT+1
C
C**** START TIME INTERVAL FOR THE MICROSTRUCTURAL PROBLEM
C
        CALL INTERVS(DISTOS,ELDATS,ELPRES,ELVARS,ELMATS,
     .               HTLODS,IFFIXS,PRESCS,LNODSS,LPNTNS,
     .               MATNOS,PROELS,PROPSS,RLOADS,RLOAHS,
     .               FICTOS,TFICTS,ADVELS,
     .               HEADSS,TLOADS,COORDS,TEMPIS,FPCHAS,PREASS,
     .               DISPLS,LACTIS,WORK1S)
C
C**** INDEX TO PRINT INITIAL CONDITIONS
C
        IPRCOS=IPRCOS+1
C
C**** LOOP OVER TIME STEPS (FOR BOTH PROBLEMS)
C
        DO ISTEPS = KSTEPS,NSTEPS
C
C**** ASSIGN EQUAL TIME STEPS FOR BOTH PROBLEMS
C
         ISTEPT=ISTEPS
C
C**** SOME CONTROLS
C
         IF(NSTEPS.NE.NSTEPT)
     .    CALL RUNENDS('ERROR: CONTROL NUMBER OF TIME STEPS')
         IF(NSSTEPS*DTIMES.NE.NSSTEPT*DTIMET)
     .    CALL RUNENDS('ERROR: CHANGE SUBINCREMENTATION & TIME INCREM.')
C
         CALL COUPLE0S(LNODSS,MATNOS,PROELS,PROPSS,COORDS,HTLODS,
     .                 IFFIXS,PRESCS,RLOADS,RLOAHS,FICTOS,TFICTS,
     .                 DISITS,DISPRS,DISTOS,HEADSS,REFORS,TLOADS,
     .                 LPNTNS,ELDATS,ELPRES,ELVARS,ELMATS,DISPLS,
     .                 PWORKS,PREASS,TGAPSS,TEMPIS,ADVELS,FPCHAS,
     .                 LACTIS,
     .                 WORK1S,
     .                 LNODST,MATNOT,PROELT,PROPST,COORDT,HTLODT,
     .                 IFFIXT,PRESCT,RLOADT,RLOAHT,FICTOT,TFICTT,
     .                 DISITT,DISPRT,DISTOT,HEADST,REFORT,TLOADT,
     .                 LPNTNT,ELDATT,ELPRET,ELVART,ELMATT,DISPLT,
     .                 PWORKT,PREAST,TGAPST,TEMPIT,ADVELT,FPCHAT,
     .                 LACTIT,
     .                 WORK1T)
C
C***********************************************************************
C
C**** STAGGERED SOLUTION
C
C***********************************************************************
C
C     ISTAGGS=0    >>    UNIDIRECTIONAL COUPLING IN THE CURRENT STEP
C                        (INCREMENTAL STAGGERED SOLUTION)
C
C     ISTAGGS=1    >>    BIDIRECTIONAL COUPLING IN THE CURRENT STEP
C                        (ITERATIVE STAGGERED SOLUTION)
C
C***********************************************************************
C
         IF(ISTAGGS.EQ.1) THEN
C
C**** SOME CONTROL
C
          IF(NSSTEPS.GT.1)
     .     CALL RUNENDS('ERROR: ISTAGGS=1 & NSSTEPS >1 IS INCOMPATIBLE')
          IF(NSSTEPT.GT.1)
     .     CALL RUNENDT('ERROR: ISTAGGS=1 & NSSTEPT >1 IS INCOMPATIBLE')
C
          CALL COUSTA0S(LNODSS,MATNOS,PROELS,PROPSS,COORDS,HTLODS,
     .                  IFFIXS,PRESCS,RLOADS,RLOAHS,FICTOS,TFICTS,
     .                  DISITS,DISPRS,DISTOS,HEADSS,REFORS,TLOADS,
     .                  LPNTNS,ELDATS,ELPRES,ELVARS,ELMATS,DISPLS,
     .                  PWORKS,PREASS,TGAPSS,TEMPIS,ADVELS,FPCHAS,
     .                  LACTIS,
     .                  WORK1S,
     .                  LNODST,MATNOT,PROELT,PROPST,COORDT,HTLODT,
     .                  IFFIXT,PRESCT,RLOADT,RLOAHT,FICTOT,TFICTT,
     .                  DISITT,DISPRT,DISTOT,HEADST,REFORT,TLOADT,
     .                  LPNTNT,ELDATT,ELPRET,ELVART,ELMATT,DISPLT,
     .                  PWORKT,PREAST,TGAPST,TEMPIT,ADVELT,FPCHAT,
     .                  LACTIT,
     .                  WORK1T)
C
         ENDIF
C
C***********************************************************************
C
C**** FINISH
C
C***********************************************************************
C
C**** CHECK FOR ANALYSIS TERMINATION FOR THE THERMAL PROBLEM
C
         IF(NSSTEPT.EQ.1)
     .    CALL FINISHT(DISPRT,DISTOT,ELPRET,ELVART,HEADST,HTLODT,
     .                 IFFIXT,REFORT,RLOADT,TLOADT,PWORKT,TEMPIT)
C
C**** CHECK FOR ANALYSIS TERMINATION FOR THE MECHANICAL PROBLEM
C
         IF(NSSTEPS.EQ.1)
     .    CALL FINISHS(DISPRS,DISTOS,ELPRES,ELVARS,HEADSS,HTLODS,
     .                 IFFIXS,REFORS,RLOADS,TLOADS,PWORKS,TEMPIS)
C
        ENDDO                                 ! END OF TIME STEPS LOOP
C
       ENDDO                                  ! END OF TIME INTERVAL 
C
      ENDIF                                   ! ievfi.eq.0
C
      END
