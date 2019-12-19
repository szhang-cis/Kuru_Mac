      SUBROUTINE MODUL2(LNODS,MATNO,PROEL,PROPS,COORD,HTLOD,
     .                  IFFIX,PRESC,RLOAD,RLOAH,FICTO,TFICT,
     .                  DISIT,DISPR,DISTO,HEADS,REFOR,TLOAD,
     .                  LPNTN,ELDAT,ELPRE,ELVAR,ELMAT,TEMPN,
     .                  DTEMP,INFRI,COFRI,PWORK,PREAS,TGAPS,
     .                  VNORM,FPCHA,LACTI,NOPRF,PRESF,PREHF,
     .                  VANIS,
     .                  WORK1,
     .                  LNODST,MATNOT,PROELT,PROPST,COORDT,HTLODT,
     .                  IFFIXT,PRESCT,RLOADT,RLOAHT,FICTOT,TFICTT,
     .                  DISITT,DISPRT,DISTOT,HEADST,REFORT,TLOADT,
     .                  LPNTNT,ELDATT,ELPRET,ELVART,ELMATT,DISPLT,
     .                  PWORKT,PREAST,TGAPST,TEMPIT,ADVELT,FPCHAT,
     .                  LACTIT,
     .                  WORK1T)
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8, ALLOCATABLE, INTENT(INOUT) :: WORK1(:), WORK1T(:)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
      INCLUDE 'inpo_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'
      INTERFACE
        INCLUDE 'intervt.inc'
        INCLUDE 'couple2.inc'
        INCLUDE 'couini.inc'
        INCLUDE 'interv.inc'
        INCLUDE 'cousta2.inc'
      END INTERFACE
C
C**** MECHANICAL VARIABLES
C
      DIMENSION LNODS(NNODE,NELEM), MATNO(NELEM),
     .          PROEL(NPREL,NGRUP), PROPS(NPROP,NMATS),
     .          COORD(NDIME,NPOIN), HTLOD(NHLOD,NSUBF,NFUNC)
      DIMENSION IFFIX(NTOTV,2),     PRESC(NTOTV,2),
     .          RLOAD(NTOTV),       RLOAH(NTOTV,NFUNC),
     .          FICTO(NFUNC),       TFICT(NFUNC)
      DIMENSION DISIT(NTOTV,2),     DISPR(NTOTV,NDISR),
     .          DISTO(NTOTV,NDISO), HEADS(NPOIN,4),
     .          REFOR(NTOTV,2),     TLOAD(NTOTV,2)
      DIMENSION LPNTN(NPOIN),       ELDAT(NDATA),
     .          ELPRE(NPREV),       ELVAR(NSTAT),
     .          ELMAT(NMATX),       TEMPN(NPOIN,2),
     .          DTEMP(NPOIN)
      DIMENSION INFRI(NPOIN),       COFRI(NSKEW,NDIME,*),
     .          PWORK(NPOIN,2),     PREAS(NPREA,NPOIN),
     .          TGAPS(NPOIN),       VNORM(NTOTV),
     .          FPCHA(NFPCH,NPOIN), LACTI(NELEM)
      DIMENSION NOPRF(NNODE,NELEM), PRESF(NNODE,NDIME,NELEM),
     .          PREHF(NNODE,NDIME,NELEM,NFUNC),
     .          VANIS(NANIV,NANIC,NELEM)
C
C**** THERMAL VARIABLES
C
      DIMENSION LNODST(NNODET,NELEMT), MATNOT(NELEMT),
     .          PROELT(NPRELT,NGRUPT), PROPST(NPROPT,NMATST),
     .          COORDT(NDIMET,NPOINT), HTLODT(NHLODT,NSUBFT,NFUNCT)
      DIMENSION IFFIXT(NTOTVT,2),      PRESCT(NTOTVT,2),
     .          RLOADT(NTOTVT),        RLOAHT(NTOTVT,NFUNCT),
     .          FICTOT(NFUNCT),        TFICTT(NFUNCT)
      DIMENSION DISITT(NTOTVT,2),  
     .          DISPRT(NTOTVT,3),      DISTOT(NTOTVT,3),
     .          HEADST(NPOINT,4),
     .          REFORT(NTOTVT,2),      TLOADT(NTOTVT,2),
     .          LPNTNT(NPOINT)
      DIMENSION ELDATT(NDATAT),        ELPRET(NPREVT),
     .          ELVART(NSTATT),        ELMATT(NMATXT)
      DIMENSION DISPLT(NTOTVM),        PWORKT(NPOINT,3),
     .          PREAST(NPREAT,NPOINT), TGAPST(NPOINT),
     .          TEMPIT(NPOINT,2),      ADVELT(NTOTVT*NDIMET),
     .          FPCHAT(NFPCH,NPOINT),  LACTIT(NELEMT)
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
     .             HEADST,HTLODT,IFFIXT,LNODST,LPNTNT,MATNOT,
     .             PROELT,PROPST,REFORT,RLOADT,TLOADT,COORDT,
     .             DISPLT,PWORKT,PREAST,TGAPST,TEMPIT,ADVELT,
     .             FPCHAT,LACTIT,WORK1T,     1)
C
C**** START THE COMPUTATIONS FOR THE MECHANICAL PROBLEM
C
      CALL STARTS(DISTO,DISPR,ELDAT,ELPRE,ELVAR,ELMAT,
     .            HEADS,HTLOD,IFFIX,LNODS,LPNTN,MATNO,
     .            PROEL,PROPS,REFOR,RLOAD,TLOAD,COORD,
     .            TEMPN,DTEMP,PWORK,PREAS,TGAPS,VNORM,
     .            FPCHA,LACTI,PRESF,VANIS,WORK1,    1)
C
C**** ASSIGNS MECHANICAL PROPERTIES TO THERMAL PROPERTIES
C
      CALL ASSPROM(PROPS,PROPST)
C
C**** CHECKS MECHANICAL & THERMAL DATA
C
      CALL CHECKP(MATNO, PROEL, PROPS,
     .            MATNOT,PROELT,PROPST)
C
C**** CONTINUES THE COMPUTATIONS FOR THE THERMAL PROBLEM
C
      CALL STARTST(DISTOT,DISPRT,ELDATT,ELPRET,ELVART,ELMATT,
     .             HEADST,HTLODT,IFFIXT,LNODST,LPNTNT,MATNOT,
     .             PROELT,PROPST,REFORT,RLOADT,TLOADT,COORDT,
     .             DISPLT,PWORKT,PREAST,TGAPST,TEMPIT,ADVELT,
     .             FPCHAT,LACTIT,WORK1T,     2)
C
C**** CONTINUES THE COMPUTATIONS FOR THE MECHANICAL PROBLEM
C
      CALL STARTS(DISTO,DISPR,ELDAT,ELPRE,ELVAR,ELMAT,
     .            HEADS,HTLOD,IFFIX,LNODS,LPNTN,MATNO,
     .            PROEL,PROPS,REFOR,RLOAD,TLOAD,COORD,
     .            TEMPN,DTEMP,PWORK,PREAS,TGAPS,VNORM,
     .            FPCHA,LACTI,PRESF,VANIS,WORK1,    2)
C
C***********************************************************************
C
C**** ESTABLISH THE INITIAL CONDITIONS FOR THE MECHANICAL PROBLEM
C
C***********************************************************************
C
      IF(INICOU.EQ.0)
     . CALL COUINI(LNODS,MATNO,PROEL,PROPS,COORD,HTLOD,
     .             IFFIX,PRESC,RLOAD,RLOAH,FICTO,TFICT,
     .             DISIT,DISPR,DISTO,HEADS,REFOR,TLOAD,
     .             LPNTN,ELDAT,ELPRE,ELVAR,ELMAT,TEMPN,
     .             DTEMP,INFRI,COFRI,PWORK,PREAS,TGAPS,
     .             VNORM,FPCHA,LACTI,NOPRF,PRESF,PREHF,
     .             VANIS,WORK1,NEQNS,
     .             LNODST,MATNOT,PROELT,PROPST,COORDT,HTLODT,
     .             IFFIXT,PRESCT,RLOADT,RLOAHT,FICTOT,TFICTT,
     .             DISITT,DISPRT,DISTOT,HEADST,REFORT,TLOADT,
     .             LPNTNT,ELDATT,ELPRET,ELVART,ELMATT,DISPLT,
     .             PWORKT,PREAST,TGAPST,TEMPIT,ADVELT,FPCHAT,
     .             LACTIT,
     .             WORK1T)
C
C***********************************************************************
C
C**** SOLUTION OF THE COUPLED PROBLEM
C
C***********************************************************************
C
C**** INITIALISE SOME MECHANICAL PARAMETERS
C
c     IF(IREST.EQ.0.AND.INITI.NE.2) THEN
c      ITIME=0
c      TTIME=0.0D+00
c     ENDIF
C
C**** LOOP ON TIME INTERVALS (FOR BOTH PROBLEMS)
C
      DO WHILE (.TRUE.)                                
C
C**** START TIME INTERVAL FOR THE THERMAL PROBLEM
C
       CALL INTERVT(DISTOT,ELDATT,ELPRET,ELVART,ELMATT,
     .              HTLODT,IFFIXT,PRESCT,LNODST,LPNTNT,
     .              MATNOT,PROELT,PROPST,RLOADT,RLOAHT,
     .              FICTOT,TFICTT,ADVELT,
     .              HEADST,TLOADT,COORDT,TEMPIT,FPCHAT,PREAST,
     .              DISPLT,LACTIT,WORK1T)
C
C**** INDEX TO PRINT INITIAL CONDITIONS
C
       IPRCOT=IPRCOT+1
C
C**** START TIME INTERVAL FOR THE MECHANICAL PROBLEM
C
       CALL INTERV(DISTO,ELDAT,ELPRE,ELVAR,ELMAT,HTLOD,IFFIX,
     .             PRESC,LNODS,LPNTN,MATNO,PROEL,PROPS,RLOAD,
     .             RLOAH,FICTO,TFICT,INFRI,COFRI,COORD,HEADS,
     .             TLOAD,REFOR,PWORK,TGAPS,PREAS,LACTI,NOPRF,
     .             PREHF,VANIS,WORK1,NEQNS)
C
C**** INDEX TO PRINT INITIAL CONDITIONS
C
       IPRCO=IPRCO+1
C
C**** LOOP OVER TIME STEPS (FOR BOTH PROBLEMS)
C
       DO ISTEP = KSTEP,NSTEP
C
C**** ASSIGN EQUAL TIME STEPS FOR BOTH PROBLEMS
C
        ISTEPT=ISTEP
C
C**** SOME CONTROLS
C
        IF(NSTEP.NE.NSTEPT)
     .   CALL RUNEND('ERROR: CONTROL NUMBER OF TIME STEPS')
        IF((DABS(NSSTEP*DTIME-NSSTEPT*DTIMET)).GT.1.0D-10)
     .   CALL RUNEND('ERROR: CHANGE SUBINCREMENTATION AND TIME INCREM.')
C
        CALL COUPLE2(LNODS,MATNO,PROEL,PROPS,COORD,HTLOD,
     .               IFFIX,PRESC,RLOAD,RLOAH,FICTO,TFICT,
     .               DISIT,DISPR,DISTO,HEADS,REFOR,TLOAD,
     .               LPNTN,ELDAT,ELPRE,ELVAR,ELMAT,TEMPN,
     .               DTEMP,INFRI,COFRI,PWORK,PREAS,TGAPS,
     .               VNORM,FPCHA,LACTI,NOPRF,PRESF,PREHF,
     .               VANIS,WORK1,NEQNS,
     .               LNODST,MATNOT,PROELT,PROPST,COORDT,HTLODT,
     .               IFFIXT,PRESCT,RLOADT,RLOAHT,FICTOT,TFICTT,
     .               DISITT,DISPRT,DISTOT,HEADST,REFORT,TLOADT,
     .               LPNTNT,ELDATT,ELPRET,ELVART,ELMATT,DISPLT,
     .               PWORKT,PREAST,TGAPST,TEMPIT,ADVELT,FPCHAT,
     .               LACTIT,
     .               WORK1T)
C
C***********************************************************************
C
C**** STAGGERED SOLUTION
C
C***********************************************************************
C
C     ISTAGG=0    >>    UNIDIRECTIONAL COUPLING IN THE CURRENT STEP
C                       (INCREMENTAL STAGGERED SOLUTION)
C
C     ISTAGG=1    >>    BIDIRECTIONAL COUPLING IN THE CURRENT STEP
C                       (ITERATIVE STAGGERED SOLUTION)
C
C***********************************************************************
C
        IF(ISTAGG.EQ.1) THEN
C
C**** SOME CONTROL
C
         IF(NSSTEP.GT.1)
     .    CALL RUNEND('ERROR: ISTAGG=1 & NSSTEP > 1 IS INCOMPATIBLE')
         IF(NSSTEPT.GT.1)
     .    CALL RUNENDT('ERROR: ISTAGG=1 & NSSTEPT > 1 IS INCOMPATIBLE')
C
         CALL COUSTA2(LNODS,MATNO,PROEL,PROPS,COORD,HTLOD,
     .                IFFIX,PRESC,RLOAD,RLOAH,FICTO,TFICT,
     .                DISIT,DISPR,DISTO,HEADS,REFOR,TLOAD,
     .                LPNTN,ELDAT,ELPRE,ELVAR,ELMAT,TEMPN,
     .                DTEMP,INFRI,COFRI,PWORK,PREAS,TGAPS,
     .                VNORM,FPCHA,LACTI,NOPRF,PRESF,PREHF,
     .                VANIS,WORK1,NEQNS,
     .                LNODST,MATNOT,PROELT,PROPST,COORDT,HTLODT,
     .                IFFIXT,PRESCT,RLOADT,RLOAHT,FICTOT,TFICTT,
     .                DISITT,DISPRT,DISTOT,HEADST,REFORT,TLOADT,
     .                LPNTNT,ELDATT,ELPRET,ELVART,ELMATT,DISPLT,
     .                PWORKT,PREAST,TGAPST,TEMPIT,ADVELT,FPCHAT,
     .                LACTIT,
     .                WORK1T)
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
     .   CALL FINISHT(DISPRT,DISTOT,ELPRET,ELVART,HEADST,HTLODT,
     .                IFFIXT,REFORT,RLOADT,TLOADT,PWORKT,TEMPIT)
C
C**** CHECK FOR ANALYSIS TERMINATION FOR THE MECHANICAL PROBLEM
C
        IF(NSSTEP.EQ.1)
     .   CALL FINISH(DISPR,DISTO,ELPRE,ELVAR,HEADS,HTLOD,IFFIX,
     .               REFOR,RLOAD,TLOAD,INFRI,COFRI,LACTI)
C
       ENDDO                                 ! END OF TIME STEPS LOOP
C
      ENDDO                                  ! END OF TIME INTERVAL 
C
      END
