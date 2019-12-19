      SUBROUTINE ITERATT(DISITT,DISPRT,DISTOT,ELDATT,ELPRET,ELVART,
     .                   ELMATT,HEADST,IFFIXT,PRESCT,LNODST,MATNOT,
     .                   PROELT,PROPST,REFORT,RLOADT,FICTOT,TLOADT,
     .                   DISPLT,PWORKT,PREAST,TGAPST,COORDT,TEMPIT,
     .                   ADVELT,FPCHAT,LACTIT,HTLODT,WORK1T)
C***********************************************************************
C
C**** THIS ROUTINE PERFORMS THE ITERATIVE CORRECTION ON THE RESPONSE
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
      INCLUDE 'nuec_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'
C
      DIMENSION DISITT(NTOTVT,2),      DISPRT(NTOTVT,3),
     .          DISTOT(NTOTVT,3),      ELDATT(NDATAT),
     .          ELPRET(NPREVT),        ELVART(NSTATT),
     .          ELMATT(NMATXT),        HEADST(NPOINT,4),
     .          IFFIXT(NTOTVT,2),      LNODST(NNODET,NELEMT),
     .          MATNOT(NELEMT),        PROELT(NPRELT,NGRUPT),
     .          PROPST(NPROPT,NMATST), REFORT(NTOTVT,2),
     .          RLOADT(NTOTVT),        TLOADT(NTOTVT,2),
     .          WORK1T(*),             DISPLT(NTOTVM),
     .          PWORKT(NPOINT,3),      COORDT(NDIMET,NPOINT)
      DIMENSION PREAST(NPREAT,NPOINT), TGAPST(NPOINT)
      DIMENSION PRESCT(NTOTVT,2),      FICTOT(NFUNCT)
      DIMENSION TEMPIT(NPOINT,2),      ADVELT(NTOTVT*NDIMET),
     .          FPCHAT(NFPCH,NPOINT),  LACTIT(NELEMT),
     .          HTLODT(NHLODT,NSUBFT,NFUNCT)
C
C**** UPDATE ITERATION COUNTER
C
      IITERT=IITERT+1
C
C**** ASSEMBLE GLOBAL COUPLED STIFFNESS MATRIX, IF NECESSARY
C
      CALL ALGORST     ! select solution flags KSTIF, KRESL
C
      IF(NMEMO7.EQ.0) THEN
       IF(KSTIFT.EQ.1)                ! evaluate a new stiffness matrix
     .  CALL STIFMXT(ELDATT,ELPRET,ELVART,ELMATT,LNODST,MATNOT,PROELT,
     .               PROPST,WORK1T,DISTOT(1,2),DISITT(1,1),COORDT,
     .               ADVELT,TEMPIT(1,1),PREAST(1,1),TGAPST,DISTOT(1,1),
     .               TEMPIT(1,2),DISPLT,FPCHAT,LACTIT)
      ENDIF
C
      IF(NMEMO6.EQ.0) THEN
       IF(KRESLT.EQ.1)        ! assemble global coupled stiffness matrix
     .  CALL ASELMTT(MATNOT,PROELT,PROPST,
     .               ELMATT(IMATXT(1)),ELMATT(IMATXT(2)),
     .               ELMATT(IMATXT(3)),ELMATT(IMATXT(4)),
     .               ELMATT(IMATXT(5)),ELMATT(IMATXT(6)))
      ENDIF
C
C**** INITIALISE THE ITERA.TEMPER. TAKING INTO ACCOUNT FIXITY CONDITIONS
C
      CALL INIDIST(DISITT,IFFIXT,PRESCT,FICTOT,RLOADT)
C
C**** SOLVE EQUATIONS FOR ITERATIVE TEMPERATURES
C
      CALL ESEPAST(DISITT,DISPRT,IFFIXT,LNODST,REFORT,RLOADT,TLOADT,
     .             WORK1T,
     .             MATNOT,PROELT,PROPST,ELDATT,ELPRET,
     .             ELVART,ELMATT,
     .             DISTOT,COORDT,
     .             ADVELT,TEMPIT,PREAST(1,1),TGAPST,DISPLT,FPCHAT,
     .             LACTIT)
C
C**** CALCULATES THERMAL RESIDUAL VECTOR
C
      CALL RESIDUT(DISITT,DISPRT,DISTOT,ELDATT,ELPRET,ELVART,ELMATT,
     .             HEADST,IFFIXT,LNODST,MATNOT,PROELT,PROPST,REFORT,
     .            TLOADT,DISPLT,PWORKT,PREAST(1,1),TGAPST,COORDT,TEMPIT,
     .             ADVELT,FPCHAT,LACTIT,HTLODT,WORK1T)
C
C**** CHECK FOR CONVERGENCE
C
      CALL CONVERT(DISITT,DISPRT,DISTOT,REFORT,TLOADT,IFFIXT)
C
C**** OUTPUT RESULTS
C
      CALL OUTPUTT(DISTOT,ELDATT,ELPRET,ELVART,ELMATT,HEADST,
     .             IFFIXT,LNODST,MATNOT,PROELT,PROPST,TLOADT,
     .             COORDT,FPCHAT,PREAST(NPRE1T+1,1),DISPLT,LACTIT,
     .             WORK1T)
C
      RETURN
      END
