      SUBROUTINE STARTST(DISTOT,DISPRT,ELDATT,ELPRET,ELVART,ELMATT,
     .                   HEADST,HTLODT,IFFIXT,LNODST,LPNTNT,MATNOT,
     .                   PROELT,PROPST,REFORT,RLOADT,TLOADT,COORDT,
     .                   DISPLT,PWORKT,PREAST,TGAPST,TEMPIT,ADVELT,
     .                   FPCHAT,LACTIT,WORK1T,INDEXS)
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
      INCLUDE 'addi_omt.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
      INCLUDE 'nued_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION MATNOT(NELEMT),        LNODST(NNODET,NELEMT),
     .          PROELT(NPRELT,NGRUPT), PROPST(NPROPT,NMATST),
     .          ELDATT(NDATAT),        ELPRET(NPREVT),
     .          ELVART(NSTATT),        ELMATT(NMATXT),
     .          COORDT(NDIMET,NPOINT), WORK1T(*)
      DIMENSION DISTOT(NTOTVT,*),      DISPRT(NTOTVT,*), 
     .          HEADST(NPOINT,*),      HTLODT(NHLODT,NSUBFT,*), 
     .          IFFIXT(NTOTVT,*),      LPNTNT(*),
     .          REFORT(*),             RLOADT(*), 
     .          TLOADT(NTOTVT,*)
      DIMENSION DISPLT(NTOTVM),        PWORKT(NPOINT,3),
     .          PREAST(NPREAT,NPOINT), TGAPST(NPOINT),
     .          TEMPIT(NPOINT,2),      ADVELT(NTOTVT*NDIMET),
     .          FPCHAT(NFPCH,NPOINT),  LACTIT(NELEMT)
C
      IF(IRESTT.EQ.0) THEN
C***********************************************************************
C
C**** THIS IS A FRESH RUN
C
C***********************************************************************
       IF(INDEXS.EQ.1) THEN
C
C**** READ IN MOST OF THE PROBLEM DATA
C
        CALL INPDATT(ELDATT,ELPRET,ELVART,ELMATT,LNODST,MATNOT,PROELT,
     .               PROPST,COORDT,WORK1T)
C
C**** RENUMBER NODES
C
        CALL RENUMNT(KRENUT,LNODST,LPNTNT,NELEMT,NNODET,NPOINT,WORK1T)  
C
C**** WRITE TO DATA BASE
C
        IF(NFURES.EQ.2)
     .   CALL RSTAR2T(LPNTNT,LNODST,MATNOT,PROELT,PROPST,    1)
C
C**** INITIALIZE GLOBAL AND ELEMENTAL VARIABLES
C
        CALL ZEROTET(ELPRET,ELVART,DISTOT,HEADST,TLOADT,RLOADT,
     .               REFORT,DISPLT,PWORKT,PREAST,TGAPST,TEMPIT,
     .               ADVELT,FPCHAT,LACTIT)
C
C**** READ IN INITIAL GLOBAL AND ELEMENTAL VARIABLES
C
        IF(INITIT.EQ.1)
     .   CALL PREVOST(DISTOT,ELDATT,ELPRET,ELVART,ELMATT,HEADST,LNODST,
     .                MATNOT,PROELT,PROPST,TEMPIT(1,1),WORK1T)
        IF(INITIT.EQ.2)
     .   CALL RSTAR5T(TEMPIT,DISTOT,ELPRET,ELVART,TLOADT)
       ELSE
C
C**** CALCULATE THE CONSTANT MATRICES AND ARRAYS OF THE PROBLEM
C
        CALL SETMTXT(ELDATT,ELPRET,ELVART,ELMATT,LNODST,MATNOT,PROELT,
     .               PROPST,COORDT,TEMPIT(1,1),FPCHAT,DISPLT,
     .               WORK1T)
C
C**** CALCULATES INITIAL PHASE-CHANGE FUNCTION FOR MICRO TC FLOW PROB.
C
        IF(ITERMEF.GT.0.AND.IMICR.EQ.1) THEN
         CALL OUTPUTT(DISTOT,ELDATT,ELPRET,ELVART,ELMATT,HEADST,
     .                IFFIXT,LNODST,MATNOT,PROELT,PROPST,TLOADT,
     .                COORDT,FPCHAT,PREAST(NPRE1T+1,1),DISPLT,LACTIT,
     .                WORK1T)
        ENDIF
C
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
C**** READ GENERAL DATA FROM DATA BASE
C
        CALL RSTAR2T(LPNTNT,LNODST,MATNOT,PROELT,PROPST,    2)
C
C**** READ LAST CONVERGED VALUES FROM DATA BASE
C
        CALL RSTAR3T(DISPRT,DISTOT,ELPRET,ELVART,HEADST,HTLODT,IFFIXT,
     .               REFORT,RLOADT,TLOADT,    2)
       ENDIF                  ! indexs.eq.1
C
      ENDIF                  ! irestt.eq.0
C
C***********************************************************************
C
C**** DUMP GPCOD TO TAPE FOR POSTPROCESS IF NECESSARY
C
      IF(INDEXS.EQ.1) RETURN
C
      IF(NWPOST.EQ.1) THEN
       IF(NMEMO1.EQ.0) THEN
        CALL OUTPOST(ELDATT,LNODST,MATNOT,PROELT,PROPST,
     .               WORK1T(ISTART(1)),WORK1T(ISTART(2)),
     .               ELPRET,ELVART,ELMATT,DISPLT,WORK1T)
       ELSE
        CALL OUTPOST(ELDATT,LNODST,MATNOT,PROELT,PROPST,
     .               COORDT,WORK1T(ISTART(1)),
     .               ELPRET,ELVART,ELMATT,DISPLT,WORK1T)
       ENDIF
      ENDIF
C
      RETURN
      END
