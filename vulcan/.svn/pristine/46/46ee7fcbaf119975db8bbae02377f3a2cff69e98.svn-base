      SUBROUTINE MODULT0(LNODST,MATNOT,PROELT,PROPST,COORDT,HTLODT,
     .                   IFFIXT,PRESCT,RLOADT,RLOAHT,FICTOT,TFICTT,
     .                   DISITT,DISPRT,DISTOT,HEADST,REFORT,TLOADT,
     .                   LPNTNT,ELDATT,ELPRET,ELVART,ELMATT,DISPLT,
     .                   PWORKT,PREAST,TGAPST,TEMPIT,ADVELT,FPCHAT,
     .                   LACTIT,
     .                   WORK1T)
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8, ALLOCATABLE, INTENT(INOUT)  :: WORK1T(:)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'   ! thermal-mechanical
      INCLUDE 'nued_om.f'   ! thermal-microstructural
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'

      INTERFACE
        INCLUDE 'intervt0.inc'
        INCLUDE 'tmstept.inc'
      END INTERFACE

C
C**** THERMAL VARIABLES
C
      DIMENSION LNODST(NNODET,NELEMT), MATNOT(NELEMT),
     .          PROELT(NPRELT,NGRUPT), PROPST(NPROPT,NMATST),
     .          COORDT(NDIMET,NPOINT), HTLODT(NHLODT,NSUBFT,NFUNCT)
      DIMENSION IFFIXT(NTOTVT,2),      PRESCT(NTOTVT,2),
     .          RLOADT(NTOTVT),        RLOAHT(NTOTVT,NFUNCT),
     .          FICTOT(NFUNCT),        TFICTT(NFUNCT)
      DIMENSION DISITT(NTOTVT,2),      DISPRT(NTOTVT,3),
     .          DISTOT(NTOTVT,3),      HEADST(NPOINT,4),
     .          REFORT(NTOTVT,2),      TLOADT(NTOTVT,2)
      DIMENSION LPNTNT(NPOINT),        ELDATT(NDATAT),
     .          ELPRET(NPREVT),        ELVART(NSTATT),
     .          ELMATT(NMATXT),        DISPLT(NTOTVM),
     .          PWORKT(NPOINT,3)
      DIMENSION PREAST(NPREAT,NPOINT), TGAPST(NPOINT),
     .          TEMPIT(NPOINT,2),      ADVELT(NTOTVT*NDIMET),
     .          FPCHAT(NFPCH,NPOINT),  LACTIT(NELEMT)
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
      CALL STARTST(DISTOT,DISPRT,ELDATT,ELPRET,ELVART,ELMATT,
     .             HEADST,HTLODT,IFFIXT,LNODST,LPNTNT,MATNOT,
     .             PROELT,PROPST,REFORT,RLOADT,TLOADT,COORDT,
     .             DISPLT,PWORKT,PREAST,TGAPST,TEMPIT,ADVELT,
     .             FPCHAT,LACTIT,WORK1T,     1)
C
C**** CONTINUES THE COMPUTATIONS FOR THE THERMAL PROBLEM
C
      CALL STARTST(DISTOT,DISPRT,ELDATT,ELPRET,ELVART,ELMATT,
     .             HEADST,HTLODT,IFFIXT,LNODST,LPNTNT,MATNOT,
     .             PROELT,PROPST,REFORT,RLOADT,TLOADT,COORDT,
     .             DISPLT,PWORKT,PREAST,TGAPST,TEMPIT,ADVELT,
     .             FPCHAT,LACTIT,WORK1T,     2)
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
C**** LOOP OVER TIME STEPS (FOR THERMAL PROBLEMS)
C
       DO ISTEPT = KSTEPT,NSTEPT
C
C**** ADVANCE ON TIME AND PREDICT RESPONSE FOR THE THERMAL PROBLEM
C
        DO ISSTEPT=1,NSSTEPT             !  SUBINCREMENTS
C
         CALL TMSTEPT(DISITT,DISPRT,DISTOT,ELDATT,ELPRET,
     .               ELVART,ELMATT,HEADST,HTLODT,IFFIXT,
     .               PRESCT,LNODST,MATNOT,PROELT,PROPST,
     .               REFORT,RLOADT,RLOAHT,FICTOT,TFICTT,
     .               TLOADT,DISPLT,PWORKT,PREAST,TGAPST,
     .               COORDT,TEMPIT,ADVELT,FPCHAT,LACTIT,
     .               LPNTNT,WORK1T)
C
C**** LOOP ON ITERATIONS FOR THE THERMAL PROBLEM
C
         DO WHILE
     .    ((NCHEKT.NE.0).AND.(IITERT.LE.MITERT))
C
C**** ITERATIVE CORRECTION FOR THE THERMAL PROBLEM
C
          CALL ITERATT(DISITT,DISPRT,DISTOT,ELDATT,ELPRET,
     .                 ELVART,ELMATT,HEADST,IFFIXT,PRESCT,
     .                 LNODST,MATNOT,PROELT,PROPST,REFORT,
     .                 RLOADT,FICTOT,TLOADT,DISPLT,PWORKT,
     .                 PREAST,TGAPST,COORDT,TEMPIT,ADVELT,
     .                 FPCHAT,LACTIT,HTLODT,WORK1T)
C
         ENDDO
C
         IF(NCHEKT.NE.0)
     .    CALL RUNENDT(' NOT CONVERGED                    ') 
C
C**** CHECK FOR ANALYSIS TERMINATION FOR THE THERMAL PROBLEM
C
         IF(NSSTEPT.GT.1)
     .    CALL FINISHT(DISPRT,DISTOT,ELPRET,ELVART,HEADST,HTLODT,
     .                 IFFIXT,REFORT,RLOADT,TLOADT,PWORKT,TEMPIT)
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
     .   CALL FINISHT(DISPRT,DISTOT,ELPRET,ELVART,HEADST,HTLODT,
     .                IFFIXT,REFORT,RLOADT,TLOADT,PWORKT,TEMPIT)
C
       ENDDO                                 ! END OF TIME STEPS LOOP
C
      ENDDO                                  ! END OF TIME INTERVAL 
C
      END
