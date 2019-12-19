      SUBROUTINE INTERVT0(DISTOT,ELDATT,ELPRET,ELVART,ELMATT,HTLODT,
     .                    IFFIXT,PRESCT,LNODST,LPNTNT,MATNOT,PROELT,
     .                    PROPST,RLOADT,RLOAHT,FICTOT,TFICTT,ADVELT,
     .                    HEADST,TLOADT,COORDT,TEMPIT,FPCHAT,PREAST,
     .                    DISPLT,LACTIT,WORK1T)
C***********************************************************************
C
C**** THIS ROUTINE SETS UP THE TYPE OF ALGORITHM AND REFERENCE LOAD
C     TO BE USED IN THE CURRENT TIME INTERVAL
C     FOR THE FIRST INTERVAL, IT ALSO EVALUATES THE CONSTANT ELEMENT
C     MATRICES OF THE PROBLEM
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8, ALLOCATABLE  :: WORK1T(:)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'   ! thermal-mechanical
      INCLUDE 'nuef_om.f'   ! thermal-flow
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'

      INCLUDE 'soladdt.inc'

C
      DIMENSION MATNOT(NELEMT),          LNODST(NNODET,NELEMT),
     .          PROELT(NPRELT,NGRUPT),   PROPST(NPROPT,NMATST),
     .          ELDATT(NDATAT),          ELPRET(NPREVT),
     .          ELVART(NSTATT),          ELMATT(NMATXT)
      DIMENSION DISTOT(*),
     .          HTLODT(NHLODT,NSUBFT,*), IFFIXT(NTOTVT,*),
     .          LPNTNT(*),               RLOADT(*)       
      DIMENSION PRESCT(NTOTVT,2),        RLOAHT(NTOTVT,NFUNCT),
     .          FICTOT(NFUNCT),          TFICTT(NFUNCT)
      DIMENSION ADVELT(NTOTVT*NDIMET)
      DIMENSION HEADST(NPOINT,4),        TLOADT(NTOTVT,2),
     .          COORDT(NDIMET,NPOINT),   TEMPIT(NPOINT,2),
     .          FPCHAT(NFPCH,NPOINT),    PREAST(NPREAT,NPOINT),
     .          DISPLT(NTOTVM),          LACTIT(NELEMT)
C
C**** LOOK FOR THE BEGINING OF THE INTERVAL DATA CARD  
C
      CALL INTREDT0(TFICTT)
C
C**** OPEN RESTART FILE
C
      IF(NFURES.EQ.2)
     . OPEN(UNIT=LURSTT,FILE=CKT,STATUS='UNKNOWN',ACCESS='DIRECT',
     .      FORM='UNFORMATTED',RECL=LENRCT)
C
C**** READ THE NEW LOAD, NEW BOUNDARY CONDITIONS AND SOLUTION STRATEGY
C     ( UNLESS THIS IS A "CONTINUE RESTART" )
C
      IF(IRESTT.NE.1.OR.(IRESTT.EQ.1.AND.ISKIPT.NE.0))
     . CALL INTLODT(ELDATT,ELPRET,ELVART,ELMATT,HTLODT,IFFIXT,PRESCT,
     .              LNODST,MATNOT,PROELT,PROPST,RLOADT,RLOAHT,ADVELT,
     .              TEMPIT,COORDT,FPCHAT,DISPLT,LACTIT,
     .              WORK1T)
C
C**** PRINTS INITIAL VALUES
C
      IF(IPRCOT.EQ.0)
     . CALL OUTPUTT(DISTOT,ELDATT,ELPRET,ELVART,ELMATT,HEADST,
     .              IFFIXT,LNODST,MATNOT,PROELT,PROPST,TLOADT,
     .              COORDT,FPCHAT,PREAST(NPRE1T+1,1),DISPLT,LACTIT,
     .              WORK1T)
C
C**** VIRTUAL MEMORY ALLOCATION FOR SOLVER
C        ( when new boundary or when re-starting or new active elements)
C
      IF(NEWBOT.EQ.1.OR.IRESTT.NE.0.OR.NEWACT.EQ.1)
     . CALL SOLADDT(IFFIXT,LNODST,LPNTNT,PRESCT,
     .              WORK1T(IFIXYT( 1)),WORK1T(IFIXYT( 2)),
     .              WORK1T(IFIXYT( 3)),WORK1T(IFIXYT( 4)),
     .              WORK1T(IFIXYT( 5)),WORK1T(IFIXYT( 6)),
     .              WORK1T(IFIXYT( 7)),WORK1T)
C
      IF(ITIMET.EQ.1.OR.IRESTT.NE.0) THEN
C                             ! 1st interval or when re-starting
C
C**** DEAL WITH VIRTUAL MEMORY ALLOCATION 
C
       CALL ADDDATT(DISTOT,ELDATT,ELPRET,ELVART)
       LTOTLT=LPRINT+LWOR1T+LSOLVT+LDABAT
       WRITE(LUREST,900) LPRINT,LWOR1T,LSOLVT,LDABAT,LTOTLT
       WRITE(LUPRIT,900) LPRINT,LWOR1T,LSOLVT,LDABAT,LTOTLT
C
C**** COMPUTE CONSTANT ELEMENT MATRICES
C
       CALL CONMTXT(ELDATT,ELPRET,ELVART,ELMATT,LNODST,MATNOT,PROELT,
     .              PROPST,WORK1T)
C
      ENDIF
C
C**** CLOSE RESTART FILE
C
      IF(NFURES.EQ.2) THEN
       IF(KSAVET.EQ.-1) THEN
        CLOSE(LURSTT,STATUS='DELETE')
       ELSE
        CLOSE(LURSTT,STATUS='KEEP')
       ENDIF
      ENDIF
C
      RETURN
  900 FORMAT(1H1,//5X,'SUMMARY OF MEMORY REQUIREMENTS :',/,
     .             5X,'==============================  ',/,
     .            15X,'REQUIRED PERMANENT MEMORY  =',I8,' WORDS (R*8)'/
     .            15X,'REQUIRED TEMPORARY MEMORY  =',I8,' WORDS (R*8)'/
     .            15X,'REQUIRED SOLVER    MEMORY  =',I8,' WORDS (R*8)'/
     .            15X,'ALLOCAT. DATA BASE MEMORY  =',I8,' WORDS (R*8)'/
     .            15X,'ALLOCAT. TOTAL     MEMORY  =',I8,' WORDS (R*8)'/)
      END
