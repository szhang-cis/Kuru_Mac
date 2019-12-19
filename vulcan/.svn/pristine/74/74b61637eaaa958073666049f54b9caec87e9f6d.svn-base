      SUBROUTINE RSTAR4T(TEMPIT,DISTOT,ELPRET,ELVART,TLOADT)
C***********************************************************************
C
C**** THIS ROUTINE REFORMS THE FOLLOWING TASKS:
C
C     NEW RUN OR RESTART: WRITES THE RESULTS OF THE CURRENT TIME STEP
C                         INTO *.fan file (future analysis file)
C
C     OPTIONS IN THE INPUT DATA FILE (see check0t.f):
C     1) START,NOINITIAL,FUTURE_ANALYSIS=FORM1                 (new run)
C     2) START,INITIAL,FUTURE_ANALYSIS=FORM1                   (new run)
C     3) START,PREVIOUS_ANALYSIS,time of the previous analysis, \
C        FUTURE_ANALYSIS=FORM1                                 (restart)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'
C
      DIMENSION TEMPIT(NTOTVT,2),             DISTOT(NTOTVT,3)
      DIMENSION ELPRET(NPREVT),               ELVART(NSTATT)
      DIMENSION TLOADT(NTOTVT,2)
C
      CALL CPUTIMT(TIME1T)
C
C**** OPEN FUTURE ANALYSIS FILE
C
      OPEN(UNIT=LUFANT,FILE=COT,STATUS='OLD',ACCESS='DIRECT',
     .     FORM='UNFORMATTED',RECL=LENRCT)
C
C**** READ FROM PROCESS AREA AND WRITE TO CONVERGED AREA ELPRE & ELVAR
C
      DO IELEMT=1,NELEMT
C
C**** READ ELPRE & ELVAR ( CURRENT )
C
       IF(NMEMO.EQ.1)
     .  CALL DATBAST(ELPRET,    2,    2)
       IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
     .  CALL DATBAST(ELVART,    3,    2)
C
C**** WRITE ELPRE & ELVAR
C
       IF(NMEMO.EQ.1)
     .  CALL DATRSTT(ELPRET,    2,    1)
       IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
     .  CALL DATRSTT(ELVART,    3,    1)
      ENDDO
C
C**** WRITE ALL MATRICES
C
      CALL DATRSTT(DISTOT,   4,     1)
C
      CALL DATRSTT(TEMPIT,   5,     1)               ! instead of disprt
C
      CALL DATRSTT(TLOADT,   8,     1)
C
C**** CLOSE FUTURE ANALYSIS FILE
C
      CLOSE(LUFANT)
C
      CALL CPUTIMT(TIME2T)
      CPURST=CPURST+(TIME2T-TIME1T)
      RETURN
C
      END
