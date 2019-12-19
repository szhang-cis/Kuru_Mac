      SUBROUTINE STOREST(DISPRT,DISTOT,ELPRET,ELVART,HEADST,HTLODT,
     .                   IFFIXT,REFORT,RLOADT,TLOADT,TEMPIT,ILATET)
C***********************************************************************
C
C**** THIS ROUTINE STORES THE RESULTS OF THE CURRENT TIME STEP
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
C
      DIMENSION HTLODT(NHLODT,NSUBFT,NFUNCT), IFFIXT(NTOTVT,2), 
     .          RLOADT(NTOTVT),               DISPRT(NTOTVT,3),
     .          DISTOT(NTOTVT,3),             HEADST(NPOINT,4), 
     .          REFORT(NTOTVT),               TLOADT(NTOTVT,2),
     .          ELPRET(NPREVT),               ELVART(NSTATT),
     .          TEMPIT(NTOTVT,2)
C
C=======================================================================
C
C**** RESTART FILE
C
C=======================================================================
C
      IF(NFURES.EQ.1) THEN
       CALL RSTAR4T(TEMPIT,DISTOT,ELPRET,ELVART,TLOADT)
C
       WRITE(LUPRIT,*)
       WRITE(LUPRIT,*) 'FUTURE RESTART INFORMATION'
       WRITE(LUPRIT,*) 'RESULTS STORED FOR TIME:',TTIMET
       WRITE(LUPRIT,*)
       CLOSE(LUPRIT)
       OPEN(UNIT=LUPRIT,FILE=CFT,STATUS='OLD',ACCESS='APPEND')
      ENDIF                      ! nfures.eq.1
C
C**** DECIDE ON SAVING THESE RESULTS
C
      IF(NFURES.EQ.2) THEN
C
       ISAVET=0
       IF(MOD(ISTEPT,NOUTPT(1)).EQ.0) ISAVET=1
       IF(ILATET.EQ. 1)               ISAVET=1
       IF(KSAVET.EQ.-1)               ISAVET=0
C
C**** DUMP TO DATA BASE & SET UP APPROPRIATE COUNTERS
C
       IF(ISAVET.EQ.1) THEN
        IF(KTSTET.EQ.0) KTSTET=1
        IF(KSAVET.EQ.1) KTSTET=KTSTET+1
        CALL RSTAR3T(DISPRT,DISTOT,ELPRET,ELVART,HEADST,HTLODT,IFFIXT,
     .               REFORT,RLOADT,TLOADT,    1)
        IF(KSAVET.EQ.1) NRECCT=NRECCT+NLENCT
       ENDIF
C
      ENDIF                      ! nfures.eq.2
C
C=======================================================================
C
C**** PROCESS FILE
C
C=======================================================================
C
C**** INTERCHANGE CURRENT AND CONVERGED VALUES OF ELPRE & ELVAR
C
      DO INDEXT=1,5
       IDUMMT=IDATPT(2,INDEXT)
       IDATPT(2,INDEXT)=IDATPT(4,INDEXT)
       IDATPT(4,INDEXT)=IDUMMT
       IDUMMT=IDATPT(3,INDEXT)
       IDATPT(3,INDEXT)=IDATPT(5,INDEXT)
       IDATPT(5,INDEXT)=IDUMMT
      ENDDO
C
C**** WRITE DISTO TO PROCESS FILE
C
      CALL DATBAST(DISTOT,  12,    1)
C
C=======================================================================
C
C**** RESULTS FILE
C
C=======================================================================
C
      WRITE(LUPRIT,900)
  900 FORMAT(50X,'>>>>>>> STEP WRITTEN UP',/)
ctm
ctm   cambio en el sistema operativo del CONVEX (febrero de 1993)
ctm
c      CLOSE(LUREST)
c      OPEN(UNIT=LUREST,STATUS='OLD',ACCESS='APPEND')
ctm
C
C**** INITIALISE FLAGS IREST,ISKIP
C
      IRESTT=0
      ISKIPT=0
C
      RETURN
      END
