      SUBROUTINE INTREDT(TFICTT)
C***********************************************************************
C
C**** THIS ROUTINE SETS UP THE TYPE OF ALGORITHM AND REFERENCE LOAD
C     TO BE USED IN THE CURRENT TIME INTERVAL
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
C
      DIMENSION TFICTT(NFUNCT)
C
      CHARACTER XCOMMT*9,YCOMMT*123
C
C**** NEW RUN: LOOK FOR 'INTE'RVAL CARD OR FOR END OF ANALYSIS ('STOP')
C
      IF(IRESTT.EQ.0) THEN
       CALL INTRD0T(DTIMET,KSTEPT,NSTEPT,NSSTEPT)
       GO TO 1000
      ENDIF
C
C**** RESTART RUN:
C
      IF(IRESTT.EQ.1) THEN
C
C**** FIRST OF ALL READ ALL DATA UP TO THE END OF THE GENERAL DATA
C     AND COUNT THE NUMBER OF "INTERVAL DATA BLOCKS" UNTIL
C     THE ONE TO RESTART IS FOUND
C
       KOUNTT=0
  200  CONTINUE
       READ(LUDATT,900) XCOMMT,YCOMMT
       IF(XCOMMT.NE.'INTERVAL_') GOTO 200
       KOUNTT=KOUNTT+1
       IF(KOUNTT.LT.ITIMET) GOTO 200
C
C**** NOW LOOK FOR THE LAST "END_INTERVAL" CARD
C
  210  CONTINUE
       READ(LUDATT,900) XCOMMT,YCOMMT
       IF(XCOMMT.NE.'END_INTER') GOTO 210
C
C**** IF THE LAST INTERVAL HAS BEEN FINISHED ( NSTEP < KSTEP )
C     IN A "RESTART CONTINUE" THE NEXT CARDS ARE SUPPOSED TO BE 
C     THE DATA FOR A NEW INTERVAL
C
       IF(ISKIPT.EQ.0) THEN
        KSTEPT=ISTEPT+1
        IF(KSTEPT.LE.NSTEPT) GO TO 1000
C
        IRESTT=2
C
        CALL INTRD0T(DTIMET,KSTEPT,NSTEPT,NSSTEPT)
        GO TO 1000
C
       ENDIF
C
C**** IN A "RESTART SKIP" THE NEXT CARDS ARE SUPPOSED TO BE 
C     THE NEW STRATEGY FOR THE NEXT INTERVAL 
C
       IF(ISKIPT.EQ.1)THEN
C
        CALL INTRD0T(DTIMET,KSTEPT,NSTEPT,NSSTEPT)
        GO TO 1000
C
       ENDIF
C
       GOTO 1000
C
      ENDIF
C
C**** END OF SEARCH
C
 1000 CONTINUE
C
C**** IF A NEW INTERVAL IS STARTING UPDATE SOME VALUES
C
      IF(KSTEPT.EQ.1) THEN
       ITIMET=ITIMET+1
       DO IFUNCT=1,NFUNCT
        TFICTT(IFUNCT)=0.0
       ENDDO 
       KUNLDT=0
      ENDIF
C
C**** SOME CONTROL (always coupled problems)
C
      IF(NSSTEPT.GT.1.AND.NFURES.GT.0)
     . CALL RUNMENT('WARNING: INCREMENTATION AND FUTURE ANALYSIS COULD
     . NOT BE COMPATIBLE')
C
      RETURN
  900 FORMAT(A9,A123)
      END