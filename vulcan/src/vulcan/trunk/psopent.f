      SUBROUTINE PSOPENT
C***********************************************************************
C
C**** THIS ROUTINE OPENS THE FOLLOWING FILE:
C
C     - FILE 210 (30 LINUX)  POSPROCESS FILE 
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
      INCLUDE 'auxl_omt.f'
C
      IF(KPOSTT.EQ.0) THEN
        NWPOST=0
        RETURN
      ENDIF
C
      NWPOST=1
C
      CALL CPUTIMT(TIME1T)
C
      IF(IRESTT.NE.0)  GOTO 100
C
C**** N E W   R U N  *** OPEN NEW POSTPROCESS FILE  . . . . . 
C
      IERORT=3
      IF(KPOSTT.EQ.1) THEN
C
       if(nmachi.eq.1) then
        OPEN(UNIT=LUPOST,FILE=CJT,STATUS='NEW',FORM='UNFORMATTED')
       else
        OPEN(UNIT=LUPOST,FILE=CJT,STATUS='UNKNOWN',FORM='UNFORMATTED')
       endif
       CLOSE(LUPOST)
      ENDIF
      IERORT=0
      GOTO 1000
C
C**** R E S T A R T  R U N ***   OPEN OLD POSTPROCESS FILE . . . . .
C
 100  CONTINUE
      IERORT=3
      IF(KPOSTT.EQ.1) THEN
        OPEN(UNIT=LUPOST,FILE=CJT,STATUS='OLD',FORM='UNFORMATTED',
     .       ACCESS='APPEND',ERR=101)
        CLOSE(LUPOST)
        NWPOST=0
      ENDIF
      IERORT=0
      GO TO 1000
C
  101 CONTINUE
      OPEN(UNIT=LUPOST,FILE=CJT,STATUS='NEW',FORM='UNFORMATTED')
      CLOSE(LUPOST)
      IERORT=0
C
 1000 CONTINUE
C
      CALL CPUTIMT(TIME2T)
      CPURST=CPURST+(TIME2T-TIME1T)
C
      IF(IERORT.NE.0)THEN
        IF(IERORT.EQ.3) WRITE(LUREST,903)
        CALL RUNENDT('  ERROR IN OPENING FILES           ')
      ENDIF
C
      RETURN
C
  903 FORMAT(' ERROR IN OPENING POSTPR. FILE 210 (30 LINUX)')     
      END
