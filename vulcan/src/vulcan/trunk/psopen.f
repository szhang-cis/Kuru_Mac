      SUBROUTINE PSOPEN
C***********************************************************************
C
C**** THIS ROUTINE OPENS THE FOLLOWING FILE:
C
C     - FILE 10  POSPROCESS FILE 
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'auxl_om.f'
C
      IF(KPOST.EQ.0) THEN
        NWPOS=0
        RETURN
      ENDIF
C
      NWPOS=1
C
      CALL CPUTIM(TIME1)
C
      IF(IREST.NE.0)  GOTO 100
C
C**** N E W   R U N  *** OPEN NEW POSTPROCESS FILE  . . . . . 
C
      IEROR=3
      IF(KPOST.EQ.1) THEN
C
       if(nmachim.eq.1) then
        OPEN(UNIT=LUPOS,FILE=CJ,STATUS='NEW',FORM='UNFORMATTED')
       else
        OPEN(UNIT=LUPOS,FILE=CJ,STATUS='UNKNOWN',FORM='UNFORMATTED')
       endif
       CLOSE(LUPOS)
      ENDIF
      IEROR=0
      GOTO 1000
C
C**** R E S T A R T  R U N ***   OPEN OLD POSTPROCESS FILE . . . . .
C
 100  CONTINUE
      IEROR=3
      IF(KPOST.EQ.1) THEN
        OPEN(UNIT=LUPOS,FILE=CJ,STATUS='OLD',FORM='UNFORMATTED',
     .       ACCESS='APPEND',ERR=101)
        CLOSE(LUPOS)
        NWPOS=0
      ENDIF
      IEROR=0
      GO TO 1000
C
  101 CONTINUE
      OPEN(UNIT=LUPOS,FILE=CJ,STATUS='NEW',FORM='UNFORMATTED')
      CLOSE(LUPOS)
      IEROR=0
C
 1000 CONTINUE
C
      CALL CPUTIM(TIME2)
      CPURS=CPURS+(TIME2-TIME1)
C
      IF(IEROR.NE.0)THEN
        IF(IEROR.EQ.3) WRITE(LURES,903)
        CALL RUNEND('  ERROR IN OPENING FILES           ')
      ENDIF
C
      RETURN
C
  903 FORMAT(' ERROR IN OPENING POSTPR. FILE 10')     
      END
