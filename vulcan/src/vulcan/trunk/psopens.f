      SUBROUTINE PSOPENS
C***********************************************************************
C
C**** THIS ROUTINE OPENS THE FOLLOWING FILE:
C
C     - FILE 210 (20 LINUX)  POSPROCESS FILE
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_oms.f'
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'auxl_oms.f'
C
      IF(KPOSTS.EQ.0) THEN
        NWPOSS=0
        RETURN
      ENDIF
C
      NWPOSS=1
C
      CALL CPUTIMS(TIME1S)
C
      IF(IRESTS.NE.0)  GOTO 100
C
C**** N E W   R U N  *** OPEN NEW POSTPROCESS FILE  . . . . . 
C
      IERORS=3
      IF(KPOSTS.EQ.1) THEN
C
       if(nmachis.eq.1) then
        OPEN(UNIT=LUPOSS,FILE=CJS,STATUS='NEW',FORM='UNFORMATTED')
       else
        OPEN(UNIT=LUPOSS,FILE=CJS,STATUS='UNKNOWN',FORM='UNFORMATTED')
       endif
       CLOSE(LUPOSS)
      ENDIF
      IERORS=0
      GOTO 1000
C
C**** R E S T A R T  R U N ***   OPEN OLD POSTPROCESS FILE . . . . .
C
 100  CONTINUE
      IERORS=3
      IF(KPOSTS.EQ.1) THEN
        OPEN(UNIT=LUPOSS,FILE=CJS,STATUS='OLD',FORM='UNFORMATTED',
     .       ACCESS='APPEND',ERR=101)
        CLOSE(LUPOSS)
        NWPOSS=0
      ENDIF
      IERORS=0
      GO TO 1000
C
  101 CONTINUE
      OPEN(UNIT=LUPOSS,FILE=CJS,STATUS='NEW',FORM='UNFORMATTED')
      CLOSE(LUPOSS)
      IERORS=0
C
 1000 CONTINUE
C
      CALL CPUTIMS(TIME2S)
      CPURSS=CPURSS+(TIME2S-TIME1S)
C
      IF(IERORS.NE.0)THEN
        IF(IERORS.EQ.3) WRITE(LURESS,903)
        CALL RUNENDS('  ERROR IN OPENING FILES           ')
      ENDIF
C
      RETURN
C
  903 FORMAT(' ERROR IN OPENING POSTPR. FILE 210 (30 LINUX)')     
      END
