      SUBROUTINE FUNLO8(FACTI,HTLOD,TTIME)
C***********************************************************************
C
C**** THIS ROUTINE MODELS  LINEAR RISING + CONSTANT PLATEAU + LINEAR
C     FALLING (EQUAL RISING & FALLING TIME INTERVAL ARE ASSUMED)
C
C         . . . .
C       .         .
C     .             .
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION HTLOD(*)
C
      FACTI=0.0
      IF(TTIME.LE.HTLOD(2)) RETURN       ! NOT STARTED YET
      IF(TTIME.GT.HTLOD(3)) RETURN       ! ALREADY FINISHED
C
      IF(TTIME.GT.HTLOD(2).AND.TTIME.LE.HTLOD(4))        ! RISING
     . FACTI=(TTIME-HTLOD(2))*(HTLOD(5)/(HTLOD(4)-HTLOD(2)))
C
      TIMEX=HTLOD(3)-HTLOD(4)+HTLOD(2)   ! staring falling time
      IF(TTIME.GT.HTLOD(4).AND.TTIME.LE.TIMEX)           ! PLATEAU
     . FACTI=HTLOD(5)
C
      IF(TTIME.GT.TIMEX.AND.TTIME.LE.HTLOD(3))           ! FALLING
     . FACTI=-(TTIME-TIMEX)*(HTLOD(5)/(HTLOD(3)-TIMEX))+HTLOD(5)
C
      RETURN
      END
