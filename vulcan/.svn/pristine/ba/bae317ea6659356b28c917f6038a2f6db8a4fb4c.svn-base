      SUBROUTINE FUNLO5(FACTI,HTLOD,TTIME)
C***********************************************************************
C
C**** THIS ROUTINE MODELS AN EXPONENTIAL FUNCTION
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION HTLOD(*)
C
      FACTI=0.0
      IF(TTIME.LE.HTLOD(2)) RETURN       ! NOT STARTED YET
C
      RTIME=TTIME-HTLOD(2)
      FACTI=HTLOD(3)+(HTLOD(4)-HTLOD(3))*EXP(-HTLOD(5)*RTIME)-HTLOD(4)
C      
      RETURN
      END
