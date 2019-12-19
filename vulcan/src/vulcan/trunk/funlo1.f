      SUBROUTINE FUNLO1(FACTI,HTLOD,TTIME)
C***********************************************************************
C
C**** THIS ROUTINE MODELS LINEAR RISING + CONSTANT PLATEAU
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION HTLOD(*)
C
      FACTI=0.0
      IF(TTIME.LE.HTLOD(2)) RETURN       ! NOT STARTED YET
      IF(TTIME.GT.HTLOD(3)) RETURN       ! ALREADY FINISHED
      IF(TTIME.GT.HTLOD(4)) THEN         ! TIME RISE FINISHED
       FACTI=HTLOD(5)
      ELSE                               ! STILL RISING
       FACTI=(TTIME-HTLOD(2))*(HTLOD(5)/(HTLOD(4)-HTLOD(2)))
      ENDIF
C      
      RETURN
      END
