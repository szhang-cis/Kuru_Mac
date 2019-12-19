
      SUBROUTINE FUNLO2(FACTI,HTLOD,TTIME)
C***********************************************************************
C
C****THIS ROUTINE MODELS PARABOLIC RISING + CONSTANT PLATEAU
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION HTLOD(*)
C
      FACTI=0.0
      IF(TTIME.LE.HTLOD(2)) RETURN       ! NOT STARTED YET
      IF(TTIME.GT.HTLOD(3)) RETURN       ! ALREADY FINISHED
      IF(TTIME.GT.HTLOD(4)) THEN         ! TIME RISE FINISHED
        FACTI=HTLOD(5)
      ELSE                               ! STILL RISING
        RTIME=TTIME-HTLOD(2)
        CONST=HTLOD(4)-HTLOD(2)
        FACTI=RTIME*(2.*CONST-RTIME)*(HTLOD(5)/CONST**2)
      ENDIF
C      
      RETURN
      END
