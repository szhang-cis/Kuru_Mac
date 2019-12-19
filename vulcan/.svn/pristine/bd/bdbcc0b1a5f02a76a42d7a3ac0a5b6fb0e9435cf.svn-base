      SUBROUTINE FUNLO3(FACTI,HTLOD,TTIME)
C***********************************************************************
C
C**** THIS ROUTINE MODELS SINUSOIDAL WAVES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION HTLOD(*)
C
      TWOPI=6.283185307179586
C
      FACTI=0.0
      IF(TTIME.LE.HTLOD(2)) RETURN       ! NOT STARTED YET
      IF(TTIME.GT.HTLOD(3)) RETURN       ! ALREADY FINISHED
C
      FACTI=HTLOD(5)*DSIN((TWOPI/HTLOD(4))*(TTIME-HTLOD(2)))
C      
      RETURN
      END
