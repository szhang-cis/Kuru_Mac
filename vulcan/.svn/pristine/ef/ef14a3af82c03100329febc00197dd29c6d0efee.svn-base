      SUBROUTINE EFFJNT(SGTOT,YIELD,NSTRE)
C************************************************************************
C
C****ROUTINE TO EVALUATE THE EQUIVALENT EFFECTIVE STRESS FOR A
C    DRUCKER-PRAGER FRICTIONAL MATERIAL
C
C************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SGTOT(*)
C
      SHEA1=SGTOT(3)
      SHEA2=0.0
      IF(NSTRE.GT.4) SHEA2=SGTOT(5)
      YIELD=SQRT(SHEA1**2+SHEA2**2)
C
      RETURN
      END
