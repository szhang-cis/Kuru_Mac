      SUBROUTINE YLDJNT(PROPS,SGTOT,COMPA)
C************************************************************************
C
C****ROUTINE TO EVALUATE THE YIELD SHEAR STRESS FOR A
C    DRUCKER-PRAGER FRICTIONAL MATERIAL
C
C************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PROPS(*), SGTOT(*)
      DATA RADIA/1.7453292519943296D-02/
C
      SNORM=SGTOT(1)
      COHES=PROPS(13)
      FRICT=PROPS(14)*RADIA
C
      COMPA=COHES-SNORM*DTAN(FRICT)
C
      RETURN
      END
