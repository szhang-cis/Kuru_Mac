      SUBROUTINE VECPRO(N,V1,V2,V3)
C***********************************************************************
C
C**** TRIDIMENSIONAL VECTORIAL PRODUCT OF TWO VECTORS  V3 = V1 X V2
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION V1(N),V2(N),V3(N)
C
      V3(1) = V1(2)*V2(3) - V1(3)*V2(2)
      V3(2) = V1(3)*V2(1) - V1(1)*V2(3)
      V3(3) = V1(1)*V2(2) - V1(2)*V2(1)
      RETURN
      END
