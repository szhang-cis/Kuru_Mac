      SUBROUTINE VECADD(N,V1,V2,V3)
C***********************************************************************
C
C**** VECTOR ADITION:    V3(I) = V2(I) + V1(I)   I=1..N
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION V1(N),V2(N),V3(N)
C
      DO I=1,N
       V3(I) = V2(I) + V1(I)
      END DO
C
      RETURN
      END
