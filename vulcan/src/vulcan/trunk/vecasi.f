      SUBROUTINE VECASI(N,V1,V2)
C***********************************************************************
C
C**** VECTOR ASSIGN:    V2(I) = V1(I)   I=1..N
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION V1(N),V2(N)
C
      DO I=1,N
       V2(I)=V1(I)
      END DO
      RETURN
      END
