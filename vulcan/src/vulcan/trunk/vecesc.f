      SUBROUTINE VECESC(N,A,V1)
C***********************************************************************
C
C**** THIS ROUTINE ESCALATES VECTOR V1  A * V1(I) -> V1(I)   I=1..N
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION V1(N)
C
      DO I=1,N
       V1(I)=A*V1(I)
      END DO
C
      END
