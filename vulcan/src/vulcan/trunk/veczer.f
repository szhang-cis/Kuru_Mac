      SUBROUTINE VECZER(N,V)
C***********************************************************************
C
C**** THIS ROUTINE INITIALISES ONE VECTOR TO ZERO
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION V(N)
C
      DO I=1,N
       V(I)=0.0D0
      END DO
      RETURN
      END
