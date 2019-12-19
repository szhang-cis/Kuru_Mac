      SUBROUTINE PROMA1(A,B,C,N1,N2,N3)
C********************************************************************
C
C***THIS ROUTINE EVALUATES A MATRIX PRODUCT
C                                
C           A(I,J) = B(I,K) * C(K,J)  
C
C********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N1,N2), B(N1,N3), C(N3,N2)
C
      DO 10 I=1,N1
      DO 10 J=1,N2
        A(I,J)=0.0
      DO 10 K=1,N3
        A(I,J)=A(I,J)+B(I,K)*C(K,J)
   10 CONTINUE
C
      RETURN
      END
