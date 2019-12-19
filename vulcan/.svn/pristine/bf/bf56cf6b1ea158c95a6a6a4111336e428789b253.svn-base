      SUBROUTINE PCGMUL(A,F,W,N)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N,N), F(N), W(N)
C
C$DIR SCALAR
      DO I=1,N
        W(I)=0.0
C$DIR SCALAR
        DO J=1,N
          W(I)=W(I)+A(I,J)*F(J)
        ENDDO
        F(I)=W(I)
      ENDDO
C
      RETURN
      END
