      SUBROUTINE INVERT(A,NMAX,NDM)
C***********************************************************************
C
C**** THIS ROUTINE PERFORMS THE INVERSION OF A NDM*NDM SQUARE MATRIX 
C     OR JUST PART OF IT (NMAX*NMAX)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(NDM,NDM)
      DO 200 N = 1,NMAX
      D = A(N,N)
      DO 100 J = 1,NMAX
100   A(N,J) = -A(N,J)/D
      DO 150 I = 1,NMAX
      IF(N.EQ.I) GO TO 150
      DO 140 J = 1,NMAX
      IF(N.NE.J) A(I,J) = A(I,J) +A(I,N)*A(N,J)
140   CONTINUE
150   A(I,N) = A(I,N)/D
      A(N,N) = 1.0/D
200   CONTINUE
      RETURN
      END
