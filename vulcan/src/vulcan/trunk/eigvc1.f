
      SUBROUTINE EIGVC1(S,XL,V,COM)
C***********************************************************************
C
C****THIS ROUTINE CALCULATES ONE EIGENVECTOR
C
C             INPUT:    XL     EIGENVALUE
C                       S (6)  ORIGINAL MATRIX
C             OUTPUT:   V (3)  EIGENVECTOR
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION S(6),V(3)
      DATA TOL/1.0D-07/
C
      X1=S(1)-XL
      X2=S(2)-XL
      X3=S(4)-XL
      X4=S(3)
      X5=S(6)
      X6=S(5)
C
      DET=X1*X2-X4*X4
      IF(ABS(DET/COM).GT.TOL) THEN
         V1=(X4*X5-X6*X2)/DET
         V2=(X6*X4-X5*X1)/DET
         XN=SQRT(V1*V1+V2*V2+1.0)
         XN=1./XN
         V(1)=V1*XN
         V(2)=V2*XN
         V(3)=XN
         RETURN
      END IF
C
      DET=X4*X5-X6*X2
      IF(ABS(DET/COM).GT.TOL) THEN
         V1=(X6*X4-X5*X1)/DET
         XN=SQRT(V1*V1+1.0)
         XN=1./XN
         V(1)=XN
         V(2)=V1*XN
         V(3)=0.0
         RETURN
      END IF
C
      DET=X1*X5-X4*X6
      IF(ABS(DET/COM).GT.TOL) THEN
         V(1)=0.0
         V(2)=1.0
         V(3)=0.0
         RETURN
      END IF
C
      DET=X4*X3-X5*X6
      IF(ABS(DET/COM).GT.TOL) THEN
         V1=(X5*X5-X2*X3)/DET
         XN=SQRT(V1*V1+1.0)
         XN=1./XN
         V(1)=V1*XN
         V(2)=XN
         V(3)=0.0
         RETURN
      END IF
C
      DET=X2*X3-X5*X5
      IF(ABS(DET/COM).GT.TOL) THEN
         V(1)=1.0
         V(2)=0.0
         V(3)=0.0
         RETURN
      END IF
C
      V(1)=0.0
      V(2)=1.0
      V(3)=0.0
      RETURN
      END
