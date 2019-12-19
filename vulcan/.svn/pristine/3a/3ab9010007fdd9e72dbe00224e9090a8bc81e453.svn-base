

      SUBROUTINE VECPR0(A,B,C)
C***********************************************************************
C
C****THIS ROUTINE CALCULATES  AN ORTHONORMAL BASE FROM A
C
C            INPUT :      A(3)        BASE VECTOR
C            OUTPUT:      B(3),C(3)   VECTORS NORMAL TO A
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(3),B(3),C(3)
      DATA TOL/1.0D-07/
C
      XN=A(2)*A(2)+A(3)*A(3)
      IF(XN.GT.TOL) THEN
         B(1)=0.0
         XN=SQRT(XN)
         B(2)= A(3)/XN
         B(3)=-A(2)/XN
         C(1)= A(2)*B(3)-A(3)*B(2)
         C(2)=-A(1)*B(3)
         C(3)= A(1)*B(2)
      ELSE
         B(1)=0.0
         B(2)=1.0
         B(3)=0.0
         C(1)=0.0
         C(2)=0.0
         C(3)=1.0
      END IF
      RETURN
      END
