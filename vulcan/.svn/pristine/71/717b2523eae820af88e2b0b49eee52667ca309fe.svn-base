      SUBROUTINE PDIREC(S,XL,V)
C***********************************************************************
C
C****THIS ROUTINE COMPUTE THE EIGENVECTORS
C
C         INPUT :     XL (3)   EIGENVALUES
C                     S  (6)   ORIGINAL MATRIX
C         OUTPUT:     V(3,3)   EIGENVECTORS (NORMALISED)
C
C    N.B.    S  = ( Sxx , Syy , Sxy , Szz , Sxz , Syz )
C            XL = (  S1 ,  S2 ,  S3 )
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION S(6),XL(3),V(3,3)
      DATA TOL/1.0D-07/
C
      COMP=XL(1)+XL(2)+XL(3)
      COMP2=COMP*COMP
      A3=ABS((XL(1)-XL(2))/COMP)
      A2=ABS((XL(1)-XL(3))/COMP)
      A1=ABS((XL(2)-XL(3))/COMP)
C
C      CASE-A: THREE EQUAL EIGENVALUES
C
      IF(A1.LT.TOL.AND.A2.LT.TOL) THEN
         V(1,1)=1.0
         V(2,2)=1.0
         V(3,3)=1.0
         V(1,2)=0.0
         V(1,3)=0.0
         V(2,1)=0.0
         V(2,3)=0.0
         V(3,1)=0.0
         V(3,2)=0.0
         RETURN
      END IF
C
C      CASE-B: TWO EQUAL EIGENVALUES
C
      IF(A1.LT.TOL) THEN
         CALL EIGVC1(S,XL(1),V(1,1),COMP2)
         CALL VECPR0(V(1,1),V(1,2),V(1,3))
         RETURN
      END IF
      IF(A2.LT.TOL) THEN
         CALL EIGVC1(S,XL(2),V(1,2),COMP2)
         CALL VECPR0(V(1,2),V(1,3),V(1,1))
         RETURN
      END IF
      IF(A3.LT.TOL) THEN
         CALL EIGVC1(S,XL(3),V(1,3),COMP2)
         CALL VECPR0(V(1,3),V(1,1),V(1,2))
         RETURN
      END IF
C
C      CASE-C: THREE DIFFERENT EIGENVALUES
C
      CALL EIGVC1(S,XL(1),V(1,1),COMP2)
      CALL EIGVC1(S,XL(2),V(1,2),COMP2)
      DOT=V(1,1)*V(1,2)+V(2,1)*V(2,2)+V(3,1)*V(3,2)
      V(1,2)=V(1,2)-DOT*V(1,1)
      V(2,2)=V(2,2)-DOT*V(2,1)
      V(3,2)=V(3,2)-DOT*V(3,1)
      XN=SQRT(V(1,2)*V(1,2)+V(2,2)*V(2,2)+V(3,2)*V(3,2))
      XN=1./XN
      V(1,2)=V(1,2)*XN
      V(2,2)=V(2,2)*XN
      V(3,2)=V(3,2)*XN
      V(1,3)=V(2,1)*V(3,2)-V(3,1)*V(2,2)
      V(2,3)=V(3,1)*V(1,2)-V(1,1)*V(3,2)
      V(3,3)=V(1,1)*V(2,2)-V(2,1)*V(1,2)
      RETURN
      END
