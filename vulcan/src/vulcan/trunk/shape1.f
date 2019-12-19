      SUBROUTINE SHAPE1(DERIV,S,NDIME,NNODE,NNODX,SHAPE)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES SHAPE FUNCTIONS AND THEIR DERIVATIVES FOR
C     LINEAR AND QUADRATIC ISOPARAMETRIC 1-D ELEMENTS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DERIV(NDIME,*),SHAPE(*)
C
C**** 2 NODED ELEMENT
C
      IF(NNODE.EQ.2) THEN
       SHAPE(1+NNODX)=0.5D0*(1.0D0-S)
       SHAPE(2+NNODX)=0.5D0*(1.0D0+S)
C
       DERIV(1,1+NNODX)=-0.5D0
       DERIV(1,2+NNODX)= 0.5D0
C
       RETURN
      ENDIF
C
C**** 3 NODED ELEMENT
C
      IF(NNODE.EQ.3) THEN
       SHAPE(1)=S*(-1.0D0+S)*0.5D0
       SHAPE(2)=S*( 1.0D0+S)*0.5D0
       SHAPE(3)=1.0D0-S*S
C
       DERIV(1,1)= -0.5D0+S
       DERIV(1,2)=  0.5D0+S
       DERIV(1,3)= -2.0D0*S
C
       RETURN
      ENDIF
C
      END
