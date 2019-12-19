      SUBROUTINE PRIRCG(NTYPE,RCGTT,ALAMB,ILAEQ)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE PRINCIPAL STRETCHES OF THE DEFORMATION
C     GRADIENT F BY MEANS OF THE EIGENVALUES OF THE RIGHT CAUCHY TENSOR
C     C = F^T : F (THE LAST ARE THE SQUARE OF THE FORMER)
C
C     RCGTT = ( C11 , C22 , C12 , C33 , C13 , C23 )
C
C     ILAEQ: index for equal principal stretches
C            1=three different principal stretches
C            2=two equal principal stretches
C            3=three equal principal stretches
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION RCGTT(*), ALAMB(*)
      DIMENSION X(6), A(6)
C
      A(1)=1.0D0              ! A = Kronecker Delta
      A(2)=1.0D0
      A(3)=0.0D0
      A(4)=1.0D0
      A(5)=0.0D0
      A(6)=0.0D0
      TP=6.283185307179586D0
      S3=1.732050807568877D0
      ZM=1.0D-06              ! the value 1.0D-10 is too strict
C
      IF(NTYPE.EQ.5) THEN     ! 1-D
       ALAMB(1)=RCGTT(1)
       ALAMB(2)=1.0D0
       ALAMB(3)=1.0D0
C
       ILAEQ=2
       IF(DABS(ALAMB(1)-ALAMB(2)).LT.ZM) ILAEQ=3
C
       IF(ALAMB(1).LE.ZM)
     .  CALL RUNEND('ERROR: ZERO OR NEGATIVE PRINCIPAL STRETCH')
C
       ALAMB(1)=DSQRT(ALAMB(1))
       RETURN
      ENDIF
C
      IF(NTYPE.EQ.4) THEN     ! 3-D
       DO I=1,6
        X(I)=0.0D0
       ENDDO
C
       TP3=TP/3.0D0
C
       XJ1=RCGTT(1)+RCGTT(2)+RCGTT(4)
       SM=XJ1/3.0D0
C
       DO I=1,6
        X(I)=RCGTT(I)-A(I)*SM
       ENDDO
C
       XJ2=(X(1)*X(1)+X(2)*X(2)+X(4)*X(4))*0.5D00+
     .      X(3)*X(3)+X(5)*X(5)+X(6)*X(6)
C
       IF(XJ2.GT.0.0D0) THEN
        R=2.0D0*DSQRT(XJ2/3.0D0)
        XJ3=X(1)*(X(2)*X(4)-X(6)*X(6))-
     .      X(3)*(X(3)*X(4)-X(5)*X(6))+
     .      X(5)*(X(3)*X(6)-X(2)*X(5))
        ARG= -((3.0D0)*XJ3)/(R*XJ2)
        IF(ARG.GT.1.0D0) THEN
         ARG=1.0D0
        ELSE IF(ARG.LT.-1.0D0) THEN
         ARG=-1.0D0
        ENDIF
        T=DASIN(ARG)/3.0D0
        ALAMB(1)=R*DSIN(T+      TP3)+SM
        ALAMB(2)=R*DSIN(T          )+SM
        ALAMB(3)=R*DSIN(T+2.0D0*TP3)+SM
       ELSE
        ALAMB(1)=RCGTT(1)
        ALAMB(2)=RCGTT(2)
        ALAMB(3)=RCGTT(4)
       END IF
      ELSE                    ! 2-D: plane stress/strain or axisymmetric
       C11=RCGTT(1)
       C22=RCGTT(2)
       C12=RCGTT(3)
       ALAMB(1)=(C11+C22)/2.0+DSQRT(0.25D0*(C11-C22)*(C11-C22)+C12*C12)
       ALAMB(2)=(C11+C22)/2.0-DSQRT(0.25D0*(C11-C22)*(C11-C22)+C12*C12)
       ALAMB(3)=RCGTT(4)
      END IF
C
      IF(ALAMB(1).LE.ZM.OR.ALAMB(2).LE.ZM.OR.ALAMB(3).LE.ZM)
     . CALL RUNEND('ERROR: ZERO OR NEGATIVE PRINCIPAL STRETCH')
C
      ILAEQ=1
      IF(DABS(ALAMB(1)-ALAMB(2)).LT.ZM) ILAEQ=ILAEQ+1
      IF(DABS(ALAMB(1)-ALAMB(3)).LT.ZM) ILAEQ=ILAEQ+1
      IF(DABS(ALAMB(2)-ALAMB(3)).LT.ZM) ILAEQ=ILAEQ+1
      IF(ILAEQ.EQ.4) ILAEQ=3
C
      IF(ILAEQ.EQ.2) THEN ! ordered in such a way that lambda_2=lambda_3
       A1=ALAMB(1)
       A2=ALAMB(2)
       A3=ALAMB(3)
       IF(DABS(ALAMB(1)-ALAMB(2)).LT.ZM) THEN
        A1=ALAMB(3)
        A2=ALAMB(1)
        A3=ALAMB(2)
       ENDIF
       IF(DABS(ALAMB(1)-ALAMB(3)).LT.ZM) THEN
        A1=ALAMB(2)
        A2=ALAMB(1)
        A3=ALAMB(3)
       ENDIF
       ALAMB(1)=A1
       ALAMB(2)=A2
       ALAMB(3)=A3
      ENDIF
C
      ALAMB(1)=DSQRT(ALAMB(1))
      ALAMB(2)=DSQRT(ALAMB(2))
      ALAMB(3)=DSQRT(ALAMB(3))
C
      RETURN
      END
