      SUBROUTINE PRINST(S,D)
C***********************************************************************
C
C****THIS ROUTINE CALCULATES THE EIGEN-VALUES OF S(6) ( STRAINS )  
C    AND PUTS THE RESULTS IN D(3)
C
C    N.B.    S = ( Sxx , Syy , Sxy , Szz , Sxz , Syz )
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION S(*),D(*)
      DIMENSION X(6), A(6), B(6)
      DATA A/1.0D00,1.0D00,0.0D00,1.0D00,0.0D00,0.0D00/
      DATA B/1.0D00,1.0D00,0.5D00,1.0D00,0.5D00,0.5D00/
      DATA TP,S3,ZM/6.283185307179586,1.732050807568877,1.0E-10/
C
      TP3=TP/3.0D00
C
      XJ1=S(1)+S(2)+S(4)
      SM=XJ1/3.0D00
C
      DO I=1,6
        X(I)=B(I)*S(I)-A(I)*SM
      ENDDO
C
      XJ2=(X(1)**2+X(2)**2+X(4)**2)*0.5D00+
     &     X(3)**2+X(5)**2+X(6)**2
C
      IF(XJ2.GT.0.0D00) THEN
        R=2.0D00*DSQRT(XJ2/3.0D00)
        XJ3=X(1)*(X(2)*X(4)-X(6)*X(6))-
     &      X(3)*(X(3)*X(4)-X(5)*X(6))+
     &      X(5)*(X(3)*X(6)-X(2)*X(5))
        ARG= -((3.00D00)*XJ3)/(R*XJ2)
        IF(ARG.GT.1.0) THEN
          ARG=1.0
        ELSE IF(ARG.LT.-1.0) THEN
          ARG=-1.0
        ENDIF
        T=DASIN(ARG)/3.0D00
        D(1)=R*DSIN(T+       TP3)+SM
        D(2)=R*DSIN(T           )+SM
        D(3)=R*DSIN(T+2.0D00*TP3)+SM
      ELSE
        D(1)=S(1)
        D(2)=S(2)
        D(3)=S(4)
      END IF
C
      RETURN
      END
