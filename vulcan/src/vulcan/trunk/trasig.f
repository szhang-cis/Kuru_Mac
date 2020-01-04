
      SUBROUTINE TRASIG(S,T,NS,ND,FLAG)
C************************************************************************
C
C***ROUTINE TO TRANSFORM STRESS TENSOR 
C
C************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION S(NS), T(ND,ND)
      DIMENSION D(6,6)
      CHARACTER*6 FLAG
C
C*** 2D CASE
C
      IF(ND.EQ.2) THEN
C
      IF(FLAG.EQ.'G_TO_L') THEN               ! GLOBAL TO LOCAL
C
      D(1,1)=S(1)*T(1,1)+S(3)*T(1,2)
      D(2,1)=S(3)*T(1,1)+S(2)*T(1,2)
C
      D(1,2)=S(1)*T(2,1)+S(3)*T(2,2)
      D(2,2)=S(3)*T(2,1)+S(2)*T(2,2)
C
      S(1)=T(1,1)*D(1,1)+T(1,2)*D(2,1)
      S(2)=T(2,1)*D(1,2)+T(2,2)*D(2,2)
      S(3)=T(1,1)*D(1,2)+T(1,2)*D(2,2)
C
      ELSE                           ! LOCAL TO GLOBAL
C
      D(1,1)=S(1)*T(1,1)+S(3)*T(2,1)
      D(2,1)=S(3)*T(1,1)+S(2)*T(2,1)
C
      D(1,2)=S(1)*T(1,2)+S(3)*T(2,2)
      D(2,2)=S(3)*T(1,2)+S(2)*T(2,2)
C
      S(1)=T(1,1)*D(1,1)+T(2,1)*D(2,1)
      S(2)=T(1,2)*D(1,2)+T(2,2)*D(2,2)
      S(3)=T(1,1)*D(1,2)+T(2,1)*D(2,2)
C
      ENDIF
C
C
      ELSE
C
C*** 3D CASE
C
      IF(FLAG.EQ.'G_TO_L') THEN               ! GLOBAL TO LOCAL
C
      D(1,1)=S(1)*T(1,1)+S(3)*T(1,2)+S(5)*T(1,3)
      D(2,1)=S(3)*T(1,1)+S(2)*T(1,2)+S(6)*T(1,3)
      D(3,1)=S(5)*T(1,1)+S(6)*T(1,2)+S(4)*T(1,3)
C
      D(1,2)=S(1)*T(2,1)+S(3)*T(2,2)+S(5)*T(2,3)
      D(2,2)=S(3)*T(2,1)+S(2)*T(2,2)+S(6)*T(2,3)
      D(3,2)=S(5)*T(2,1)+S(6)*T(2,2)+S(4)*T(2,3)
C
      D(1,3)=S(1)*T(3,1)+S(3)*T(3,2)+S(5)*T(3,3)
      D(2,3)=S(3)*T(3,1)+S(2)*T(3,2)+S(6)*T(3,3)
      D(3,3)=S(5)*T(3,1)+S(6)*T(3,2)+S(4)*T(3,3)
C
      S(1)=T(1,1)*D(1,1)+T(1,2)*D(2,1)+T(1,3)*D(3,1)
      S(2)=T(2,1)*D(1,2)+T(2,2)*D(2,2)+T(2,3)*D(3,2)
      S(3)=T(1,1)*D(1,2)+T(1,2)*D(2,2)+T(1,3)*D(3,2)
      S(4)=T(3,1)*D(1,3)+T(3,2)*D(2,3)+T(3,3)*D(3,3)
      S(5)=T(1,1)*D(1,3)+T(1,2)*D(2,3)+T(1,3)*D(3,3)
      S(6)=T(2,1)*D(1,3)+T(2,2)*D(2,3)+T(2,3)*D(3,3)
C
      ELSE                           ! LOCAL TO GLOBAL
C
      D(1,1)=S(1)*T(1,1)+S(3)*T(2,1)+S(5)*T(3,1)
      D(2,1)=S(3)*T(1,1)+S(2)*T(2,1)+S(6)*T(3,1)
      D(3,1)=S(5)*T(1,1)+S(6)*T(2,1)+S(4)*T(3,1)
C
      D(1,2)=S(1)*T(1,2)+S(3)*T(2,2)+S(5)*T(3,2)
      D(2,2)=S(3)*T(1,2)+S(2)*T(2,2)+S(6)*T(3,2)
      D(3,2)=S(5)*T(1,2)+S(6)*T(2,2)+S(4)*T(3,2)
C
      D(1,3)=S(1)*T(1,3)+S(3)*T(2,3)+S(5)*T(3,3)
      D(2,3)=S(3)*T(1,3)+S(2)*T(2,3)+S(6)*T(3,3)
      D(3,3)=S(5)*T(1,3)+S(6)*T(2,3)+S(4)*T(3,3)
C
      S(1)=T(1,1)*D(1,1)+T(2,1)*D(2,1)+T(3,1)*D(3,1)
      S(2)=T(1,2)*D(1,2)+T(2,2)*D(2,2)+T(3,2)*D(3,2)
      S(3)=T(1,1)*D(1,2)+T(2,1)*D(2,2)+T(3,1)*D(3,2)
      S(4)=T(1,3)*D(1,3)+T(2,3)*D(2,3)+T(3,3)*D(3,3)
      S(5)=T(1,1)*D(1,3)+T(2,1)*D(2,3)+T(3,1)*D(3,3)
      S(6)=T(1,2)*D(1,3)+T(2,2)*D(2,3)+T(3,2)*D(3,3)
C
      ENDIF
C
      ENDIF
C
      RETURN
      END