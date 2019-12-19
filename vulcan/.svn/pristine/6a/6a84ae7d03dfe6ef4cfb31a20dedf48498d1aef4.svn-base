      SUBROUTINE TRANS3(B,T,IFLA1,IFLA2)
C**************************************************************
C
C****THIS ROUTINE BUILDS A 3x3 TRANSFORMATION MATRIX
C    BASED ON THE 2x2 DIRECTOR COSINE MATRIX
C
C**************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*6 IFLA1,IFLA2
C
      DIMENSION B(2,*), T(3,*)
C
      C=B(1,1)
      S=B(1,2)
C
      T(1,1)=C*C
      T(1,2)=S*S
      T(2,1)=S*S
      T(2,2)=C*C
      T(3,3)=C*C-S*S
C
C     CASE 1 :  'STRAIN' 'G_TO_L'
C
      IF(IFLA1.EQ.'STRAIN'.AND.IFLA2.EQ.'G_TO_L') THEN
C
      T(1,3)=C*S
      T(2,3)=-C*S
C
      T(3,1)=-2.0D00*C*S
      T(3,2)= 2.0D00*C*S
C
      ENDIF
C
C
C     CASE 2 :  'STRESS' 'G_TO_L'
C
      IF(IFLA1.EQ.'STRESS'.AND.IFLA2.EQ.'G_TO_L') THEN
C
      T(1,3)=2.0D00*C*S
      T(2,3)=-2.0D00*C*S
C
      T(3,1)=-C*S
      T(3,2)=C*S
C
      ENDIF
C
C
C     CASE 3 :  'STRESS' 'L_TO_G'
C
      IF(IFLA1.EQ.'STRESS'.AND.IFLA2.EQ.'L_TO_G') THEN
C
      T(1,3)=-2.0D00*C*S
      T(2,3)=2.0D00*C*S
C
      T(3,1)=C*S
      T(3,2)=-C*S
C
      ENDIF
C
C
C     CASE 4 :  'STRAIN' 'L_TO_G'
C
      IF(IFLA1.EQ.'STRAIN'.AND.IFLA2.EQ.'L_TO_G') THEN
C
      T(1,3)=-C*S
      T(2,3)=C*S
C
      T(3,1)=2.0D00*C*S
      T(3,2)=-2.0D00*C*S
C
      ENDIF
C
      RETURN
      END
