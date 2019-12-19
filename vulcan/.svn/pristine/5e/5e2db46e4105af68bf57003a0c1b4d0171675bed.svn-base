      SUBROUTINE SQTOTR(TRIAN,QUADR,CONST,NDIME)
C***********************************************************************
C
C**** THIS ROUTINE LOADS A SQUARE MATRIX INTO A TRIANGULAR ONE
C    
C             TRIAN = TRIAN + QUADR(NDIME,NDIME) * CONST
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION QUADR(NDIME,*), TRIAN(*)
C
      IKOVA=0
      DO 10 IDIME=1    ,NDIME
      DO 10 JDIME=IDIME,NDIME
      IKOVA=IKOVA+1
   10 TRIAN(IKOVA)=TRIAN(IKOVA)+QUADR(IDIME,JDIME)*CONST
C
      RETURN
      END
