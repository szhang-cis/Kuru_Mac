      SUBROUTINE SQTOTA(TRIAN,QUADR,CONST,NDIME,KSYMM)
C***********************************************************************
C
C**** THIS ROUTINE LOADS:
C      
C     - A SQUARE MATRIX INTO A TRIANGULAR ONE FOR SYMMETRIC PROBLEMS
C       (KSYMM.EQ.1)
C    
C             TRIAN = TRIAN + QUADR(NDIME,NDIME) * CONST
C
C     - A SQUARE MATRIX INTO A SQUARE ONE FOR UNSYMMETRIC PROBLEMS
C       (KSYMM.EQ.0)
C
C             TRIAN = TRIAN + QUADR(NDIME,NDIME) * CONST
C    
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION QUADR(NDIME,*), TRIAN(*)
C
      IKOVA=0
      DO IDIME=1,NDIME
       INDEX=IDIME
       IF(KSYMM.EQ.0) INDEX=1
       DO JDIME=INDEX,NDIME
       IKOVA=IKOVA+1
        TRIAN(IKOVA)=TRIAN(IKOVA)+QUADR(IDIME,JDIME)*CONST
       ENDDO
      ENDDO
C
      RETURN
      END
