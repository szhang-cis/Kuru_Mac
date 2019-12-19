      SUBROUTINE SOLCONT(RATAT)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE "COSOLI" COMMON PARAMETERS FOR THE
C     RESIDUAL FORCE SUBINCREMENTATION
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      IF(RATAT.GT.1000.0) THEN
       NSOL1=10
       TSOL1=0.7
       NSOL2=10
       TSOL2=0.5
       NSOL3=20
       TSOL3=0.5
      ENDIF
C
      IF(RATAT.LT.1000.0.AND.RATAT.GT.500.0) THEN
       NSOL1=5
       TSOL1=0.7
       NSOL2=10
       TSOL2=0.5
       NSOL3=20
       TSOL3=0.3
      ENDIF
C
      IF(RATAT.LT.500.0.AND.RATAT.GT.200.0) THEN
       NSOL1=4
       TSOL1=0.7
       NSOL2=8
       TSOL2=0.5
       NSOL3=16
       TSOL3=0.3
      ENDIF
C
      IF(RATAT.LT.200.0.AND.RATAT.GT.100.0) THEN
       NSOL1=3
       TSOL1=0.6
       NSOL2=6
       TSOL2=0.4
       NSOL3=12
       TSOL3=0.2
      ENDIF
C
      IF(RATAT.LT.100.0) THEN
       NSOL1=2
       TSOL1=0.7
       NSOL2=4
       TSOL2=0.5
       NSOL3=8
       TSOL3=0.3
      ENDIF
C
      RETURN
      END
