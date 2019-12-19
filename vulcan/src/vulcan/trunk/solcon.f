      SUBROUTINE SOLCON(RATAT)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE "COSOLI" COMMON PARAMETERS FOR THE
C     RESIDUAL FORCE SUBINCREMENTATION
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      IF(RATAT.GT.1000.0) THEN
       NSOL1M=10
       TSOL1M=0.7
       NSOL2M=10
       TSOL2M=0.5
       NSOL3M=20
       TSOL3M=0.5
      ENDIF
C
      IF(RATAT.LT.1000.0.AND.RATAT.GT.500.0) THEN
       NSOL1M=5
       TSOL1M=0.7
       NSOL2M=10
       TSOL2M=0.5
       NSOL3M=20
       TSOL3M=0.3
      ENDIF
C
      IF(RATAT.LT.500.0.AND.RATAT.GT.200.0) THEN
       NSOL1M=4
       TSOL1M=0.7
       NSOL2M=8
       TSOL2M=0.5
       NSOL3M=16
       TSOL3M=0.3
      ENDIF
C
      IF(RATAT.LT.200.0.AND.RATAT.GT.100.0) THEN
       NSOL1M=3
       TSOL1M=0.6
       NSOL2M=6
       TSOL2M=0.4
       NSOL3M=12
       TSOL3M=0.2
      ENDIF
C
      IF(RATAT.LT.100.0) THEN
       NSOL1M=2
       TSOL1M=0.7
       NSOL2M=4
       TSOL2M=0.5
       NSOL3M=8
       TSOL3M=0.3
      ENDIF
C
      RETURN
      END
