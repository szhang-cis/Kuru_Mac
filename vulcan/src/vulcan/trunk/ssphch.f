      SUBROUTINE SSPHCH(SSPCF,DSSPC,NELASF,NELASA,DELSS,TEMPG,DTEMG,
     .                  FRANU,TIMEN,PROPS)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE SOLID-SOLID PHASE-CHANGE FUNCTION
C     
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION PROPS(*)
C
c      IF(TEMPG.GE.TEMSS) THEN
c       SSPCF=1.0D+00
c       DSSPC=0.0D+00
c       NELASF=0
c       NELASA=1
C
c       FRANU=0.0
c       TIMEN=0.0
C
c       RETURN
c      ENDIF
C
c      IF(TEMPG.LT.TEMSS) THEN
C
C**** COMPUTES NUCLEATION FRACTION
C
c       DFRAN=(TEMSS-TEMPG)**2/(AKCOO*TEMPG)
c       IF(DFRAN.LT.0.0) DFRAN=0.0
c       FRANU=FRANU+DFRAN
c       TIMEN=TIMEN+DTIME
C
C**** COMPUTES SOLID-SOLID PHASE-CHANGE FUNCTION
C
c       FFRAN=0.0
c       IF(FRANU.GT.1.0) FFRAN=1.0
c       SSPCF=1.0-(1.0-EXP(-AKC11*EXP(-AKC22/TEMPG)*TIMEN**ANC11))*FFRAN
C
C**** COMPUTES SOLID-SOLID PHASE-CHANGE FUNCTION RATE
C
c       G1=1.0
c       G2=1.0
c       DSSPC=G1*G2*FFRAN
C
c       NELASF=0
c       NELASA=1
c       TOLNEA=1.0D-04
c       IF(SSPCF.LT.1.0)    NELASF=1
c       IF(SSPCF.LE.TOLNEA) NELASA=0
C
c       RETURN
c      ENDIF


      return
C
      END
