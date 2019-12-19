      SUBROUTINE DARECO(DBASE,DAMAG,TAUAM,TAUAP)
C***********************************************************************
C
C**** THIS SUBROUTINE RECOVERS INTERNAL VARIABLES & THEIR CONJUGATES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'prob_om.f'
C
      DIMENSION DBASE(*)
C
      DAMAG=0.0D0
      IF(IDAMG.EQ.0) RETURN
C
      IF(IDAMG.LE.20) THEN   ! damage models governed by pl. or viscopl.
c      IF(IDAMG.EQ.1) DAMAG=DBASE(1)
c      IF(IDAMG.EQ.4) DAMAG=DBASE(1)
c      IF(IDAMG.EQ.5) DAMAG=DBASE(1)
c      IF(IDAMG.EQ.6) DAMAG=DBASE(1)
       DAMAG=DBASE(1)
      ELSE                   ! damage models governed by a damage crit.
       IF(IDAMG.EQ.21) THEN  ! concrete damage model
        DAMAG=0.0D0
        TAUAM=DBASE(1)
        TAUAP=DBASE(2)
       ENDIF
      ENDIF
C
      RETURN
      END
