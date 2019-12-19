      SUBROUTINE RECACOT(SOUR1T,SOUR2T,DTEMPT)
C***********************************************************************
C
C**** THIS ROUTINE CONTROLS THE TEMPERATURE DERIVATIVE OF THE
C     PHASE-CHANGE FUNCTION (RECALESCENSE EFFECT), ONLY FOR 
C     MICROSTRUCTURAL MODELS WITH RATE PHASE-CHANGE FORMULATIONS
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** CONTROLS df_pc/dT
C
      SPEAPT=SOUR2T-SOUR1T                           ! SOUR1T always 0.0
      IF(SPEAPT.LT.0.0D0) SOUR2T=0.0D0
C
      RETURN
      END
