      SUBROUTINE RUNMENT(MESSAGET)
C***********************************************************************
C
C**** THIS ROUTINE WRITES MESSAGES
C
C***********************************************************************
      INCLUDE 'prob_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'
C
      CHARACTER*80 MESSAGET
C
      WRITE(LUPRIT,900) MESSAGET
      WRITE(LUREST,900) MESSAGET
C
      RETURN
C
  900 FORMAT(///,1X,A80)
      END
