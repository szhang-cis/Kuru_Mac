      SUBROUTINE RUNMENS(MESSAGES)
C***********************************************************************
C
C**** THIS ROUTINE WRITES MESSAGES
C
C***********************************************************************
      INCLUDE 'prob_oms.f'
      INCLUDE 'auxl_oms.f'
      INCLUDE 'inpo_oms.f'
C
      CHARACTER*35 MESSAGES
C
      IF(IEVFI.EQ.1) WRITE(LUPRIS,900) MESSAGES
      WRITE(LURESS,900) MESSAGES
C
      RETURN
C
  900 FORMAT(1H1,///,5X,A35)
      END
