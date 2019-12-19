      SUBROUTINE INTRD0(DTIME,KSTEP,NSTEP,NSSTEP)
C*********************************************************************
C
C**** THIS ROUTINE READS IN "INTERVAL_DATA" CARDS
C
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'inpo_om.f'
      INCLUDE 'prob_om.f'
C
      DTIME=1.
      NSTEP=1
      KSTEP=1
      NSSTEP=1
C
      NPRIN=0
      ITAPE=LUDAT
  100 CALL LISTEN('INTRD0',NPRIN,ITAPE)
C
C**** 'INTER'VAL CARD IS FOUND
C
      IF(WORDS(1).EQ.'INTER') THEN
       IF(PARAM(1).GT.0.) NSTEP=INT(PARAM(1))
       IF(PARAM(2).GT.0.) DTIME=PARAM(2)
       IF(PARAM(3).GT.0.) NSSTEP=INT(PARAM(3))
       RETURN
      ENDIF
C
C**** 'STOP' CARD IS FOUND
C
      IF(WORDS(1).EQ.'STOP')
     . CALL RUNEND('  * * *   END OF ANALISYS   * * *  ')
C
      GO TO 100
C
      END
