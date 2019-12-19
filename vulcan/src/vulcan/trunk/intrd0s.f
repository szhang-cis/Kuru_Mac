      SUBROUTINE INTRD0S(DTIMES,KSTEPS,NSTEPS,NSSTEPS)
C*********************************************************************
C
C**** THIS ROUTINE READS IN "INTERVAL_DATA" CARDS
C
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'inpo_oms.f'
      INCLUDE 'prob_oms.f'
C
      DTIMES=1.
      NSTEPS=1
      KSTEPS=1
      NSSTEPS=1
C
      NPRINS=0
      ITAPES=LUDATS
  100 CALL LISTENS('INTRD0S',NPRINS,ITAPES)
C
C**** 'INTER'VAL CARD IS FOUND
C
      IF(WORDSS(1).EQ.'INTER') THEN
       IF(PARAMS(1).GT.0.) NSTEPS=INT(PARAMS(1))
       IF(PARAMS(2).GT.0.) DTIMES=PARAMS(2)
       IF(PARAMS(3).GT.0.) NSSTEPS=INT(PARAMS(3))
       RETURN
      ENDIF
C
C**** 'STOP' CARD IS FOUND 
C                           
      IF(WORDSS(1).EQ.'STOP')
     . CALL RUNENDS('  * * *   END OF ANALISYS   * * *  ')
C
      GO TO 100
C
      END
