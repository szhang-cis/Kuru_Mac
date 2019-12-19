      SUBROUTINE INTRD0T0(DTIMET,KSTEPT,NSTEPT,NSSTEPT)
C*********************************************************************
C
C**** THIS ROUTINE READS IN "INTERVAL_DATA" CARDS
C
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'inpo_omt.f'
      INCLUDE 'prob_omt.f'
C
      DTIMET=1.
      NSTEPT=1
      KSTEPT=1
      NSSTEPT=1
C
      NPRINT=0
      ITAPET=LUDATT
  100 CALL LISTENT('INTRD0T0',NPRINT,ITAPET)
C
C**** 'INTER'VAL CARD IS FOUND
C
      IF(WORDST(1).EQ.'INTER') THEN
       IF(PARAMT(1).GT.0.) NSTEPT=INT(PARAMT(1))
       IF(PARAMT(2).GT.0.) DTIMET=PARAMT(2)
       IF(PARAMT(3).GT.0.) NSSTEPT=INT(PARAMT(3))
       RETURN
      ENDIF
C
C**** 'STOP' CARD IS FOUND 
C                           
      IF(WORDST(1).EQ.'STOP')
     . CALL RUNENDT('  * * *   END OF ANALISYS   * * *  ')
C
      GO TO 100
C
      END
