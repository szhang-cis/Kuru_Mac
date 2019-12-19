      SUBROUTINE INTRD0T(DTIMET,KSTEPT,NSTEPT,NSSTEPT)
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
      COMMON/NOENDRUN/IRUNEN
C
      DTIMET=1.
      NSTEPT=1
      KSTEPT=1
      NSSTEPT=1
      IRUNEN=0
C
      NPRINT=0
      ITAPET=LUDATT
  100 CALL LISTENT('INTRD0T',NPRINT,ITAPET)
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
C     IT MEANS: THE THERMAL PROBLEM FINISHES
C               THE MECHANICAL PROBLEM CONTINUES
C
      IF(WORDST(1).EQ.'STOP') THEN
       CALL RUNENTT('  * * *   END OF ANALISYS   * * *  ')
       IRUNEN=1
       RETURN
      ENDIF
C
      GO TO 100
C
      END
