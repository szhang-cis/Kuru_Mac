      SUBROUTINE EXPONE(CONSS,FTULT,REFER,TSTRA,TSTRE,YOUNG)
      IMPLICIT REAL*8 (A-H,O-Z)
C**********************************************************************
C
C****THIS ROUTINE FOLLOWS THE LOAD-UNLOAD-RELOAD CURVE
C
C**********************************************************************
C  
      DEFO1=FTULT/YOUNG
      EXPUP=-(TSTRA-DEFO1)/CONSS
      IF(EXPUP.GT.0.5E+02) EXPUP=50.
C
      VALU1=TSTRA*REFER                   !  UNLOAD-RELOAD VALUE
      VALU2=FTULT*EXP(EXPUP)              !  LOAD VALUE
C
      TSTRE=MIN(FTULT,VALU1,VALU2)
      IF(TSTRE.LT.0.0) TSTRE=0.0
C
      RETURN
      END 
