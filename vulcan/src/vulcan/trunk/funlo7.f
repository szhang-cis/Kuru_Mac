      SUBROUTINE FUNLO7(FACTI,HTLOD,TTIME)
C***********************************************************************
C
C**** THIS ROUTINE MODELS A SAW_TOOTH SIGNAL (ONLY POSITIVE VALUES)
C
C         .
C       .   .
C     .       . . . . .
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION HTLOD(*)
C
      FACTI=0.0
      IF(TTIME.LE.HTLOD(2)) RETURN       ! NOT STARTED YET
      IF(TTIME.GT.HTLOD(3)) RETURN       ! ALREADY FINISHED
C
      PERI0=HTLOD(4)
      PERI1=PERI0*0.25
      PERI2=PERI0*0.75
      PERI3=PERI0*0.50
C
      RTIME=TTIME-HTLOD(2)
      KCICL=INT(RTIME/HTLOD(4))
      TAUTI=DABS(RTIME-KCICL*HTLOD(4))
C
      IF(TAUTI.LE.PERI1) 
     .  FACTI=HTLOD(5)*TAUTI/PERI1
C
      IF(TAUTI.GT.PERI1.AND.TAUTI.LT.PERI3)
     .  FACTI=HTLOD(5)*(1.-2.*(TAUTI-PERI1)/PERI3)
C
      IF(TAUTI.GE.PERI3) 
     .  FACTI=0.0
C      
      RETURN
      END
