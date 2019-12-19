      SUBROUTINE STEPLA(AINTE,KGAUS,MINTE)
C***********************************************************************
C
C****THIS ROUTINE CALCULATES THE NUMBER OF STEPS FOR ELASTO-PLASTIC 
C    EQUILIBRIUM
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      MINTE=AINTE*10+5
      
      IF(MINTE.LT.500)RETURN
C
c     WRITE(7,900) KGAUS,MINTE
  900 FORMAT(3X,'SUB.FACTOR:THE CALCULATED MINTE AT GAUSS PT.',I5,
     . ' WAS EQUAL TO',I5,' THEREFORE IT HAS BEEN SETTED TO 500',/)
      MINTE=500.
C
      RETURN
      END
