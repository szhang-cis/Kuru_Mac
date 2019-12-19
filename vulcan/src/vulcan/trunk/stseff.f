      SUBROUTINE STSEFF(CONSS,FTULT,STNCU,STSVL,TSTRA,TSTRE,YOUNG)
C********************************************************************
C
C****THIS ROUTINE CALCULATES THE STRESS ACROSS THE CRACK
C
C   Input parameters:
C
C      CONSS      Softening parameter
C      FTULT      Tensile strength
C      TSTRA      Current total strain normal to crack
C      YOUNG      Young's modulus
C
C   I/O parameters:
C
C      STNCU      Current strain normal to crack
C      STSVL      Secant stiffness modulus for crack
C      TSTRE      Current stress normal to crack
C
C********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      IF(TSTRA.LT.0.0) RETURN             ! THE CRACK IS CLOSED
C
C***FIND CURRENT STRESS VALUE
C
      REFER=YOUNG
      IF(STNCU.EQ.0.0)   STSVL=YOUNG      ! JUST OPENED
      IF(STSVL.LT.YOUNG) REFER=STSVL
C
      CALL EXPONE(CONSS,FTULT,REFER,TSTRA,TSTRE,YOUNG)
C
C***UPDATE STIFFNESS VALUE
C 
      IF(TSTRA.NE.0.0D00) REFER=TSTRE/TSTRA
      IF(REFER.LT.STSVL)  STSVL=REFER
C
      STNCU=TSTRA
C
      RETURN
      END
