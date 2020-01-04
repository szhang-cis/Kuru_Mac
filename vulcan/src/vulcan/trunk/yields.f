      SUBROUTINE YIELDS(AVECT,DEFPI,EHIST,HARDS,NCRIT,NSTR1,NTYPE,
     .                  PROPS,SGTOT,VSTOT)
C***********************************************************************
C
C****THIS ROUTINE DETERMINES THE GRADIENT TO THE POTENTIAL FUNCTION
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AVECT(*), DEVIA(6), EHIST(*), PROPS(*), SGTOT(*),
     .          VECA3(6)
C
C***CALCULATE STRESS INVARIANTS
C
      CALL INVARS(SGTOT,DEVIA,NTYPE,PMEAN,QINVA,THETA,VARJ2,
     .            VARJ3)
C
C***DETERMINE THE GRADIENT COEFFICIENTS
C
      ITEST=1
      CALL GRACOE(CONS1,CONS2,CONS3,DEFPI,EHIST,FAC11,FAC22,
     .            FAC33,FAC12,FAC13,FAC23,HARDS,IFLAG,ITEST,
     .            NCRIT,PMEAN,PROPS,QINVA,THETA,VSTOT)
C
C***COMPUTE GRADIENT VECTOR
C
      CALL GRAVEC(AVECT,CONS1,CONS2,CONS3,DEVIA,NSTR1,VARJ2,
     .            VECA3)
C
      RETURN
      END