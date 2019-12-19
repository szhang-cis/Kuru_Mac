      SUBROUTINE  PLEQF5(CAPAP,PREY,PEND,SIKPI,SIKMA,ICODI)
C************************************************************
C***  CURVA PARA COMPACTACION c-rho  (PARABOLA): codigo -1-
C***  CURVA PARA COMPACTACION c-capa (RECTA)   : codigo -2-
C***********************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      IF(ICODI.EQ.1)THEN
        PREY=SIKMA+(SIKPI-SIKMA)*(CAPAP)**2.0
        PEND=2.0*(SIKPI-SIKMA)*CAPAP
      ELSE
        PREY=SIKMA+(SIKPI-SIKMA)*CAPAP
        PEND=SIKPI-SIKMA
      ENDIF
      RETURN
      END
