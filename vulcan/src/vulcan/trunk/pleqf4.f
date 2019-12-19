      SUBROUTINE  PLEQF4(CAPAP,PEND,PREY,SIKMA,HARD0,NICUR,GEN)
C************************************************************
C***  FUNCION DE  C- CAPA --> EXPONENCIAL
C             DA  C- EPSI --> LINEAL
C***********************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PREY=0.1D-10
      PEND=0.0
      IF((CAPAP.LT.1.0D0).AND.(GEN.GT.0.0D0))THEN
       PECOM=-0.5D0*(SIKMA*SIKMA)/GEN
       IF(NICUR.EQ.0)THEN
        HARD0=PECOM
       ELSE
        IF(HARD0.LT.PECOM)THEN
c        PRINT 60, PECOM
         CALL RUNEND(' !!!ERROR!!! --  H  debe ser mayor que: ')
        ENDIF
       ENDIF
       PREY=DSQRT((SIKMA*SIKMA)+(2.0D0*HARD0*CAPAP*GEN))
       PEND=HARD0*GEN/PREY
      ENDIF
      RETURN
  60  FORMAT(' !!!ERROR!!! -- En la pendiente: H>',1E14.6)
      END
