      SUBROUTINE  PLEQF1(CAPAP,PEND,PREY,SIKMA,NICUR,HARD0,GEN)
C************************************************************
C***  FUNCION DE  C- CAPA --> LINEAL
C             DA  C- EPSI --> EXPONENCIAL
C***********************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PREY=0.1D-10
      PEND=0.0
      IF(CAPAP.LT.1.0D0)THEN
       IF(NICUR.EQ.0)THEN
        PREY=SIKMA*(1.0-CAPAP)
        PEND=-SIKMA
       ELSE
         IF(HARD0.LT.(-SIKMA))THEN
c         PRINT 60, (-SIKMA)
          CALL RUNEND(' !!!ERROR!!! --  H  debe ser mayor que: ')
         ENDIF
        PREY=SIKMA+HARD0*CAPAP*GEN
        PEND=HARD0
       ENDIF
      ENDIF
      RETURN
  60  FORMAT(' !!!ERROR!!! -- En la pendiente: H>',1E14.6)
      END
