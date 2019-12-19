      SUBROUTINE  PLEQF3(CAPAP,PEND,PREY,SIKMA,SIKPI,CAPPI)
C************************************************************
C***  FUNCION DE  C- CAPA --> POLINOMICA CUADRATICA Y CUBICA
C             DA  C- EPSI --> POLINOMICA
C***********************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PREY=-0.1D-10
      PEND=0.0D0
      IF(CAPAP.LT.1.0D0)THEN
        IF(CAPAP.LE.CAPPI)THEN
C*         PARABOLA CUADRATICA (Endurecimiento)
            COEF=(-SIKPI+SIKMA)/(CAPPI*CAPPI)
            PREY=(COEF*(-CAPPI+CAPAP)*(-CAPPI+CAPAP))+SIKPI
            PEND=2.0*COEF*(-CAPPI+CAPAP)
        ELSE
C*         POLINOMIO DE TERCER GRADO (Ablandamiento)
            COEF1=CAPAP-CAPPI
            COEF2=1.0D0-CAPPI
            COEF3=COEF1/COEF2
            PREY=(1.0-(3.0*COEF3*COEF3)+(2.0*COEF3*COEF3*COEF3))*SIKPI
            PEND=((-6.0*COEF1/(COEF2*COEF2))+(6.0*(COEF1*COEF1)/
     .              (COEF2*COEF2*COEF2)))*(SIKPI)
         ENDIF
        ENDIF
      RETURN
      END
