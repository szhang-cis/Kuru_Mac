      SUBROUTINE  PLEQF2(CAPAP,PEND,PREY,SIKMA,SIKPI,CAPPI)
C************************************************************
C***  FUNCION DE  C- CAPA --> EXPONENCIAL
C             DA  C- EPSI --> POLINOMICA
C***********************************************************
 
      IMPLICIT REAL*8 (A-H,O-Z)

      PREY=-0.1D-10
      PEND=0.0D0
      IF(CAPAP.LT.1.0D0)THEN
       RO=DSQRT(1.0D0-(SIKMA/SIKPI))
       ALFA=DLOG((1.0D0-(1.0D0-RO)*(1.0D0-RO))/
     .            ((3.0D0-RO)*(1.0D0+RO)*(CAPPI)))
       ALFA=DEXP(ALFA/(1.0D0-CAPPI))
       FI=((1.0D0-RO)*(1.0D0-RO))+((3.0D0-RO)*(1.0D0+RO)*CAPAP
     .     *((ALFA)**(1.0D0-CAPAP)))
       PREY=SIKPI*(2.0D0*DSQRT(FI)-FI)
       PEND=SIKPI*((1.0D0/(DSQRT(FI)))-1.0D0)
     .          *(3.0D0-RO)*(1.0D0+RO)*((ALFA)**(1.0D0-CAPAP))
     .          *(1.0D0-DLOG(ALFA)*CAPAP)
       ENDIF 
 
      RETURN
      END
