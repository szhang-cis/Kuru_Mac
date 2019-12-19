      SUBROUTINE PLANSI(ANGSI,CAPAP,HANGS,PROPS,ANGFI,
     .                  DFICA)
C***********************************************************************
C
C**** ESTA RUTINA CALCULA LA VARIABLE INTERNA:
C               ANGSI
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION PROPS(*)
C
C**** OBTENCION DE LA DILATANCIA UNIAXIAL (caso de funcion explicita)
C
      ANGFM=PROPS(14)*(1.0-CAPAP)
      SNFMA=DSIN(ANGFM)
      ANGSM=PROPS(23)
      SNSMA=DSIN(ANGSM)
      ANFCV=DASIN((SNFMA-SNSMA)/(1.0D0-SNFMA*SNSMA))
      ANGSI=0.0D0
      IF(ANGFI.GT.ANFCV)ANGSI=ANGFI-ANFCV
      DSIFI=1.0D0
      HANGS=DSIFI*DFICA

      RETURN
      END
