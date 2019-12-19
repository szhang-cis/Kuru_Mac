      SUBROUTINE PLANFI(ANGFI,CAPAP,HANGF,PROPS)

C***********************************************************************
C
C**** ESTA RUTINA CALCULA LA VARIABLE INTERNA:
C               ANGFI
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION PROPS(*)
C
C**** OBTENCION DE LA FRICCION UNIAXIAL (caso de funcion explicita)
C
      CAPAL=0.8D0        !  0.0 < CAPAL < 1
      ANGFM=PROPS(14)*(1.0-CAPAP)
      SNFMA=DSIN(ANGFM)
      UVALU=2.0D0*SNFMA*(DSQRT(CAPAP*CAPAL)/(CAPAP+CAPAL))
      ANGFI=DASIN(UVALU)
      UPRIM=(2.0D0*SNFMA/(DSQRT(CAPAP*CAPAL)*(CAPAP+CAPAL)*
     .    (CAPAP+CAPAL)))*(0.5D0*CAPAL*(CAPAP+CAPAL)-(CAPAP*CAPAL))
      HANGF=UPRIM/(DSQRT(1.0D0-UVALU*UVALU))

      RETURN
      END
