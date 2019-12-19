
      SUBROUTINE CRKLCL(DMATX,NCRAK,NSTRS,PROPS,STNCU,STSVL)
C*********************************************************************
C
C****THIS ROUTINE BUILDS THE D-MATRIX IN THE MATERIAL SYSTEM
C    ( FOR ELASTO-BRITTLE MATERIALS )
C
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DMATX(*), PROPS(*), STNCU(*), STSVL(*)
      DIMENSION E(3), G(3)
      DATA MCRAK/3/
C
      YOUNG=PROPS(2)
      SHEAR=YOUNG*0.5
      SHFAC=PROPS(34)
C
      DO 10 INDEX=1,3
      E(INDEX)=YOUNG
   10 G(INDEX)=SHEAR
C
      DO 20 ICRAK=1,NCRAK
      REFER=STSVL(ICRAK)
   20 IF(REFER.GT.0.0.AND.REFER.LT.YOUNG) E(ICRAK)=REFER
C
      G(1)=SHEAR*SHFAC
      G(2)=G(1)
      IF(NCRAK.GE.2) G(3)=G(1)
C
C***BUILD MATRIX
C
      DMATX(1)=E(1)
      DMATX(2)=E(2)
      DMATX(3)=G(1)
C
      IF(NSTRS.EQ.3) RETURN
C
      DMATX(4)=E(3)
C
      IF(NSTRS.EQ.4) RETURN
C
      DMATX(5)=G(2)
      DMATX(6)=G(3)
C
      RETURN
      END
