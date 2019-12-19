      SUBROUTINE BRIBA1(EHIST,KGAUS,NCRIT,NTYPE,PROPS,SGTOT)
C***********************************************************************
C
C****THIS ROUTINE DETERMINE THE ELLIPSE PASSING THROUGH A GIVEN
C    STRESS POINT
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DEVIA(6), EHIST(*), PROPS(*), SGTOT(*)
      DCOEF=PROPS(31)
      RATIO=PROPS(34)
      FACTO=PROPS(35)
C
      CALL INVARS(SGTOT,DEVIA,NTYPE,PMEAN,QINVA,THETA,VARJ2,
     .            VARJ3)
      IF(NCRIT.EQ.4)GOTO 10
C
C***CALCULATE RATIO
C
      CALL VARTHE(CAPTH,CINTE,FINCL,IFLAG,PROPS,THETA)
      RATIO=FINCL*FACTO 
C
C***CALCULATE MAJOR ELLIPSE SEMI-AXIS
C
   10 CONTINUE
      DCOPM=DCOEF-PMEAN
      ACOEF=FACTO*FACTO-1.
      IF(DABS(ACOEF).LT.1.E-03)ACOEF=0.
      BCOEF=FACTO*DCOPM
      CCOEF=DCOPM*DCOPM+QINVA*QINVA/(RATIO*RATIO)
C
      IF(ACOEF.NE.0.)HPARA=(BCOEF+DSQRT(BCOEF*BCOEF-ACOEF*CCOEF))/ACOEF
      IF(ACOEF.EQ.0.)HPARA=CCOEF/(2*BCOEF)
      IF(HPARA.LE.0.)GOTO 20
      BETAS=DCOEF-FACTO*HPARA
      EQUMN=BETAS-HPARA
C
      EHIST(3)=HPARA
      EHIST(4)=BETAS
      EHIST(5)=EQUMN
C
      RETURN
C
   20 CONTINUE
c     WRITE(7,900)HPARA,KGAUS
  900 FORMAT(10X,'ELIPSE MAJOR SEMI-AXIS LT. 0. I.E. HPARA =',E15.6,
     .           ' AT GAUSS POINT NO. ',I5)
C
      CALL RUNEND('ERROR IN BRIBA1                    ')
      END
