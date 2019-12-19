
      SUBROUTINE NOPOIS(NCRAK,NSTRS,PROPS,SIGMA,STRAN)
C********************************************************************
C
C****THIS ROUTINE RECOMPUTES PSEUDO-ELASTIC STRESSES 
C      NEGLECTING POISSON'S RATIO AND 
C      APPLYING SHEAR RETENTION FACTOR
C
C********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PROPS(*), SIGMA(*), STRAN(*)
C
      YOUNG=PROPS(2)
      POISS=PROPS(3)
      SHFAC=PROPS(34)
C
      SHEAR=SHFAC*YOUNG*0.5D00
C
      SIGMA(1)=STRAN(1)*YOUNG
      SIGMA(2)=STRAN(2)*YOUNG
      SIGMA(3)=STRAN(3)*SHEAR
C
      IF(NSTRS.EQ.6) THEN
        SIGMA(4)=STRAN(4)*YOUNG
        SIGMA(5)=STRAN(5)*SHEAR
        SIGMA(6)=STRAN(6)*SHEAR
      ELSE IF(NSTRS.GT.3) THEN
        SIGMA(4)=STRAN(4)*YOUNG
      ENDIF
C
      RETURN
      END
