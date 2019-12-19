      SUBROUTINE PRISTN(NSTR1,STRAN,STPRI)
C***********************************************************************
C
C****THIS ROUTINE CALCULATES THE PRINCIPAL STRAIN VALUES OF
C    STRAN(NSTR1)  AND PUTS THE RESULTS IN STPRI(3)
C
C    N.B.    STRAN = ( Sxx , Syy , Sxy , Szz , Sxz , Syz )
C            STPRI = ( S11 , S22 , S33 )         for NSTR1=6
C            STPRI = ( S11 , S22 , Szz )         for NSTR1=4
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION STRAN(*), STPRI(*)
C
      IF(NSTR1.LT.6) THEN
        S11=STRAN(1)
        S22=STRAN(2)
        S12=STRAN(3)/2.0D00
        STPRI(1)=(S11+S22)/2.0+SQRT(0.25*(S11-S22)**2+S12**2)
        STPRI(2)=(S11+S22)/2.0-SQRT(0.25*(S11-S22)**2+S12**2)
        STPRI(3)=STRAN(4)
      ELSE
        CALL PRINST(STRAN,STPRI)
      ENDIF
C
      RETURN
      END   
