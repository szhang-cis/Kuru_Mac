      SUBROUTINE PRIVAL(NSTR1,SIGMA,SGPRI)
C***********************************************************************
C
C****THIS ROUTINE CALCULATES THE PRINCIPAL STRESS VALUES OF
C    SIGMA(NSTR1)  AND PUTS THE RESULTS IN SGPRI(3)
C
C    N.B.    SIGMA = ( Sxx , Syy , Sxy , Szz , Sxz , Syz )
C            SGPRI = ( S11 , S22 , S33 )         for NSTR1=6
C            SGPRI = ( S11 , S22 , Szz )         for NSTR1=4
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SIGMA(*), SGPRI(*)
C
      IF(NSTR1.LT.6) THEN
        S11=SIGMA(1)
        S22=SIGMA(2)
        S12=SIGMA(3)
        SGPRI(1)=(S11+S22)/2.0+SQRT(0.25*(S11-S22)**2+S12**2)
        SGPRI(2)=(S11+S22)/2.0-SQRT(0.25*(S11-S22)**2+S12**2)
        SGPRI(3)=SIGMA(4)
      ELSE
        CALL PRINSI(SIGMA,SGPRI)
      ENDIF
C
      RETURN
      END   
