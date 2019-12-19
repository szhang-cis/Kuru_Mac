
      SUBROUTINE DIRCOS(NDIME,NSTR1,SIGMA,SGPRI,BETAM)
C***********************************************************************
C
C****THIS ROUTINE COMPUTES THE EIGENVECTORS OF THE STRESS TENSOR
C
C    N.B.    SIGMA = ( Sxx , Syy , Sxy , Szz , Sxz , Syz )
C            SGPRI = ( S11 , S22 , S33 )         for NSTR1=6
C            SGPRI = ( S11 , S22 , Szz )         for NSTR1=4
C            SGPRI = ( S11 , S22 , 0.0 )         for NSTR1=3
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION BETAM(NDIME,*), SIGMA(*), SGPRI(*)
      DIMENSION DUMMY(3,3)
      DATA TWOPI,ERROR/6.283185307179586,1.0E-15/
C
      IF(NDIME.EQ.2) THEN
C
        STR11=SIGMA(1)
        STR22=SIGMA(2)
        STR12=SIGMA(3)
C
C       FIND ANGLE FOR MAX. & MIN. STRESSES
C
        DENTA=STR11-STR22
        IF(ABS(STR12).LE.ERROR) THEN
          ANGLE=0.0D00
          IF(DENTA.LT.0.0D00) ANGLE=TWOPI/4.0D00
        ELSE
          IF(ABS(DENTA).LE.ERROR) THEN
            ANGLE=TWOPI/8.0D00
            IF(STR12.LT.0.0D00) ANGLE=-TWOPI/8.0D00
          ELSE
            ANGLE=0.5D00*DATAN(2.0D00*STR12/DENTA)
          ENDIF
        ENDIF
C
        CA=COS(ANGLE)
        SA=SIN(ANGLE)
        SGPRI(1)=CA*CA*STR11+SA*SA*STR22+2.*SA*CA*STR12
        SGPRI(2)=SA*SA*STR11+CA*CA*STR22-2.*SA*CA*STR12
C
        IF(SGPRI(1).LT.SGPRI(2)) THEN
          ANGLE=ANGLE+TWOPI/4.0D00
          VALUE=SGPRI(1)
          SGPRI(1)=SGPRI(2)
          SGPRI(2)=VALUE
        ENDIF
C
C       BUILD TRANSFORMATION MATRIX
C
        CA=COS(ANGLE)
        SA=SIN(ANGLE)
C
        BETAM(1,1)= COS(ANGLE)
        BETAM(2,2)= COS(ANGLE)
        BETAM(2,1)=-SIN(ANGLE)
        BETAM(1,2)= SIN(ANGLE)
      ELSE
        CALL PDIREC(SIGMA,SGPRI,DUMMY)
        BETAM(1,1)=DUMMY(1,1)
        BETAM(1,2)=DUMMY(2,1)
        BETAM(1,3)=DUMMY(3,1)
        BETAM(2,1)=DUMMY(1,2)
        BETAM(2,2)=DUMMY(2,2)
        BETAM(2,3)=DUMMY(3,2)
        BETAM(3,1)=DUMMY(1,3)
        BETAM(3,2)=DUMMY(2,3)
        BETAM(3,3)=DUMMY(3,3)
      ENDIF
C
      RETURN
      END
