      SUBROUTINE PRIDIR(NTYPE,STRSG,STRSP)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE PRINCIPAL DIRECTION OF THE 
C     STRESS TENSOR ( 2D & 3D ) 
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION STRSG(*), STRSP(*), VICTO(3,3)
C
      DATA TWOPI,ZM/6.283185307179586,1.0E-15/
C
      IF(NTYPE.EQ.4)GOTO 1
C
C**** 2-D CASE
C
        STR11=STRSG(1)
        STR22=STRSG(2)
        STR12=STRSG(3)
C
        DENTA=STR11-STR22
        IF(ABS(STR12).LE.ZM) THEN
          ANGLE=0.0D00
          IF(DENTA.LT.0.0D00) ANGLE=TWOPI/4.0D00
        ELSE
          IF(ABS(DENTA).LE.ZM) THEN
            ANGLE=TWOPI/8.0D00
            IF(STR12.LT.0.0D00) ANGLE=-TWOPI/8.0D00
          ELSE
            ANGLE=0.5D00*DATAN(2.0D00*STR12/DENTA)
          ENDIF
        ENDIF
C
        CA=DCOS(ANGLE)
        SA=DSIN(ANGLE)
        PSTR1=CA*CA*STR11+SA*SA*STR22+2.0D+00*SA*CA*STR12
        PSTR2=SA*SA*STR11+CA*CA*STR22-2.0D+00*SA*CA*STR12
        IF(PSTR1.GE.PSTR2) THEN
          STRSP(1)=PSTR1
          STRSP(2)=PSTR2
          STRSP(3)=ANGLE
        ELSE
          ANGLE=ANGLE+TWOPI/4.0D+00
          STRSP(1)=PSTR2
          STRSP(2)=PSTR1
          STRSP(3)=ANGLE
        ENDIF
C
      RETURN
C
C**** 3-D CASE
C
    1 CONTINUE
C     CALL PRINSI(STRSG,STRSP)
C     CALL PDIREC(STRSG,STRSP,VICTO)
      RETURN
      END
