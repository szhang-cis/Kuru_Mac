      SUBROUTINE ROTSYS(BETAM,SIGMA,TRANS,NDIME,NSTRE,NSTRS)
C***********************************************************************
C
C**** THIS ROUTINE ROTATES THE LOCAL MATERIAL AXES
C
C     N.B.     BETA0 :   GLOBAL   --> LOCAL 1
C              BETA2 :   LOCAL 1  --> LOCAL 2
C              BETAM :   GLOBAL   --> LOCAL 2
C
C              (BETAM) = (BETA2).(BETA0)
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION BETAM(NDIME,*), SIGMA(*), TRANS(6,*)
      DIMENSION BETA0(3,3), BETA2(3,3)
C
      DATA TWOPI/6.283185307179586/
C
      STR11=SIGMA(2)
      STR22=SIGMA(4)
      STR12=SIGMA(6)
C
C***FIND ANGLE FOR MAX. & MIN. STRESSES
C
      ERROR=(ABS(STR11)+ABS(STR22))*1.0E-10
      DENTA=STR11-STR22
      IF(ABS(STR12).LE.ERROR) THEN
        ANGLE=0.0D+00
        IF(DENTA.LT.0.0) ANGLE=TWOPI/4.
      ELSE
        IF(ABS(DENTA).LE.ERROR) THEN
          ANGLE=TWOPI/8.
          IF(STR12.LT.0.0) ANGLE=-TWOPI/8.
        ELSE
          ANGLE=0.5*DATAN(2.*STR12/DENTA)
        ENDIF
      ENDIF
C
C***FIND MAX. & MIN. VALUE OF STRESSES
C
      CA=COS(ANGLE)
      SA=SIN(ANGLE)
      PSTR1=CA*CA*STR11+SA*SA*STR22+2.*SA*CA*STR12
      PSTR2=SA*SA*STR11+CA*CA*STR22-2.*SA*CA*STR12
      IF(PSTR1.LT.PSTR2) ANGLE=ANGLE+TWOPI/4.
C
      IF(ANGLE.EQ.0.0) RETURN
C
C***ROTATE MATERIAL AXES
C
C
C...store (beta0) : global  --> local 1
C
      DO 1 I=1,3
      DO 1 J=1,3
      BETA2(I,J)=0.0
      BETA0(I,J)=BETAM(I,J)
   1  BETAM(I,J)=0.0
C
C...build  (beta2) : local 1 --> local 2
C
      BETA2(1,1)= 1.
      BETA2(2,2)= COS(ANGLE)
      BETA2(3,3)= COS(ANGLE)
      BETA2(2,3)= SIN(ANGLE)
      BETA2(3,2)=-SIN(ANGLE)
C
C...build  (betam) : global  --> local 2
C
      DO 2 I=1,3
      DO 2 J=1,3
      DO 2 K=1,3
   2  BETAM(I,J)=BETAM(I,J)+BETA2(I,K)*BETA0(K,J)
C
C***TRANSFORM STRESSES TO ROTATED SYSTEM
C
      CALL TRASIG(SIGMA,BETA2,NSTRE,NDIME,'G_TO_L')
C
      RETURN
      END
