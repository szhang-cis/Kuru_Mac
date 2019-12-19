      SUBROUTINE JACOBSS(CARTD,DERIV,DETJM,ELCOD,GPCOD,IELEM,NDIME,
     .                   NNODE,SHAPE,XJACM,LURES,LUPRI)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE JACOBIAN MATRIX AND THE CARTESIAN
C     SHAPE FUNCTION DERIVATIVES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON/JACOBSSA/IERORS,KERORS
C
      DIMENSION CARTD(NDIME,*), DERIV(NDIME,*), ELCOD(NDIME,*),
     .          GPCOD(*),       SHAPE(*)      , XJACM(NDIME,*)
      DIMENSION XJACI(3,3)
C
C**** CALCULATE COORDINATES OF SAMPLING POINT
C
      DO 10 IDIME=1,NDIME
      GPCOD(IDIME)=0.0
      DO 10 INODE=1,NNODE
      GPCOD(IDIME)=GPCOD(IDIME)+
     .                   ELCOD(IDIME,INODE)*SHAPE(INODE)
   10 CONTINUE
C
C**** CREATE JACOBIAN MATRIX
C
      CALL PROMA2(XJACM,DERIV,ELCOD,NDIME,NNODE)
C
C**** CALCULATE THE DETERMINANT AND INVERS OF JACOBIAN MATRIX
C
      CALL INVMTX(XJACM,XJACI,DETJM,NDIME)
C
C**** CHECK IF ITS POSITIVENESS
C
      IF(DETJM.LE.0.0) THEN
        IERORT=IERORT+1
        KERORT=KERORT+1
        WRITE(LURES,600) IELEM,
     .      ((ELCOD(IDIME,INODE),IDIME=1,NDIME),INODE=1,NNODE)
        WRITE(LUPRI,601) IELEM
        RETURN
      ENDIF
C
C**** CALCULATE CARTESIAN DERIVATIVES
C
      CALL PROMA1(CARTD,XJACI,DERIV,NDIME,NNODE,NDIME)
C
      RETURN
  600 FORMAT(/,'ZERO OR NEGATIVE AREA IN ELEMENT NUMBER ',I5,
     .        /,10(3(F15.5,3X),/))
  601 FORMAT(/,'ZERO OR NEGATIVE AREA IN ELEMENT NUMBER ',I5)
      END
