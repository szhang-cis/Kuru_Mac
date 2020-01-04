      SUBROUTINE MIDPOI(COORD,IELEM,LNODS,NDIME,NNODE)
C***********************************************************************
C
C**** THIS ROUTINE GENERATES THE CENTRAL POINT FOR 9-NODED 2D
C     ISOPARAMETRIC ELEMENTS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION COORD(NDIME,*), LNODS(*)
      DIMENSION DERIV(2,8),     ELCOD(2,8), SHAPE(8)
C
      MNODE=8
      ETASP=0.0
      EXISP=0.0
C
C***EVALUATE THE COORDINATES OF THE ELEMENT NODALS POINT
C
      DO 10 INODE=1,MNODE
      KLOCA=LNODS(INODE)
      DO 10 IDIME=1,NDIME
      ELCOD(IDIME,INODE)=COORD(IDIME,KLOCA)
   10 CONTINUE
C
C***COMPUTE THE SHAPE FUNCTION
C
      CALL SHAFUN(DERIV,EXISP,ETASP,EZETA,NDIME,    8,    0,    0,SHAPE)
C
C***COMPUTE THE GLOBAL COORDINATES
C
      KLOCA=LNODS(9)
      DO 20 IDIME=1,NDIME
      COORD(IDIME,KLOCA)=0.0
      DO 20 INODE=1,NNODE
      COORD(IDIME,KLOCA)=COORD(IDIME,KLOCA)+SHAPE(INODE)*
     .                   ELCOD(IDIME,INODE)
   20 CONTINUE
C
      RETURN
      END