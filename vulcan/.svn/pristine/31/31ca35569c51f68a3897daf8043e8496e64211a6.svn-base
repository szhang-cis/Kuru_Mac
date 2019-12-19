      SUBROUTINE NODES1(COORD,LNODS,NDIME,IELEM,NNODE)
C***********************************************************************
C
C**** THIS ROUTINE INTERPOLATES THE MID-SIDE NODE FOR STRAIGHT 1D 
C     ELEMENTS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION COORD(NDIME,*), LNODS(*)
C
      IF(NNODE.NE.3) RETURN
C
      LNODE=LNODS(3)
      TOTAL=DABS(COORD(1,LNODE))
      IF(TOTAL.EQ.0.0)
     .  COORD(1,LNODE)=(COORD(1,1)+COORD(1,2))*0.5
C
      RETURN
      END
