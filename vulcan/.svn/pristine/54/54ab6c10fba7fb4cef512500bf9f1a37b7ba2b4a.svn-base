      SUBROUTINE RULGAU(NGAUS,POSGP,WEIGP)
C***********************************************************************
C
C****THIS ROUTINE SETS UP THE GAUSS-LEGENDRE INTEGRATION CONSTANTS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION POSGP(*), WEIGP(*)
      IF(NGAUS-2) 1,2,4
C
    1 POSGP(1)=0.0
      WEIGP(1)=2.0
C
      RETURN
C
    2 POSGP(1)=-0.577350269189626
      WEIGP(1)=1.0
      GOTO 6
C
    4 POSGP(1)=-0.774596669241483
      POSGP(2)=0.0
      WEIGP(1)=0.555555555555556
      WEIGP(2)=0.888888888888889
C
    6 KGAUS=NGAUS/2
      DO 8 IGASH=1,KGAUS
      JGASH=NGAUS+1-IGASH
      POSGP(JGASH)=-POSGP(IGASH)
      WEIGP(JGASH)=WEIGP(IGASH)
    8 CONTINUE
C
      RETURN
      END
