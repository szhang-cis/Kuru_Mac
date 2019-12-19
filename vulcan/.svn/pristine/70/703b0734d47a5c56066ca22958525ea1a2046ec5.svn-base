      SUBROUTINE SHAPR3(S,T,Z,NDIME,NNODE,NGAUR,NQUTR,SHAPE)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES "REDUCED" SHAPE FUNCTIONS FOR
C     LINEAR AND QUADRATIC ISOPARAMETRIC 3-D ELEMENTS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION SHAPE(*)
C
      C1=1.0D0/0.577350269189626D0
C
C**** 4 NODED TETRAHEDRAL ELEMENT
C
      IF(NNODE.EQ.4) THEN
       NGAUR=1
       SHAPE(1)=1.0D0
       RETURN
      ENDIF
C
C**** 5 NODED TETRAHEDRAL ELEMENT
C
      IF(NNODE.EQ.5) THEN
c      NGAUR=1
c      SHAPE(1)=1.0D0

       NGAUR=4
       SHAPE(1)=0.25D0
       SHAPE(2)=0.25D0
       SHAPE(3)=0.25D0
       SHAPE(4)=0.25D0

c      NGAUR=4
c      A=(3.0D0+DSQRT(3.0D0))/2.0D0
c      B=DSQRT(3.0D0)/(2.0D0*A)
c      C=1.0D0/(2.0D0*A)
c      SHAPE(1)=B
c      SHAPE(2)=C
c      SHAPE(3)=C
c      SHAPE(4)=C

c      NGAUR=4
c      SHAPE(1)=1.0D0-S-T-Z
c      SHAPE(2)=S
c      SHAPE(3)=T
c      SHAPE(4)=Z

c      NGAUR=1
c      SHAPE(1)=256.0D0*(1.0D0-S-T-Z)*S*T*Z

c      NGAUR=4
c      A1=1.0D0-S-T-Z
c      A2=S
c      A3=T
c      A4=Z
C
c      SHAPE(   1) =(2.0D0*A1-1.0D0)*A1
c      SHAPE(   2) =(2.0D0*A2-1.0D0)*A2
c      SHAPE(   3) =(2.0D0*A3-1.0D0)*A3
c      SHAPE(   4) =(2.0D0*A4-1.0D0)*A4

c      NGAUR=3
c      S5S=256.0D0*(T*Z-2.0D0*S*T*Z-T*T*Z-T*Z*Z)
c      S5T=256.0D0*(S*Z-2.0D0*S*T*Z-S*S*Z-S*Z*Z)
c      S5Z=256.0D0*(S*T-2.0D0*S*T*Z-S*S*T-S*T*T)

c      SHAPE(1)=S5S/25.6
c      SHAPE(2)=S5T/25.6
c      SHAPE(3)=S5Z/25.6

c      NGAUR=4
c      C1=1.0D0/0.58541020D+00
c      S=S*C1
c      T=T*C1
c      Z=Z*C1
c      SHAPE(1)=1.0D0-T-S-Z
c      SHAPE(2)=S
c      SHAPE(3)=T
c      SHAPE(4)=Z

       RETURN
      ENDIF
C
C**** 10 NODED TETRAHEDRAL ELEMENT
C
      IF(NNODE.EQ.10) THEN

c      NGAUR=1
c      SHAPE(1)=1.0D0

       NGAUR=4
       SHAPE(1)=1.0D0-S-T-Z
       SHAPE(2)=S
       SHAPE(3)=T
       SHAPE(4)=Z

c      NGAUR=10
c      A1=1.0D0-S-T-Z
c      A2=S
c      A3=T
c      A4=Z
C
c      SHAPE(   1) =(2.0D0*A1-1.0D0)*A1
c      SHAPE(   2) =(2.0D0*A2-1.0D0)*A2
c      SHAPE(   3) =(2.0D0*A3-1.0D0)*A3
c      SHAPE(   4) =(2.0D0*A4-1.0D0)*A4
c      SHAPE(   5) = 4.0D0*A1*A2
c      SHAPE(   6) = 4.0D0*A2*A3
c      SHAPE(   7) = 4.0D0*A1*A3
c      SHAPE(   8) = 4.0D0*A1*A4
c      SHAPE(   9) = 4.0D0*A2*A4
c      SHAPE(  10) = 4.0D0*A3*A4

       RETURN
      ENDIF
C
C**** 8 NODED BRICK ELEMENT
C
C     Note:
C     The option NGAUR=1 makes the stiffness matrix not positive
C     definite using only one element.
C
      IF(NNODE.EQ.8) THEN
       NGAUR=1
       SHAPE(1)=1.0D0
       RETURN
      ENDIF
C
C**** 20 & 27 NODED BRICK ELEMENTS
C
      IF(NNODE.EQ.20.OR.NNODE.EQ.27) THEN
       S=S*C1
       T=T*C1
       Z=Z*C1
       SM = 0.5D0*(1.-S)
       TM = 0.5D0*(1.-T)
       ZM = 0.5D0*(1.-Z)
       SP = 0.5D0*(1.+S)
       TP = 0.5D0*(1.+T)
       ZP = 0.5D0*(1.+Z)
C
       SHAPE(1)=SM*TM*ZM
       SHAPE(2)=SP*TM*ZM
       SHAPE(3)=SP*TP*ZM
       SHAPE(4)=SM*TP*ZM
       SHAPE(5)=SM*TM*ZP
       SHAPE(6)=SP*TM*ZP 
       SHAPE(7)=SP*TP*ZP
       SHAPE(8)=SM*TP*ZP
       RETURN
      ENDIF
C
      END
