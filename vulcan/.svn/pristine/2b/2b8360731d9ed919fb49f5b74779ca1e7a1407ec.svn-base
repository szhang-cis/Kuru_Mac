      SUBROUTINE SHAPR2(S,T,NDIME,NNODE,NGAUR,NQUTR,SHAPE)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES "REDUCED" SHAPE FUNCTIONS FOR
C     LINEAR AND QUADRATIC ISOPARAMETRIC 2-D ELEMENTS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION SHAPE(*)
C
C**** 3 NODED ELEMENT
C
      IF(NNODE.EQ.3) THEN
       NGAUR=1
       SHAPE(1)=1.0D0
       RETURN
      ENDIF
C
C**** 4 NODED ELEMENT
C
      IF(NNODE.EQ.4) THEN
c      IF(NQUTR.EQ.1) THEN
        NGAUR=1
        SHAPE(1)=1.0D0
c      ENDIF

c      NGAUR=1
c      SHAPE(1)=27.0D0*(1.0D0-S-T)*S*T
C
c      NGAUR=3
c      SHAPE(1)=1.0D0-S-T
c      SHAPE(2)=S
c      SHAPE(3)=T
C
c      NGAUR=3
c      A1=1.0D0-S-T
c      A2=S
c      A3=T
C
c      SHAPE(1)=(2.0D0*A1-1.0D0)*A1
c      SHAPE(2)=(2.0D0*A2-1.0D0)*A2
c      SHAPE(3)=(2.0D0*A3-1.0D0)*A3
C
c      SHAPE(1)=4.0D0*A1*A2
c      SHAPE(2)=4.0D0*A2*A3
c      SHAPE(3)=4.0D0*A1*A3
C
       RETURN
      ENDIF
C
C**** 5 NODED ELEMENT
C
      IF(NNODE.EQ.5) THEN
       NGAUR=1
       SHAPE(1)=1.0D0
       RETURN
      ENDIF
C
C**** 6 NODED ELEMENT
C
      IF(NNODE.EQ.6) THEN
C
       IOPTION=1                          ! better as input
       IF(IOPTION.EQ.1) THEN
        NGAUR=1
        SHAPE(1)=1.0D0
       ENDIF
C
       IF(IOPTION.EQ.2) THEN
        NGAUR=3
        SHAPE(1)=1.0D0-S-T
        SHAPE(2)=S
        SHAPE(3)=T
       ENDIF
C
       RETURN
      ENDIF
C
C**** 8 & 9 NODED ELEMENTS
C
      IF(NNODE.EQ.8.OR.NNODE.EQ.9) THEN
C
       IOPTION=3                          ! better as input
       IF(IOPTION.EQ.1) THEN
        NGAUR=1
        SHAPE(1)=1.0D0
       ENDIF
C
       IF(IOPTION.EQ.2) THEN
        NGAUR=3
        SHAPE(1)=1.0D0
        SHAPE(2)=S
        SHAPE(3)=T
       ENDIF
C
       IF(IOPTION.EQ.3) THEN
        NGAUR=4
        C1=1.0D0/0.577350269189626D0
        S=S*C1
        T=T*C1
        ST=S*T
        SHAPE(1)=(1.0D0-T-S+ST)*0.25D0
        SHAPE(2)=(1.0D0-T+S-ST)*0.25D0
        SHAPE(3)=(1.0D0+T+S+ST)*0.25D0
        SHAPE(4)=(1.0D0+T-S-ST)*0.25D0
       ENDIF
C
       RETURN
      ENDIF
C
      END
