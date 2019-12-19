      SUBROUTINE SHAPE3(DERIV,S,T,Z,NDIME,NNODE,NQUTR,SHAPE)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES SHAPE FUNCTIONS AND THEIR DERIVATIVES FOR
C     LINEAR AND QUADRATIC ISOPARAMETRIC 3-D ELEMENTS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DERIV(NDIME,*),SHAPE(*)
C
C**** 4 NODED TETRAHEDRAL ELEMENT
C
      IF(NNODE.EQ.4) THEN
       IF(NQUTR.EQ.1) THEN
C
        SHAPE(   1) = 1.0D0-S-T-Z
        DERIV(1, 1) =-1.0D0
        DERIV(2, 1) =-1.0D0
        DERIV(3, 1) =-1.0D0
C
        SHAPE(   2) = S
        DERIV(1, 2) = 1.0D0
        DERIV(2, 2) = 0.0D0
        DERIV(3, 2) = 0.0D0
C
        SHAPE(   3) = T
        DERIV(1, 3) = 0.0D0
        DERIV(2, 3) = 1.0D0
        DERIV(3, 3) = 0.0D0
C
        SHAPE(   4) = Z
        DERIV(1, 4) = 0.0D0
        DERIV(2, 4) = 0.0D0
        DERIV(3, 4) = 1.0D0
C
       ENDIF
C
       IF(NQUTR.EQ.2) THEN
C
        SM = 0.5D0*(1.0D0-S)
        TM = 0.5D0*(1.0D0-T)
        ZM = 0.5D0*(1.0D0-Z)
        SP = 0.5D0*(1.0D0+S)
        TP = 0.5D0*(1.0D0+T)
        ZP = 0.5D0*(1.0D0+Z)
C
        SHAPE(   1) = +     SM*   ZM
        DERIV(1, 1) = -.5D0*      ZM
        DERIV(2, 1) = +.0D0
        DERIV(3, 1) = -.5D0*SM
C
        SHAPE(   2) = +     SP*TM*ZM
        DERIV(1, 2) = +.5D0   *TM*ZM
        DERIV(2, 2) = -.5D0*SP   *ZM
        DERIV(3, 2) = -.5D0*SP*TM
C
        SHAPE(   3) = +        TP*ZM
        DERIV(1, 3) = +.0D0
        DERIV(2, 3) = +.5D0*      ZM
        DERIV(3, 3) = -.5D0*   TP
C
        SHAPE(   4) =             ZP
        DERIV(1, 4) = +.0D0
        DERIV(2, 4) = +.0D0
        DERIV(3, 4) = +.5D0
C
       ENDIF
C
       RETURN
      ENDIF
C
C**** 5 NODED TETRAHEDRAL ELEMENT
C
      IF(NNODE.EQ.5) THEN
       S5=256.0D0*(1.0D0-S-T-Z)*S*T*Z
       S5S=256.0D0*(T*Z-2.0D0*S*T*Z-T*T*Z-T*Z*Z)
       S5T=256.0D0*(S*Z-2.0D0*S*T*Z-S*S*Z-S*Z*Z)
       S5Z=256.0D0*(S*T-2.0D0*S*T*Z-S*S*T-S*T*T)
C
       SHAPE(   1) = 1.0D0-S-T-Z-1.0D0/4.0D0*S5
       DERIV(1, 1) =-1.0D0      -1.0D0/4.0D0*S5S
       DERIV(2, 1) =-1.0D0      -1.0D0/4.0D0*S5T
       DERIV(3, 1) =-1.0D0      -1.0D0/4.0D0*S5Z
C
       SHAPE(   2) = S          -1.0D0/4.0D0*S5
       DERIV(1, 2) = 1.0D0      -1.0D0/4.0D0*S5S
       DERIV(2, 2) = 0.0D0      -1.0D0/4.0D0*S5T
       DERIV(3, 2) = 0.0D0      -1.0D0/4.0D0*S5Z
C
       SHAPE(   3) = T          -1.0D0/4.0D0*S5
       DERIV(1, 3) = 0.0D0      -1.0D0/4.0D0*S5S
       DERIV(2, 3) = 1.0D0      -1.0D0/4.0D0*S5T
       DERIV(3, 3) = 0.0D0      -1.0D0/4.0D0*S5Z
C
       SHAPE(   4) = Z          -1.0D0/4.0D0*S5
       DERIV(1, 4) = 0.0D0      -1.0D0/4.0D0*S5S
       DERIV(2, 4) = 0.0D0      -1.0D0/4.0D0*S5T
       DERIV(3, 4) = 1.0D0      -1.0D0/4.0D0*S5Z
C
       SHAPE(   5) =                         S5
       DERIV(1, 5) =                         S5S
       DERIV(2, 5) =                         S5T
       DERIV(3, 5) =                         S5Z
C
       RETURN
      ENDIF
C
C**** 10 NODED TETRAHEDRAL ELEMENT
C
      IF(NNODE.EQ.10) THEN
C
       A1=1.0D0-S-T-Z
       A2=S
       A3=T
       A4=Z
C
       SHAPE(   1) =(2.0D0*A1-1.0D0)*A1
       DERIV(1, 1) = 1.0D0-4.0D0*A1
       DERIV(2, 1) = 1.0D0-4.0D0*A1
       DERIV(3, 1) = 1.0D0-4.0D0*A1
C
       SHAPE(   2) =(2.0D0*A2-1.0D0)*A2
       DERIV(1, 2) = 4.0D0*A2-1.0D0
       DERIV(2, 2) = 0.0D0
       DERIV(3, 2) = 0.0D0
C
       SHAPE(   3) =(2.0D0*A3-1.0D0)*A3
       DERIV(1, 3) = 0.0D0
       DERIV(2, 3) =-1.0D0+4.0D0*A3
       DERIV(3, 3) = 0.0D0
C
       SHAPE(   4) =(2.0D0*A4-1.0D0)*A4
       DERIV(1, 4) = 0.0D0
       DERIV(2, 4) = 0.0D0
       DERIV(3, 4) = 4.0D0*A4-1.0D0
C
       SHAPE(   5) = 4.0D0*A1*A2
       DERIV(1, 5) = 4.0D0*(A1-A2)
       DERIV(2, 5) =-4.0D0*A2
       DERIV(3, 5) =-4.0D0*A2
C
       SHAPE(   6) = 4.0D0*A2*A3
       DERIV(1, 6) = 4.0D0*A3
       DERIV(2, 6) = 4.0D0*A2
       DERIV(3, 6) = 0.0D0
C
       SHAPE(   7) = 4.0D0*A1*A3
       DERIV(1, 7) =-4.0D0*A3
       DERIV(2, 7) = 4.0D0*(A1-A3)
       DERIV(3, 7) =-4.0D0*A3
C
       SHAPE(   8) = 4.0D0*A1*A4
       DERIV(1, 8) =-4.0D0*A4
       DERIV(2, 8) =-4.0D0*A4
       DERIV(3, 8) = 4.0D0*(A1-A4)
C
       SHAPE(   9) = 4.0D0*A2*A4
       DERIV(1, 9) = 4.0D0*A4
       DERIV(2, 9) = 0.0D0
       DERIV(3, 9) = 4.0D0*A2
C
       SHAPE(  10) = 4.0D0*A3*A4
       DERIV(1,10) = 0.0D0
       DERIV(2,10) = 4.0D0*A4
       DERIV(3,10) = 4.0D0*A3
C
       RETURN
      ENDIF
C
C**** TOBLERONE ELEMENT
C
      IF(NNODE.EQ.6) THEN
C
       ZM = 1.0D0-Z
       ZP = 1.0D0+Z
C
       SHAPE(   1) = 0.5D0*(1.0D0-S-T)*ZM
       DERIV(1, 1) =-0.5D0*ZM
       DERIV(2, 1) =-0.5D0*ZM
       DERIV(3, 1) =-0.5D0*(1.0D0-S-T)
C
       SHAPE(   2) = 0.5D0*S*ZM
       DERIV(1, 2) = 0.5D0*ZM
       DERIV(2, 2) = 0.0D0
       DERIV(3, 2) =-0.5D0*S
C
       SHAPE(   3) = 0.5D0*T*ZM
       DERIV(1, 3) = 0.0D0
       DERIV(2, 3) = 0.5D0*ZM
       DERIV(3, 3) =-0.5D0*T
C
       SHAPE(   4) = 0.5D0*(1.0D0-S-T)*ZP
       DERIV(1, 4) =-0.5D0*ZP
       DERIV(2, 4) =-0.5D0*ZP
       DERIV(3, 4) = 0.5D0*(1.0D0-S-T)
C
       SHAPE(   5) = 0.5D0*S*ZP
       DERIV(1, 5) = 0.5D0*ZP
       DERIV(2, 5) = 0.0D0
       DERIV(3, 5) = 0.5D0*S
C
       SHAPE(   6) = 0.5D0*T*ZP
       DERIV(1, 6) = 0.0D0
       DERIV(2, 6) = 0.5D0*ZP
       DERIV(3, 6) = 0.5D0*T
C
       RETURN
      ENDIF
C
      SM = 0.5D0*(1.0D0-S)
      TM = 0.5D0*(1.0D0-T)
      ZM = 0.5D0*(1.0D0-Z)
      SP = 0.5D0*(1.0D0+S)
      TP = 0.5D0*(1.0D0+T)
      ZP = 0.5D0*(1.0D0+Z)
C
C**** 8 NODED BRICK ELEMENT
C
      IF(NNODE.EQ.8) THEN
C
       SHAPE(   1) = +     SM*TM*ZM
       DERIV(1, 1) = -.5D0   *TM*ZM
       DERIV(2, 1) = -.5D0*SM   *ZM
       DERIV(3, 1) = -.5D0*SM*TM
C
       SHAPE(   2) = +     SP*TM*ZM
       DERIV(1, 2) = +.5D0   *TM*ZM
       DERIV(2, 2) = -.5D0*SP   *ZM
       DERIV(3, 2) = -.5D0*SP*TM
C
       SHAPE(   3) = +     SP*TP*ZM
       DERIV(1, 3) = +.5D0   *TP*ZM
       DERIV(2, 3) = +.5D0*SP   *ZM
       DERIV(3, 3) = -.5D0*SP*TP
C
       SHAPE(   4) = +     SM*TP*ZM
       DERIV(1, 4) = -.5D0   *TP*ZM
       DERIV(2, 4) = +.5D0*SM   *ZM
       DERIV(3, 4) = -.5D0*SM*TP
C
       SHAPE(   5) = +     SM*TM*ZP
       DERIV(1, 5) = -.5D0   *TM*ZP
       DERIV(2, 5) = -.5D0*SM   *ZP
       DERIV(3, 5) = +.5D0*SM*TM
C
       SHAPE(   6) = +     SP*TM*ZP 
       DERIV(1, 6) = +.5D0   *TM*ZP
       DERIV(2, 6) = -.5D0*SP   *ZP
       DERIV(3, 6) = +.5D0*SP*TM
C
       SHAPE(   7) = +     SP*TP*ZP
       DERIV(1, 7) = +.5D0   *TP*ZP
       DERIV(2, 7) = +.5D0*SP   *ZP
       DERIV(3, 7) = +.5D0*SP*TP
C
       SHAPE(   8) = +     SM*TP*ZP
       DERIV(1, 8) = -.5D0   *TP*ZP
       DERIV(2, 8) = +.5D0*SM   *ZP
       DERIV(3, 8) = +.5D0*SM*TP
C
       RETURN
      ENDIF
C
C**** 20 NODED BRICK ELEMENTS
C
      IF(NNODE.EQ.20) THEN
       SS = 0.5D0*(1.0D0-S*S)
       TT = 0.5D0*(1.0D0-T*T)
       ZZ = 0.5D0*(1.0D0-Z*Z)
       S2 = -2.0D0*S
       T2 = -2.0D0*T
       Z2 = -2.0D0*Z
C
       AA          = +SM+TM+ZM-   2.5D0
       SHAPE(   1) = +SM*TM*ZM*(AA+AA)
       DERIV(1, 1) = -   TM*ZM*(AA+SM)
       DERIV(2, 1) = -SM   *ZM*(AA+TM)
       DERIV(3, 1) = -SM*TM   *(AA+ZM)
C
       AA          = +SP+TM+ZM-   2.5D0
       SHAPE(   2) = +SP*TM*ZM*(AA+AA)
       DERIV(1, 2) = +   TM*ZM*(AA+SP)
       DERIV(2, 2) = -SP   *ZM*(AA+TM)
       DERIV(3, 2) = -SP*TM   *(AA+ZM)
C
       AA          = +SP+TP+ZM-   2.5D0
       SHAPE(   3) = +SP*TP*ZM*(AA+AA)
       DERIV(1, 3) = +   TP*ZM*(AA+SP)
       DERIV(2, 3) = +SP   *ZM*(AA+TP)
       DERIV(3, 3) = -SP*TP   *(AA+ZM)
C
       AA          = +SM+TP+ZM-   2.5D0
       SHAPE(   4) = +SM*TP*ZM*(AA+AA)
       DERIV(1, 4) = -   TP*ZM*(AA+SM)
       DERIV(2, 4) = +SM   *ZM*(AA+TP)
       DERIV(3, 4) = -SM*TP   *(AA+ZM)
C
       AA          = +SM+TM+ZP-   2.5D0
       SHAPE(   5) = +SM*TM*ZP*(AA+AA)
       DERIV(1, 5) = -   TM*ZP*(AA+SM)
       DERIV(2, 5) = -SM   *ZP*(AA+TM)
       DERIV(3, 5) = +SM*TM   *(AA+ZP)
C
       AA          = +SP+TM+ZP-   2.5D0
       SHAPE(   6) = +SP*TM*ZP*(AA+AA)
       DERIV(1, 6) = +   TM*ZP*(AA+SP)
       DERIV(2, 6) = -SP   *ZP*(AA+TM)
       DERIV(3, 6) = +SP*TM   *(AA+ZP)
C
       AA          = +SP+TP+ZP-   2.5D0
       SHAPE(   7) = +SP*TP*ZP*(AA+AA)
       DERIV(1, 7) = +   TP*ZP*(AA+SP)
       DERIV(2, 7) = +SP   *ZP*(AA+TP)
       DERIV(3, 7) = +SP*TP   *(AA+ZP)
C
       AA          = +SM+TP+ZP-   2.5D0
       SHAPE(   8) = +SM*TP*ZP*(AA+AA)
       DERIV(1, 8) = -   TP*ZP*(AA+SM)
       DERIV(2, 8) = +SM   *ZP*(AA+TP)
       DERIV(3, 8) = +SM*TP   *(AA+ZP)
C
       SHAPE(   9) = +SS*TM*ZM*2.0D0
       DERIV(1, 9) = +S2*TM*ZM
       DERIV(2, 9) = -SS   *ZM
       DERIV(3, 9) = -SS*TM
C
       SHAPE(  10) = +SP*TT*ZM*2.0D0
       DERIV(1,10) = +   TT*ZM
       DERIV(2,10) = +SP*T2*ZM
       DERIV(3,10) = -SP*TT
C
       SHAPE(  11) = +SS*TP*ZM*2.0D0
       DERIV(1,11) = +S2*TP*ZM
       DERIV(2,11) = +SS   *ZM
       DERIV(3,11) = -SS*TP
C
       SHAPE(  12) = +SM*TT*ZM*2.0D0
       DERIV(1,12) = -   TT*ZM
       DERIV(2,12) = +SM*T2*ZM
       DERIV(3,12) = -SM*TT
C
       SHAPE(  13) = +SS*TM*ZP*2.0D0
       DERIV(1,13) = +S2*TM*ZP
       DERIV(2,13) = -SS   *ZP
       DERIV(3,13) = +SS*TM
C
       SHAPE(  14) = +SP*TT*ZP*2.0D0
       DERIV(1,14) = +   TT*ZP
       DERIV(2,14) = +SP*T2*ZP
       DERIV(3,14) = +SP*TT
C
       SHAPE(  15) = +SS*TP*ZP*2.0D0
       DERIV(1,15) = +S2*TP*ZP
       DERIV(2,15) = +SS   *ZP
       DERIV(3,15) = +SS*TP
C
       SHAPE(  16) = +SM*TT*ZP*2.0D0
       DERIV(1,16) = -   TT*ZP
       DERIV(2,16) = +SM*T2*ZP
       DERIV(3,16) = +SM*TT
C
       SHAPE(  17) = +SM*TM*ZZ*2.0D0
       DERIV(1,17) = -   TM*ZZ
       DERIV(2,17) = -SM   *ZZ
       DERIV(3,17) = +SM*TM*Z2
C
       SHAPE(  18) = +SP*TM*ZZ*2.0D0
       DERIV(1,18) = +   TM*ZZ
       DERIV(2,18) = -SP   *ZZ
       DERIV(3,18) = +SP*TM*Z2
C
       SHAPE(  19) = +SP*TP*ZZ*2.0D0
       DERIV(1,19) = +   TP*ZZ
       DERIV(2,19) = +SP   *ZZ
       DERIV(3,19) = +SP*TP*Z2
C
       SHAPE(  20) = +SM*TP*ZZ*2.0D0
       DERIV(1,20) = -   TP*ZZ
       DERIV(2,20) = +SM   *ZZ
       DERIV(3,20) = +SM*TP*Z2
C
       RETURN
      ENDIF
C
C**** 27 NODED BRICK ELEMENT
C
      IF(NNODE.EQ.27) THEN                             ! to be revised !
       SL=S*(S-1.0D0)
       TL=T*(T-1.0D0)
       ZL=Z*(Z-1.0D0)
       SP=S*(S+1.0D0)
       TP=T*(T+1.0D0)
       ZP=Z*(Z+1.0D0)
       S1=2.0D0*S-1.0D0
       T1=2.0D0*T-1.0D0
       Z1=2.0D0*Z-1.0D0
       S2=1.0D0-S*S
       T2=1.0D0-T*T
       Z2=1.0D0-Z*Z
       S3=1.0D0+2.0D0*S
       T3=1.0D0+2.0D0*T
       Z3=1.0D0+2.0D0*Z
       S4=-2.0D0*S
       T4=-2.0D0*T
       Z4=-2.0D0*Z
C
       SHAPE(   1) = .125D0*SL*TL*ZL
       DERIV(1, 1) = .125D0*S1*TL*ZL
       DERIV(2, 1) = .125D0*SL*T1*ZL
       DERIV(3, 1) = .125D0*SL*TL*Z1
C
       SHAPE(   2) = .125D0*SP*TL*ZL
       DERIV(1, 2) = .125D0*S3*TL*ZL
       DERIV(2, 2) = .125D0*SP*T1*ZL
       DERIV(3, 2) = .125D0*SP*TL*Z1
C
       SHAPE(   3) = .125D0*SP*TP*ZL
       DERIV(1, 3) = .125D0*S3*TP*ZL
       DERIV(2, 3) = .125D0*SP*T3*ZL
       DERIV(3, 3) = .125D0*SP*TP*Z1
C
       SHAPE(   4) = .125D0*SL*TP*ZL
       DERIV(1, 4) = .125D0*S1*TP*ZL
       DERIV(2, 4) = .125D0*SL*T3*ZL
       DERIV(3, 4) = .125D0*SL*TP*Z1
C
       SHAPE(   5) = .125D0*SL*TL*ZP
       DERIV(1, 5) = .125D0*S1*TL*ZP
       DERIV(2, 5) = .125D0*SL*T1*ZP
       DERIV(3, 5) = .125D0*SL*TL*Z3
C
       SHAPE(   6) = .125D0*SP*TL*ZP
       DERIV(1, 6) = .125D0*S3*TL*ZP
       DERIV(2, 6) = .125D0*SP*T1*ZP
       DERIV(3, 6) = .125D0*SP*TL*Z3
C
       SHAPE(   7) = .125D0*SP*TP*ZP
       DERIV(1, 7) = .125D0*S3*TP*ZP
       DERIV(2, 7) = .125D0*SP*T3*ZP
       DERIV(3, 7) = .125D0*SP*TP*Z3
C
       SHAPE(   8) = .125D0*SL*TP*ZP
       DERIV(1, 8) = .125D0*S1*TP*ZP
       DERIV(2, 8) = .125D0*SL*T3*ZP
       DERIV(3, 8) = .125D0*SL*TP*Z3
C
       SHAPE(   9) = .25D0*S2*TL*ZL
       DERIV(1, 9) = .25D0*S4*TL*ZL
       DERIV(2, 9) = .25D0*S2*T1*ZL
       DERIV(3, 9) = .25D0*S2*TL*Z1
C
       SHAPE(  10) = .25D0*SP*T2*ZL
       DERIV(1,10) = .25D0*S3*T2*ZL
       DERIV(2,10) = .25D0*SP*T4*ZL
       DERIV(3,10) = .25D0*SP*T2*Z1
C
       SHAPE(  11) = .25D0*S2*TP*ZL
       DERIV(1,11) = .25D0*S4*TP*ZL
       DERIV(2,11) = .25D0*S2*T3*ZL
       DERIV(3,11) = .25D0*S2*TP*Z1
C
       SHAPE(  12) = .25D0*SL*T2*ZL
       DERIV(1,12) = .25D0*S1*T2*ZL
       DERIV(2,12) = .25D0*SL*T4*ZL
       DERIV(3,12) = .25D0*SL*T2*Z1
C
       SHAPE(  13) = .25D0*SL*TL*Z2
       DERIV(1,13) = .25D0*S1*TL*Z2
       DERIV(2,13) = .25D0*SL*T1*Z2
       DERIV(3,13) = .25D0*SL*TL*Z4
C
       SHAPE(  14) = .25D0*SP*TL*Z2
       DERIV(1,14) = .25D0*S3*TL*Z2
       DERIV(2,14) = .25D0*SP*T1*Z2
       DERIV(3,14) = .25D0*SP*TL*Z4
C
       SHAPE(  15) = .25D0*SP*TP*Z2
       DERIV(1,15) = .25D0*S3*TP*Z2
       DERIV(2,15) = .25D0*SP*T3*Z2
       DERIV(3,15) = .25D0*SP*TP*Z4
C
       SHAPE(  16) = .25D0*SL*TP*Z2
       DERIV(1,16) = .25D0*S1*TP*Z2
       DERIV(2,16) = .25D0*SL*T3*Z2
       DERIV(3,16) = .25D0*SL*TP*Z4
C
       SHAPE(  17) = .25D0*S2*TL*ZP
       DERIV(1,17) = .25D0*S4*TL*ZP
       DERIV(2,17) = .25D0*S2*T1*ZP
       DERIV(3,17) = .25D0*S2*TL*Z3
C
       SHAPE(  18) = .25D0*SP*T2*ZP
       DERIV(1,18) = .25D0*S3*T2*ZP
       DERIV(2,18) = .25D0*SP*T4*ZP
       DERIV(3,18) = .25D0*SP*T2*Z3
C
       SHAPE(  19) = .25D0*S2*TP*ZP
       DERIV(1,19) = .25D0*S4*TP*ZP
       DERIV(2,19) = .25D0*S2*T3*ZP
       DERIV(3,19) = .25D0*S2*TP*Z3
C
       SHAPE(  20) = .25D0*SL*T2*ZP
       DERIV(1,20) = .25D0*S1*T2*ZP
       DERIV(2,20) = .25D0*SL*T4*ZP
       DERIV(3,20) = .25D0*SL*T2*Z3
C 
       SHAPE(  21) = .5D0*S2*T2*ZL
       DERIV(1,21) = .5D0*S4*T2*ZL
       DERIV(2,21) = .5D0*S2*T4*ZL
       DERIV(3,21) = .5D0*S2*T2*Z1
C
       SHAPE(  22) = .5D0*S2*TL*Z2
       DERIV(1,22) = .5D0*S4*TL*Z2
       DERIV(2,22) = .5D0*S2*T1*Z2
       DERIV(3,22) = .5D0*S2*TL*Z4
C
       SHAPE(  23) = .5D0*SP*T2*Z2
       DERIV(1,23) = .5D0*S3*T2*Z2
       DERIV(2,23) = .5D0*SP*T4*Z2
       DERIV(3,23) = .5D0*SP*T2*Z4
C
       SHAPE(  24) = .5D0*S2*TP*Z2
       DERIV(1,24) = .5D0*S4*TP*Z2
       DERIV(2,24) = .5D0*S2*T3*Z2
       DERIV(3,24) = .5D0*S2*TP*Z4
C
       SHAPE(  25) = .5D0*SL*T2*Z2
       DERIV(1,25) = .5D0*S1*T2*Z2
       DERIV(2,25) = .5D0*SL*T4*Z2
       DERIV(3,25) = .5D0*SL*T2*Z4
C
       SHAPE(  26) = .5D0*S2*T2*ZP
       DERIV(1,26) = .5D0*S4*T2*ZP
       DERIV(2,26) = .5D0*S2*T4*ZP
       DERIV(3,26) = .5D0*S2*T2*Z3
C
       SHAPE(  27) = S2*T2*Z2
       DERIV(1,27) = S4*T2*Z2
       DERIV(2,27) = S2*T4*Z2
       DERIV(3,27) = S2*T2*Z4
C
       RETURN
      ENDIF
C
      END
