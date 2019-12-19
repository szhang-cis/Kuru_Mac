      SUBROUTINE WHAPE3(DERIV,S,T,Z,NDIME,NNODE,SHAPE,HACHE,
     .                  XJACM,IPERT,XJACI,CARTT,SHAPT)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES PERTURBATION FUNCTIONS AND THEIR
C     DERIVATIVES FOR LINEAR AND QUADRATIC ISOPARAMETRIC 3-D ELEMENTS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION DERIV(NDIME,*), SHAPE(*), HACHE(*), XJACM(NDIME,*)
C
      DIMENSION XJACI(NDIME,*), CARTT(NDIME,*), SHAPT(NDIME,*)
C
C**** 4 NODED TETRAHEDRAL ELEMENT
C
      IF(NNODE.EQ.4) THEN
       GO TO (1,2,3) IPERT
C
    1  CONTINUE
       call runendt('error in whape3')
C 
       SHAPE(   1) = 1.-S-T-Z
       DERIV(1, 1) =-1.
       DERIV(2, 1) =-1.
       DERIV(3, 1) =-1.
C
       SHAPE(   2) = S
       DERIV(1, 2) = 1.
       DERIV(2, 2) = 0.
       DERIV(3, 2) = 0.
C
       SHAPE(   3) = T
       DERIV(1, 3) = 0.
       DERIV(2, 3) = 1.
       DERIV(3, 3) = 0.
C
       SHAPE(   4) = Z
       DERIV(1, 4) = 0.
       DERIV(2, 4) = 0.
       DERIV(3, 4) = 1.
       RETURN
C
    2  CONTINUE
       call runendt('error in whape3: nnode=4')
       RETURN
C
    3  CONTINUE
C
C**** CALCULATE THE DETERMINANT AND INVERS OF JACOBIAN MATRIX
C
       CALL INVMTX(XJACM,XJACI,DETJM,NDIME)
C
       SHAPT(1, 1) =-1.
       SHAPT(2, 1) =-1.
       SHAPT(3, 1) =-1.
C
       SHAPT(1, 2) = 1.
       SHAPT(2, 2) = 0.
       SHAPT(3, 2) = 0.
C
       SHAPT(1, 3) = 0.
       SHAPT(2, 3) = 1.
       SHAPT(3, 3) = 0.
C
       SHAPT(1, 4) = 0.
       SHAPT(2, 4) = 0.
       SHAPT(3, 4) = 1.
C
C**** CALCULATE THE CARTESIAN DERIVATIVES CORRESPONDING TO "SHAPT"
C
       CALL PROMA1(CARTT,XJACI,SHAPT,NDIME,NNODE,NDIME)
       CALL PROMA3(SHAPE,CARTT,HACHE,NNODE,NDIME,    1)
C
       DERIV(1,1)=0.0
       DERIV(1,2)=0.0
       DERIV(1,3)=0.0
       DERIV(1,4)=0.0
       DERIV(2,1)=0.0
       DERIV(2,2)=0.0
       DERIV(2,3)=0.0
       DERIV(2,4)=0.0
       DERIV(3,1)=0.0
       DERIV(3,2)=0.0
       DERIV(3,3)=0.0
       DERIV(3,4)=0.0
       RETURN
      ENDIF
C
C**** 10 NODED TETRAHEDRAL ELEMENT
C
      IF(NNODE.EQ.10) THEN
       GO TO (11,12,13) IPERT
C
   11  CONTINUE
       call runendt('error in whape3')
C
       A1=1.-S-T-Z
       A2=S
       A3=T
       A4=Z
C
       SHAPE(   1) =(2.*A1-1.)*A1
       DERIV(1, 1) = 1.-4.*A1
       DERIV(2, 1) = 1.-4.*A1
       DERIV(3, 1) = 1.-4.*A1
C
       SHAPE(   2) =(2.*A2-1.)*A2
       DERIV(1, 2) = 4.*A2-1.
       DERIV(2, 2) = 0.
       DERIV(3, 2) = 0.
C
       SHAPE(   3) =(2.*A3-1.)*A3
       DERIV(1, 3) = 0.
       DERIV(2, 3) =-1.+4.*A3
       DERIV(3, 3) = 0.
C
       SHAPE(   4) =(2.*A4-1.)*A4
       DERIV(1, 4) = 0.
       DERIV(2, 4) = 0.
       DERIV(3, 4) = 4.*A4-1.
C
       SHAPE(   5) = 4.*A1*A2
       DERIV(1, 5) = 4.*(A1-A2)
       DERIV(2, 5) =-4.*A2
       DERIV(3, 5) =-4.*A2
C
       SHAPE(   6) = 4.*A2*A3
       DERIV(1, 6) = 4.*A3
       DERIV(2, 6) = 4.*A2
       DERIV(3, 6) = 0.
C
       SHAPE(   7) = 4.*A1*A3
       DERIV(1, 7) =-4.*A3
       DERIV(2, 7) = 4.*(A1-A3)
       DERIV(3, 7) =-4.*A3
C
       SHAPE(   8) = 4.*A1*A4
       DERIV(1, 8) =-4.*A4
       DERIV(2, 8) =-4.*A4
       DERIV(3, 8) = 4.*(A1-A4)
C
       SHAPE(   9) = 4.*A2*A4
       DERIV(1, 9) = 4.*A4
       DERIV(2, 9) = 0.
       DERIV(3, 9) = 4.*A2
C
       SHAPE(  10) = 4.*A3*A4
       DERIV(1,10) = 0.
       DERIV(2,10) = 4.*A4
       DERIV(3,10) = 4.*A3
       RETURN
C
   12  CONTINUE
       call runendt('error in whape3')
       RETURN
C
   13  CONTINUE
       call runendt('error in whape3')
       RETURN
      ENDIF
C
C**** TOBLERONE ELEMENT
C
      IF(NNODE.EQ.6) THEN
       GO TO (21,22,23) IPERT
C
   21  CONTINUE
       call runendt('error in whape3')
C
       ZM = 1.-Z
       ZP = 1.+Z
C
       SHAPE(   1) = 0.5*(1.-S-T)*ZM
       DERIV(1, 1) =-0.5*ZM
       DERIV(2, 1) =-0.5*ZM
       DERIV(3, 1) =-0.5*(1.-S-T)
C
       SHAPE(   2) = 0.5*S*ZM
       DERIV(1, 2) = 0.5*ZM
       DERIV(2, 2) = 0.
       DERIV(3, 2) =-0.5*S
C
       SHAPE(   3) = 0.5*T*ZM
       DERIV(1, 3) = 0.
       DERIV(2, 3) = 0.5*ZM
       DERIV(3, 3) =-0.5*T
C
       SHAPE(   4) = 0.5*(1.-S-T)*ZP
       DERIV(1, 4) =-0.5*ZP
       DERIV(2, 4) =-0.5*ZP
       DERIV(3, 4) = 0.5*(1.-S-T)
C
       SHAPE(   5) = 0.5*S*ZP
       DERIV(1, 5) = 0.5*ZP
       DERIV(2, 5) = 0.
       DERIV(3, 5) = 0.5*S
C
       SHAPE(   6) = 0.5*T*ZP
       DERIV(1, 6) = 0.
       DERIV(2, 6) = 0.5*ZP
       DERIV(3, 6) = 0.5*T
       RETURN
C
   22  CONTINUE
       call runendt('error in whape3')
       RETURN
C
   23  CONTINUE
       call runendt('error in whape3')
       RETURN
      ENDIF
C
      SM = 0.5*(1.-S)
      TM = 0.5*(1.-T)
      ZM = 0.5*(1.-Z)
      SP = 0.5*(1.+S)
      TP = 0.5*(1.+T)
      ZP = 0.5*(1.+Z)
C
C**** 8 NODED BRICK ELEMENT
C
      IF(NNODE.EQ.8) THEN
       GO TO (31,32,33) IPERT
C
   31  CONTINUE
       call runendt('error in whape3')
C
       SHAPE(   1) = +   SM*TM*ZM
       DERIV(1, 1) = -.5   *TM*ZM
       DERIV(2, 1) = -.5*SM   *ZM
       DERIV(3, 1) = -.5*SM*TM
C
       SHAPE(   2) = +   SP*TM*ZM
       DERIV(1, 2) = +.5   *TM*ZM
       DERIV(2, 2) = -.5*SP   *ZM
       DERIV(3, 2) = -.5*SP*TM
C
       SHAPE(   3) = +   SP*TP*ZM
       DERIV(1, 3) = +.5   *TP*ZM
       DERIV(2, 3) = +.5*SP   *ZM
       DERIV(3, 3) = -.5*SP*TP
C
       SHAPE(   4) = +   SM*TP*ZM
       DERIV(1, 4) = -.5   *TP*ZM
       DERIV(2, 4) = +.5*SM   *ZM
       DERIV(3, 4) = -.5*SM*TP
C
       SHAPE(   5) = +   SM*TM*ZP
       DERIV(1, 5) = -.5   *TM*ZP
       DERIV(2, 5) = -.5*SM   *ZP
       DERIV(3, 5) = +.5*SM*TM
C
       SHAPE(   6) = +   SP*TM*ZP 
       DERIV(1, 6) = +.5   *TM*ZP
       DERIV(2, 6) = -.5*SP   *ZP
       DERIV(3, 6) = +.5*SP*TM
C
       SHAPE(   7) = +   SP*TP*ZP
       DERIV(1, 7) = +.5   *TP*ZP
       DERIV(2, 7) = +.5*SP   *ZP
       DERIV(3, 7) = +.5*SP*TP
C
       SHAPE(   8) = +   SM*TP*ZP
       DERIV(1, 8) = -.5   *TP*ZP
       DERIV(2, 8) = +.5*SM   *ZP
       DERIV(3, 8) = +.5*SM*TP
       RETURN
C
   32  CONTINUE
       call runendt('error in whape3')
       RETURN
C
   33  CONTINUE
C
C**** CALCULATE THE DETERMINANT AND INVERS OF JACOBIAN MATRIX
C
       CALL INVMTX(XJACM,XJACI,DETJM,NDIME)
C
       SHAPT(1, 1) = -.5   *TM*ZM
       SHAPT(2, 1) = -.5*SM   *ZM
       SHAPT(3, 1) = -.5*SM*TM
C
       SHAPT(1, 2) = +.5   *TM*ZM
       SHAPT(2, 2) = -.5*SP   *ZM
       SHAPT(3, 2) = -.5*SP*TM
C
       SHAPT(1, 3) = +.5   *TP*ZM
       SHAPT(2, 3) = +.5*SP   *ZM
       SHAPT(3, 3) = -.5*SP*TP
C
       SHAPT(1, 4) = -.5   *TP*ZM
       SHAPT(2, 4) = +.5*SM   *ZM
       SHAPT(3, 4) = -.5*SM*TP
C
       SHAPT(1, 5) = -.5   *TM*ZP
       SHAPT(2, 5) = -.5*SM   *ZP
       SHAPT(3, 5) = +.5*SM*TM
C
       SHAPT(1, 6) = +.5   *TM*ZP
       SHAPT(2, 6) = -.5*SP   *ZP
       SHAPT(3, 6) = +.5*SP*TM
C
       SHAPT(1, 7) = +.5   *TP*ZP
       SHAPT(2, 7) = +.5*SP   *ZP
       SHAPT(3, 7) = +.5*SP*TP
C
       SHAPT(1, 8) = -.5   *TP*ZP
       SHAPT(2, 8) = +.5*SM   *ZP
       SHAPT(3, 8) = +.5*SM*TP
C
C**** CALCULATE THE CARTESIAN DERIVATIVES CORRESPONDING TO "SHAPT"
C
       CALL PROMA1(CARTT,XJACI,SHAPT,NDIME,NNODE,NDIME)
C
       DO INODE=1,NNODE
        SHAPE(INODE)=0.0
        DO IDIME=1,NDIME
         SHAPE(INODE)=SHAPE(INODE)+CARTT(IDIME,INODE)*HACHE(IDIME)
        ENDDO
       ENDDO
C
       DERIV(1,1)=0.0
       DERIV(1,2)=0.0
       DERIV(1,3)=0.0
       DERIV(1,4)=0.0
       DERIV(1,5)=0.0
       DERIV(1,6)=0.0
       DERIV(1,7)=0.0
       DERIV(1,8)=0.0
       DERIV(2,1)=0.0
       DERIV(2,2)=0.0
       DERIV(2,3)=0.0
       DERIV(2,4)=0.0
       DERIV(2,5)=0.0
       DERIV(2,6)=0.0
       DERIV(2,7)=0.0
       DERIV(2,8)=0.0
       DERIV(3,1)=0.0
       DERIV(3,2)=0.0
       DERIV(3,3)=0.0
       DERIV(3,4)=0.0
       DERIV(3,5)=0.0
       DERIV(3,6)=0.0
       DERIV(3,7)=0.0
       DERIV(3,8)=0.0
       RETURN
      ENDIF
C
C**** 20 NODED BRICK ELEMENTS
C
      IF(NNODE.EQ.20) THEN
       GO TO (41,42,43) IPERT
C
   41  CONTINUE
       call runendt('error in whape3')
C
       SS = 0.5*(1.-S*S)
       TT = 0.5*(1.-T*T)
       ZZ = 0.5*(1.-Z*Z)
       S2 = -2.*S
       T2 = -2.*T
       Z2 = -2.*Z
C
       AA          = +SM+TM+ZM-   2.5
       SHAPE(   1) = +SM*TM*ZM*(AA+AA)
       DERIV(1, 1) = -   TM*ZM*(AA+SM)
       DERIV(2, 1) = -SM   *ZM*(AA+TM)
       DERIV(3, 1) = -SM*TM   *(AA+ZM)
C
       AA          = +SP+TM+ZM-   2.5
       SHAPE(   2) = +SP*TM*ZM*(AA+AA)
       DERIV(1, 2) = +   TM*ZM*(AA+SP)
       DERIV(2, 2) = -SP   *ZM*(AA+TM)
       DERIV(3, 2) = -SP*TM   *(AA+ZM)
C
       AA          = +SP+TP+ZM-   2.5
       SHAPE(   3) = +SP*TP*ZM*(AA+AA)
       DERIV(1, 3) = +   TP*ZM*(AA+SP)
       DERIV(2, 3) = +SP   *ZM*(AA+TP)
       DERIV(3, 3) = -SP*TP   *(AA+ZM)
C
       AA          = +SM+TP+ZM-   2.5
       SHAPE(   4) = +SM*TP*ZM*(AA+AA)
       DERIV(1, 4) = -   TP*ZM*(AA+SM)
       DERIV(2, 4) = +SM   *ZM*(AA+TP)
       DERIV(3, 4) = -SM*TP   *(AA+ZM)
C
       AA          = +SM+TM+ZP-   2.5
       SHAPE(   5) = +SM*TM*ZP*(AA+AA)
       DERIV(1, 5) = -   TM*ZP*(AA+SM)
       DERIV(2, 5) = -SM   *ZP*(AA+TM)
       DERIV(3, 5) = +SM*TM   *(AA+ZP)
C
       AA          = +SP+TM+ZP-   2.5
       SHAPE(   6) = +SP*TM*ZP*(AA+AA)
       DERIV(1, 6) = +   TM*ZP*(AA+SP)
       DERIV(2, 6) = -SP   *ZP*(AA+TM)
       DERIV(3, 6) = +SP*TM   *(AA+ZP)
C
       AA          = +SP+TP+ZP-   2.5
       SHAPE(   7) = +SP*TP*ZP*(AA+AA)
       DERIV(1, 7) = +   TP*ZP*(AA+SP)
       DERIV(2, 7) = +SP   *ZP*(AA+TP)
       DERIV(3, 7) = +SP*TP   *(AA+ZP)
C
       AA          = +SM+TP+ZP-   2.5
       SHAPE(   8) = +SM*TP*ZP*(AA+AA)
       DERIV(1, 8) = -   TP*ZP*(AA+SM)
       DERIV(2, 8) = +SM   *ZP*(AA+TP)
       DERIV(3, 8) = +SM*TP   *(AA+ZP)
C
       SHAPE(   9) = +SS*TM*ZM*2.
       DERIV(1, 9) = +S2*TM*ZM
       DERIV(2, 9) = -SS   *ZM
       DERIV(3, 9) = -SS*TM
C
       SHAPE(  10) = +SP*TT*ZM*2.
       DERIV(1,10) = +   TT*ZM
       DERIV(2,10) = +SP*T2*ZM
       DERIV(3,10) = -SP*TT
C
       SHAPE(  11) = +SS*TP*ZM*2.
       DERIV(1,11) = +S2*TP*ZM
       DERIV(2,11) = +SS   *ZM
       DERIV(3,11) = -SS*TP
C
       SHAPE(  12) = +SM*TT*ZM*2.
       DERIV(1,12) = -   TT*ZM
       DERIV(2,12) = +SM*T2*ZM
       DERIV(3,12) = -SM*TT
C
       SHAPE(  13) = +SM*TM*ZZ*2.
       DERIV(1,13) = -   TM*ZZ
       DERIV(2,13) = -SM   *ZZ
       DERIV(3,13) = +SM*TM*Z2
C
       SHAPE(  14) = +SP*TM*ZZ*2.
       DERIV(1,14) = +   TM*ZZ
       DERIV(2,14) = -SP   *ZZ
       DERIV(3,14) = +SP*TM*Z2
C
       SHAPE(  15) = +SP*TP*ZZ*2.
       DERIV(1,15) = +   TP*ZZ
       DERIV(2,15) = +SP   *ZZ
       DERIV(3,15) = +SP*TP*Z2
C
       SHAPE(  16) = +SM*TP*ZZ*2.
       DERIV(1,16) = -   TP*ZZ
       DERIV(2,16) = +SM   *ZZ
       DERIV(3,16) = +SM*TP*Z2
C
       SHAPE(  17) = +SS*TM*ZP*2.
       DERIV(1,17) = +S2*TM*ZP
       DERIV(2,17) = -SS   *ZP
       DERIV(3,17) = +SS*TM
C
       SHAPE(  18) = +SP*TT*ZP*2.
       DERIV(1,18) = +   TT*ZP
       DERIV(2,18) = +SP*T2*ZP
       DERIV(3,18) = +SP*TT
C
       SHAPE(  19) = +SS*TP*ZP*2.
       DERIV(1,19) = +S2*TP*ZP
       DERIV(2,19) = +SS   *ZP
       DERIV(3,19) = +SS*TP
C
       SHAPE(  20) = +SM*TT*ZP*2.
       DERIV(1,20) = -   TT*ZP
       DERIV(2,20) = +SM*T2*ZP
       DERIV(3,20) = +SM*TT
       RETURN
C
   42  CONTINUE
       call runendt('error in whape3')
       RETURN
C
   43  CONTINUE
       call runendt('error in whape3')
       RETURN
      ENDIF
C
C**** 27 NODED BRICK ELEMENT
C
      IF(NNODE.EQ.27) THEN
       GO TO (51,52,53) IPERT
C
   51  CONTINUE
       call runendt('error in whape3')
C
       SL=S*(S-1.)
       TL=T*(T-1.)
       ZL=Z*(Z-1.)
       SP=S*(S+1.)
       TP=T*(T+1.)
       ZP=Z*(Z+1.)
       S1=2.*S-1.
       T1=2.*T-1.
       Z1=2.*Z-1.
       S2=1.-S*S
       T2=1.-T*T
       Z2=1.-Z*Z
       S3=1.+2.*S
       T3=1.+2.*T
       Z3=1.+2.*Z
       S4=-2.*S
       T4=-2.*T
       Z4=-2.*Z
C
       SHAPE(   1) = .125*SL*TL*ZL
       DERIV(1, 1) = .125*S1*TL*ZL
       DERIV(2, 1) = .125*SL*T1*ZL
       DERIV(3, 1) = .125*SL*TL*Z1
C
       SHAPE(   2) = .125*SP*TL*ZL
       DERIV(1, 2) = .125*S3*TL*ZL
       DERIV(2, 2) = .125*SP*T1*ZL
       DERIV(3, 2) = .125*SP*TL*Z1
C
       SHAPE(   3) = .125*SP*TP*ZL
       DERIV(1, 3) = .125*S3*TP*ZL
       DERIV(2, 3) = .125*SP*T3*ZL
       DERIV(3, 3) = .125*SP*TP*Z1
C
       SHAPE(   4) = .125*SL*TP*ZL
       DERIV(1, 4) = .125*S1*TP*ZL
       DERIV(2, 4) = .125*SL*T3*ZL
       DERIV(3, 4) = .125*SL*TP*Z1
C
       SHAPE(   5) = .125*SL*TL*ZP
       DERIV(1, 5) = .125*S1*TL*ZP
       DERIV(2, 5) = .125*SL*T1*ZP
       DERIV(3, 5) = .125*SL*TL*Z3
C
       SHAPE(   6) = .125*SP*TL*ZP
       DERIV(1, 6) = .125*S3*TL*ZP
       DERIV(2, 6) = .125*SP*T1*ZP
       DERIV(3, 6) = .125*SP*TL*Z3
C
       SHAPE(   7) = .125*SP*TP*ZP
       DERIV(1, 7) = .125*S3*TP*ZP
       DERIV(2, 7) = .125*SP*T3*ZP
       DERIV(3, 7) = .125*SP*TP*Z3
C
       SHAPE(   8) = .125*SL*TP*ZP
       DERIV(1, 8) = .125*S1*TP*ZP
       DERIV(2, 8) = .125*SL*T3*ZP
       DERIV(3, 8) = .125*SL*TP*Z3
C
       SHAPE(   9) = .25*S2*TL*ZL
       DERIV(1, 9) = .25*S4*TL*ZL
       DERIV(2, 9) = .25*S2*T1*ZL
       DERIV(3, 9) = .25*S2*TL*Z1
C
       SHAPE(  10) = .25*SP*T2*ZL
       DERIV(1,10) = .25*S3*T2*ZL
       DERIV(2,10) = .25*SP*T4*ZL
       DERIV(3,10) = .25*SP*T2*Z1
C
       SHAPE(  11) = .25*S2*TP*ZL
       DERIV(1,11) = .25*S4*TP*ZL
       DERIV(2,11) = .25*S2*T3*ZL
       DERIV(3,11) = .25*S2*TP*Z1
C
       SHAPE(  12) = .25*SL*T2*ZL
       DERIV(1,12) = .25*S1*T2*ZL
       DERIV(2,12) = .25*SL*T4*ZL
       DERIV(3,12) = .25*SL*T2*Z1
C
       SHAPE(  13) = .25*SL*TL*Z2
       DERIV(1,13) = .25*S1*TL*Z2
       DERIV(2,13) = .25*SL*T1*Z2
       DERIV(3,13) = .25*SL*TL*Z4
C
       SHAPE(  14) = .25*SP*TL*Z2
       DERIV(1,14) = .25*S3*TL*Z2
       DERIV(2,14) = .25*SP*T1*Z2
       DERIV(3,14) = .25*SP*TL*Z4
C
       SHAPE(  15) = .25*SP*TP*Z2
       DERIV(1,15) = .25*S3*TP*Z2
       DERIV(2,15) = .25*SP*T3*Z2
       DERIV(3,15) = .25*SP*TP*Z4
C
       SHAPE(  16) = .25*SL*TP*Z2
       DERIV(1,16) = .25*S1*TP*Z2
       DERIV(2,16) = .25*SL*T3*Z2
       DERIV(3,16) = .25*SL*TP*Z4
C
       SHAPE(  17) = .25*S2*TL*ZP
       DERIV(1,17) = .25*S4*TL*ZP
       DERIV(2,17) = .25*S2*T1*ZP
       DERIV(3,17) = .25*S2*TL*Z3
C
       SHAPE(  18) = .25*SP*T2*ZP
       DERIV(1,18) = .25*S3*T2*ZP
       DERIV(2,18) = .25*SP*T4*ZP
       DERIV(3,18) = .25*SP*T2*Z3
C
       SHAPE(  19) = .25*S2*TP*ZP
       DERIV(1,19) = .25*S4*TP*ZP
       DERIV(2,19) = .25*S2*T3*ZP
       DERIV(3,19) = .25*S2*TP*Z3
C
       SHAPE(  20) = .25*SL*T2*ZP
       DERIV(1,20) = .25*S1*T2*ZP
       DERIV(2,20) = .25*SL*T4*ZP
       DERIV(3,20) = .25*SL*T2*Z3
C
       SHAPE(  21) = .5*S2*T2*ZL
       DERIV(1,21) = .5*S4*T2*ZL
       DERIV(2,21) = .5*S2*T4*ZL
       DERIV(3,21) = .5*S2*T2*Z1
C
       SHAPE(  22) = .5*S2*TL*Z2
       DERIV(1,22) = .5*S4*TL*Z2
       DERIV(2,22) = .5*S2*T1*Z2
       DERIV(3,22) = .5*S2*TL*Z4
C
       SHAPE(  23) = .5*SP*T2*Z2
       DERIV(1,23) = .5*S3*T2*Z2
       DERIV(2,23) = .5*SP*T4*Z2
       DERIV(3,23) = .5*SP*T2*Z4
C
       SHAPE(  24) = .5*S2*TP*Z2
       DERIV(1,24) = .5*S4*TP*Z2
       DERIV(2,24) = .5*S2*T3*Z2
       DERIV(3,24) = .5*S2*TP*Z4
C
       SHAPE(  25) = .5*SL*T2*Z2
       DERIV(1,25) = .5*S1*T2*Z2
       DERIV(2,25) = .5*SL*T4*Z2
       DERIV(3,25) = .5*SL*T2*Z4
C
       SHAPE(  26) = .5*S2*T2*ZP
       DERIV(1,26) = .5*S4*T2*ZP
       DERIV(2,26) = .5*S2*T4*ZP
       DERIV(3,26) = .5*S2*T2*Z3
C
       SHAPE(  27) = S2*T2*Z2
       DERIV(1,27) = S4*T2*Z2
       DERIV(2,27) = S2*T4*Z2
       DERIV(3,27) = S2*T2*Z4
       RETURN
C
   52  CONTINUE
       call runendt('error in whape3')
       RETURN
C
   53  CONTINUE
       call runendt('error in whape3')
       RETURN
      ENDIF
      END
