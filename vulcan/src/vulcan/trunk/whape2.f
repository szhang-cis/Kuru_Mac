      SUBROUTINE WHAPE2(DERIV,S,T,NDIME,NNODE,SHAPE,HACHE,
     .                  XJACM,IPERT,XJACI,CARTT,SHAPT)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES PERTURBATIONS FUNCTIONS AND THEIR
C     DERIVATIVES FOR LINEAR AND QUADRATIC ISOPARAMETRIC 2-D ELEMENTS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION DERIV(NDIME,*), SHAPE(*), HACHE(*), XJACM(NDIME,*)
C
      DIMENSION XJACI(NDIME,*), CARTT(NDIME,*), SHAPT(NDIME,*)
C
      ST=S*T
C
C**** 3 NODED ELEMENT
C
      IF(NNODE.EQ.3) THEN
       GO TO (1,2,3) IPERT
C
    1  SHAPE(1)=27*S*T*(1.0-S-T)*HACHE(1)
       SHAPE(2)=27*S*T*(1.0-S-T)*HACHE(2)
       SHAPE(3)=27*S*T*(1.0-S-T)*HACHE(3)
C
       DERIV(1,1)=27*(T-T*T-2*S*T)*HACHE(1)
       DERIV(1,2)=27*(T-T*T-2*S*T)*HACHE(2)
       DERIV(1,3)=27*(T-T*T-2*S*T)*HACHE(3)
       DERIV(2,1)=27*(S-S*S-2*S*T)*HACHE(1)
       DERIV(2,2)=27*(S-S*S-2*S*T)*HACHE(2)
       DERIV(2,3)=27*(S-S*S-2*S*T)*HACHE(3)
       RETURN
C
    2  SHAPE(1)=27*S*T*(1.0-S-T)*HACHE(1)
       SHAPE(2)=27*S*T*(1.0-S-T)*HACHE(2)
       SHAPE(3)=27*S*T*(1.0-S-T)*HACHE(3)
C
       DERIV(1,1)=0.0
       DERIV(1,2)=0.0
       DERIV(1,3)=0.0
       DERIV(2,1)=0.0
       DERIV(2,2)=0.0
       DERIV(2,3)=0.0
       RETURN
C
    3  CONTINUE
C
C**** CALCULATE THE DETERMINANT AND INVERS OF JACOBIAN MATRIX
C
       CALL INVMTX(XJACM,XJACI,DETJM,NDIME)
C
       SHAPT(1,1)=-1.0
       SHAPT(1,2)= 1.0
       SHAPT(1,3)= 0.0
       SHAPT(2,1)=-1.0
       SHAPT(2,2)= 0.0
       SHAPT(2,3)= 1.0
C
C**** CALCULATE THE CARTESIAN DERIVATIVES CORRESPONDING TO "SHAPT"
C
       CALL PROMA1(CARTT,XJACI,SHAPT,NDIME,NNODE,NDIME)
       CALL PROMA3(SHAPE,CARTT,HACHE,NNODE,NDIME,    1)
C
       DERIV(1,1)=0.0
       DERIV(1,2)=0.0
       DERIV(1,3)=0.0
       DERIV(2,1)=0.0
       DERIV(2,2)=0.0
       DERIV(2,3)=0.0
       RETURN
      ENDIF
C
C**** 4 NODED ELEMENT
C
      IF(NNODE.EQ.4) THEN
       GO TO (11,12,13) IPERT
C
   11  A1=1.0
       SHAPE(1)=A1*(1.-S*S)*(1.-T*T)*HACHE(1)
       SHAPE(2)=A1*(1.-S*S)*(1.-T*T)*HACHE(2)
       SHAPE(3)=A1*(1.-S*S)*(1.-T*T)*HACHE(3)
       SHAPE(4)=A1*(1.-S*S)*(1.-T*T)*HACHE(4)
C
       DERIV(1,1)=-A1*2.0*S*(1.-T*T)*HACHE(1)
       DERIV(1,2)=-A1*2.0*S*(1.-T*T)*HACHE(2)
       DERIV(1,3)=-A1*2.0*S*(1.-T*T)*HACHE(3)
       DERIV(1,4)=-A1*2.0*S*(1.-T*T)*HACHE(4)
       DERIV(2,1)=-A1*2.0*(1.-S*S)*T*HACHE(1)
       DERIV(2,2)=-A1*2.0*(1.-S*S)*T*HACHE(2)
       DERIV(2,3)=-A1*2.0*(1.-S*S)*T*HACHE(3)
       DERIV(2,4)=-A1*2.0*(1.-S*S)*T*HACHE(4)
       RETURN
C
   12  A1=1.0
       SHAPE(1)=A1*(1.-S*S)*(1.-T*T)*HACHE(1)
       SHAPE(2)=A1*(1.-S*S)*(1.-T*T)*HACHE(2)
       SHAPE(3)=A1*(1.-S*S)*(1.-T*T)*HACHE(3)
       SHAPE(4)=A1*(1.-S*S)*(1.-T*T)*HACHE(4)
C
       DERIV(1,1)=0.0
       DERIV(1,2)=0.0
       DERIV(1,3)=0.0
       DERIV(1,4)=0.0
       DERIV(2,1)=0.0
       DERIV(2,2)=0.0
       DERIV(2,3)=0.0
       DERIV(2,4)=0.0
       RETURN
C
   13  CONTINUE
C
C**** CALCULATE THE DETERMINANT AND INVERS OF JACOBIAN MATRIX
C
       CALL INVMTX(XJACM,XJACI,DETJM,NDIME)
C
       SHAPT(1,1)=(-1.+T)*0.25
       SHAPT(1,2)=(+1.-T)*0.25
       SHAPT(1,3)=(+1.+T)*0.25
       SHAPT(1,4)=(-1.-T)*0.25
       SHAPT(2,1)=(-1.+S)*0.25
       SHAPT(2,2)=(-1.-S)*0.25
       SHAPT(2,3)=(+1.+S)*0.25
       SHAPT(2,4)=(+1.-S)*0.25
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
       RETURN
      ENDIF
C
C**** 6 NODED ELEMENT
C
      IF(NNODE.EQ.6) THEN
        call runendt('whape2: nnode=6 ')
        A1=1.0-S-T
        A2=S
        A3=T
C
        SHAPE(1)=(2.0*A1-1.0)*A1
        SHAPE(2)=(2.0*A2-1.0)*A2
        SHAPE(3)=(2.0*A3-1.0)*A3
        SHAPE(4)=4.0*A1*A2
        SHAPE(5)=4.0*A2*A3
        SHAPE(6)=4.0*A1*A3
C
        DERIV(1,1)=1.0-4.0*A1
        DERIV(1,2)=4.0*A2-1.0
        DERIV(1,3)=0.0
        DERIV(1,4)=4.0*(A1-A2)
        DERIV(1,5)=4.0*A3
        DERIV(1,6)=-4.0*A3
        DERIV(2,1)=1.0-4.0*A1
        DERIV(2,2)=0.0
        DERIV(2,3)=4.0*A3-1.0
        DERIV(2,4)=-4.0*A2
        DERIV(2,5)=4.0*A2
        DERIV(2,6)=4.0*(A1-A3)
C
        RETURN
      ENDIF
C
C**** 8 NODED ELEMENT
C
      IF(NNODE.EQ.8) THEN
        call runendt('whape2: nnode=8 ')
        S2=S*2.0
        T2=T*2.0
        SS=S*S
        TT=T*T
        ST=S*T
        SST=S*S*T
        STT=S*T*T
        ST2=S*T*2.0
C
        SHAPE(1)=(-1.0+ST+SS+TT-SST-STT)/4.0
        SHAPE(2)=(-1.0-ST+SS+TT-SST+STT)/4.0
        SHAPE(3)=(-1.0+ST+SS+TT+SST+STT)/4.0
        SHAPE(4)=(-1.0-ST+SS+TT+SST-STT)/4.0
        SHAPE(5)=(1.0-T-SS+SST)/2.0
        SHAPE(6)=(1.0+S-TT-STT)/2.0
        SHAPE(7)=(1.0+T-SS-SST)/2.0
        SHAPE(8)=(1.0-S-TT+STT)/2.0
C
        DERIV(1,1)=(T+S2-ST2-TT)/4.0
        DERIV(1,2)=(-T+S2-ST2+TT)/4.0
        DERIV(1,3)=(T+S2+ST2+TT)/4.0
        DERIV(1,4)=(-T+S2+ST2-TT)/4.0
        DERIV(1,5)=-S+ST
        DERIV(1,6)=(1.0-TT)/2.0
        DERIV(1,7)=-S-ST
        DERIV(1,8)=(-1.0+TT)/2.0
        DERIV(2,1)=(S+T2-SS-ST2)/4.0
        DERIV(2,2)=(-S+T2-SS+ST2)/4.0
        DERIV(2,3)=(S+T2+SS+ST2)/4.0
        DERIV(2,4)=(-S+T2+SS-ST2)/4.0
        DERIV(2,5)=(-1.0+SS)/2.0
        DERIV(2,6)=-T-ST
        DERIV(2,7)=(1.0-SS)/2.0
        DERIV(2,8)=-T+ST
C
        RETURN
      ENDIF
C
C**** 9 NODED ELEMENT
C
      IF(NNODE.EQ.9) THEN
        call runendt('whape2: nnode=9 ')
        SS=S*S
        ST=S*T
        TT=T*T
        S1=S+1.0
        T1=T+1.0
        S2=S*2.0
        T2=T*2.0
        S9=S-1.0
        T9=T-1.0
C
        SHAPE(1)=0.25*S9*ST*T9
        SHAPE(2)=0.25*S1*ST*T9
        SHAPE(3)=0.25*S1*ST*T1  
        SHAPE(4)=0.25*S9*ST*T1
        SHAPE(5)=0.5*(1.0-SS)*T*T9
        SHAPE(6)=0.5*S*S1*(1.0-TT)
        SHAPE(7)=0.5*(1.0-SS)*T*T1
        SHAPE(8)=0.5*S*S9*(1.0-TT)
        SHAPE(9)=(1.0-SS)*(1.0-TT)
C
        DERIV(1,1)=0.25*T*T9*(-1.0+S2)
        DERIV(1,2)=0.25*(1.0+S2)*T*T9
        DERIV(1,3)=0.25*(1.0+S2)*T*T1
        DERIV(1,4)=0.25*(-1.0+S2)*T*T1
        DERIV(1,5)=-ST*T9
        DERIV(1,6)=0.5*(1.0+S2)*(1.0-TT)
        DERIV(1,7)=-ST*T1
        DERIV(1,8)=0.5*(-1.0+S2)*(1.0-TT)
        DERIV(1,9)=-S2*(1.0-TT)
        DERIV(2,1)=0.25*(-1.0+T2)*S*S9
        DERIV(2,2)=0.25*S*S1*(-1.0+T2)
        DERIV(2,3)=0.25*S*S1*(1.0+T2)
        DERIV(2,4)=0.25*S*S9*(1.0+T2)
        DERIV(2,5)=0.5*(1.0-SS)*(-1.0+T2)
        DERIV(2,6)=-ST*S1
        DERIV(2,7)=0.5*(1.0-SS)*(1.0+T2)
        DERIV(2,8)=-ST*S9
        DERIV(2,9)=-T2*(1.0-SS)
C
      RETURN
      ENDIF
C
      END
