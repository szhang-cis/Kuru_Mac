      SUBROUTINE SHAPE2(DERIV,S,T,NDIME,NNODE,NQUTR,NNODX,SHAPE)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES SHAPE FUNCTIONS AND THEIR DERIVATIVES FOR
C     LINEAR AND QUADRATIC ISOPARAMETRIC 2-D ELEMENTS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DERIV(NDIME,*),SHAPE(*)
C
      ST=S*T
C
C**** 3 NODED ELEMENT
C
      IF(NNODE.EQ.3) THEN
       IF(NQUTR.EQ.1) THEN                               ! triangle
        SHAPE(1+NNODX)=1.0D0-S-T
        SHAPE(2+NNODX)=S
        SHAPE(3+NNODX)=T
C
        DERIV(1,1+NNODX)=-1.0D0
        DERIV(1,2+NNODX)= 1.0D0
        DERIV(1,3+NNODX)= 0.0D0
        DERIV(2,1+NNODX)=-1.0D0
        DERIV(2,2+NNODX)= 0.0D0
        DERIV(2,3+NNODX)= 1.0D0
       ENDIF
C
       IF(NQUTR.EQ.2) THEN                               ! quadrilateral
c       SHAPE(1)=(1.0D0-S)*(1.5D0-0.5D0*T)*0.25D0
c       SHAPE(2)=(1.0D0+S)*(1.0D0-T)*0.25D0
c       SHAPE(3)=(1.5D0+0.5D0*S)*(1.0D0+T)*0.25D0
C
c       DERIV(1,1)=(-1.5D0+0.5D0*T)*0.25D0
c       DERIV(1,2)=(+1.0D0-T)*0.25D0
c       DERIV(1,3)=(+1.0D0+T)*0.125D0
c       DERIV(2,1)=(-1.0D0+S)*0.125D0
c       DERIV(2,2)=(-1.0D0-S)*0.25D0
c       DERIV(2,3)=(+1.5D0+0.5D0*S)*0.25D0


        SHAPE(1)=(1.0D0-T-S+ST)*0.25D0   ! other option
        SHAPE(2)=(1.0D0-T+S-ST)*0.25D0
        SHAPE(3)=(1.0D0+T     )*0.50D0
C
        DERIV(1,1)=(-1.0D0+T)*0.25D0
        DERIV(1,2)=(+1.0D0-T)*0.25D0
        DERIV(1,3)=           0.00D0
        DERIV(2,1)=(-1.0D0+S)*0.25D0
        DERIV(2,2)=(-1.0D0-S)*0.25D0
        DERIV(2,3)=           0.50D0



       ENDIF
C
       RETURN
      ENDIF
C
C**** 4 NODED ELEMENT
C
      IF(NNODE.EQ.4) THEN
       IF(NQUTR.EQ.1) THEN                               ! quadrilateral
        SHAPE(1+NNODX)=(1.0D0-T-S+ST)*0.25D0
        SHAPE(2+NNODX)=(1.0D0-T+S-ST)*0.25D0
        SHAPE(3+NNODX)=(1.0D0+T+S+ST)*0.25D0
        SHAPE(4+NNODX)=(1.0D0+T-S-ST)*0.25D0
C
        DERIV(1,1+NNODX)=(-1.0D0+T)*0.25D0
        DERIV(1,2+NNODX)=(+1.0D0-T)*0.25D0
        DERIV(1,3+NNODX)=(+1.0D0+T)*0.25D0
        DERIV(1,4+NNODX)=(-1.0D0-T)*0.25D0
        DERIV(2,1+NNODX)=(-1.0D0+S)*0.25D0
        DERIV(2,2+NNODX)=(-1.0D0-S)*0.25D0
        DERIV(2,3+NNODX)=(+1.0D0+S)*0.25D0
        DERIV(2,4+NNODX)=(+1.0D0-S)*0.25D0
       ENDIF
C
       IF(NQUTR.EQ.2) THEN                               ! triangle
        S4=27.0D0*(1.0D0-S-T)*S*T
        SHAPE(1)=1.0D0-S-T-1.0D0/3.0D0*S4
        SHAPE(2)=      S  -1.0D0/3.0D0*S4
        SHAPE(3)=        T-1.0D0/3.0D0*S4
        SHAPE(4)=S4
C
        S4S=27.0D0*(T-2.0D0*S*T-T*T)
        S4T=27.0D0*(S-2.0D0*S*T-S*S)
        DERIV(1,1)=-1.0D0-1.0D0/3.0D0*S4S
        DERIV(1,2)= 1.0D0-1.0D0/3.0D0*S4S
        DERIV(1,3)=      -1.0D0/3.0D0*S4S
        DERIV(1,4)=                   S4S
        DERIV(2,1)=-1.0D0-1.0D0/3.0D0*S4T
        DERIV(2,2)=      -1.0D0/3.0D0*S4T
        DERIV(2,3)= 1.0D0-1.0D0/3.0D0*S4T
        DERIV(2,4)=                   S4T
       ENDIF
C
       RETURN
      ENDIF
C
C**** 5 NODED ELEMENT
C
      IF(NNODE.EQ.5) THEN
       S5=(1.0D0-S*S)*(1.0D0-T*T)
       SHAPE(1)=(1.0D0-T-S+ST)*0.25D0-S5*0.25D0
       SHAPE(2)=(1.0D0-T+S-ST)*0.25D0-S5*0.25D0
       SHAPE(3)=(1.0D0+T+S+ST)*0.25D0-S5*0.25D0
       SHAPE(4)=(1.0D0+T-S-ST)*0.25D0-S5*0.25D0
       SHAPE(5)=S5
C
       S5S=-2.0D0*S*(1.0D0-T*T)
       S5T=-2.0D0*T*(1.0D0-S*S)
       DERIV(1,1)=(-1.0D0+T)*0.25D0-S5S*0.25D0
       DERIV(1,2)=(+1.0D0-T)*0.25D0-S5S*0.25D0
       DERIV(1,3)=(+1.0D0+T)*0.25D0-S5S*0.25D0
       DERIV(1,4)=(-1.0D0-T)*0.25D0-S5S*0.25D0
       DERIV(1,5)=                  S5S
       DERIV(2,1)=(-1.0D0+S)*0.25D0-S5T*0.25D0
       DERIV(2,2)=(-1.0D0-S)*0.25D0-S5T*0.25D0
       DERIV(2,3)=(+1.0D0+S)*0.25D0-S5T*0.25D0
       DERIV(2,4)=(+1.0D0-S)*0.25D0-S5T*0.25D0
       DERIV(2,5)=                  S5T
C
       RETURN
      ENDIF
C
C**** 6 NODED ELEMENT
C
      IF(NNODE.EQ.6) THEN
       A1=1.0D0-S-T
       A2=S
       A3=T
C
       SHAPE(1)=(2.0D0*A1-1.0D0)*A1
       SHAPE(2)=(2.0D0*A2-1.0D0)*A2
       SHAPE(3)=(2.0D0*A3-1.0D0)*A3
       SHAPE(4)=4.0D0*A1*A2
       SHAPE(5)=4.0D0*A2*A3
       SHAPE(6)=4.0D0*A1*A3
C
       DERIV(1,1)=1.0D0-4.0D0*A1
       DERIV(1,2)=4.0D0*A2-1.0D0
       DERIV(1,3)=0.0D0
       DERIV(1,4)=4.0D0*(A1-A2)
       DERIV(1,5)=4.0D0*A3
       DERIV(1,6)=-4.0D0*A3
       DERIV(2,1)=1.0D0-4.0D0*A1
       DERIV(2,2)=0.0D0
       DERIV(2,3)=4.0D0*A3-1.0D0
       DERIV(2,4)=-4.0D0*A2
       DERIV(2,5)=4.0D0*A2
       DERIV(2,6)=4.0D0*(A1-A3)
C
       RETURN
      ENDIF
C
C**** 8 NODED ELEMENT
C
      IF(NNODE.EQ.8) THEN
       S2=S*2.0D0
       T2=T*2.0D0
       SS=S*S
       TT=T*T
       ST=S*T
       SST=S*S*T
       STT=S*T*T
       ST2=S*T*2.0D0
C
       SHAPE(1)=(-1.0D0+ST+SS+TT-SST-STT)/4.0D0
       SHAPE(2)=(-1.0D0-ST+SS+TT-SST+STT)/4.0D0
       SHAPE(3)=(-1.0D0+ST+SS+TT+SST+STT)/4.0D0
       SHAPE(4)=(-1.0D0-ST+SS+TT+SST-STT)/4.0D0
       SHAPE(5)=(1.0D0-T-SS+SST)/2.0D0
       SHAPE(6)=(1.0D0+S-TT-STT)/2.0D0
       SHAPE(7)=(1.0D0+T-SS-SST)/2.0D0
       SHAPE(8)=(1.0D0-S-TT+STT)/2.0D0
C
       DERIV(1,1)=(T+S2-ST2-TT)/4.0D0
       DERIV(1,2)=(-T+S2-ST2+TT)/4.0D0
       DERIV(1,3)=(T+S2+ST2+TT)/4.0D0
       DERIV(1,4)=(-T+S2+ST2-TT)/4.0D0
       DERIV(1,5)=-S+ST
       DERIV(1,6)=(1.0D0-TT)/2.0D0
       DERIV(1,7)=-S-ST
       DERIV(1,8)=(-1.0D0+TT)/2.0D0
       DERIV(2,1)=(S+T2-SS-ST2)/4.0D0
       DERIV(2,2)=(-S+T2-SS+ST2)/4.0D0
       DERIV(2,3)=(S+T2+SS+ST2)/4.0D0
       DERIV(2,4)=(-S+T2+SS-ST2)/4.0D0
       DERIV(2,5)=(-1.0D0+SS)/2.0D0
       DERIV(2,6)=-T-ST
       DERIV(2,7)=(1.0D0-SS)/2.0D0
       DERIV(2,8)=-T+ST
C
       RETURN
      ENDIF
C
C**** 9 NODED ELEMENT
C
      IF(NNODE.EQ.9) THEN
       SS=S*S
       ST=S*T
       TT=T*T
       S1=S+1.0D0
       T1=T+1.0D0
       S2=S*2.0D0
       T2=T*2.0D0
       S9=S-1.0D0
       T9=T-1.0D0
C
       SHAPE(1)=0.25D0*S9*ST*T9
       SHAPE(2)=0.25D0*S1*ST*T9
       SHAPE(3)=0.25D0*S1*ST*T1  
       SHAPE(4)=0.25D0*S9*ST*T1
       SHAPE(5)=0.5D0*(1.0D0-SS)*T*T9
       SHAPE(6)=0.5D0*S*S1*(1.0D0-TT)
       SHAPE(7)=0.5D0*(1.0D0-SS)*T*T1
       SHAPE(8)=0.5D0*S*S9*(1.0D0-TT)
       SHAPE(9)=(1.0D0-SS)*(1.0D0-TT)
C
       DERIV(1,1)=0.25D0*T*T9*(-1.0D0+S2)
       DERIV(1,2)=0.25D0*(1.0D0+S2)*T*T9
       DERIV(1,3)=0.25D0*(1.0D0+S2)*T*T1
       DERIV(1,4)=0.25D0*(-1.0D0+S2)*T*T1
       DERIV(1,5)=-ST*T9
       DERIV(1,6)=0.5D0*(1.0D0+S2)*(1.0D0-TT)
       DERIV(1,7)=-ST*T1
       DERIV(1,8)=0.5D0*(-1.0D0+S2)*(1.0D0-TT)
       DERIV(1,9)=-S2*(1.0D0-TT)
       DERIV(2,1)=0.25D0*(-1.0D0+T2)*S*S9
       DERIV(2,2)=0.25D0*S*S1*(-1.0D0+T2)
       DERIV(2,3)=0.25D0*S*S1*(1.0D0+T2)
       DERIV(2,4)=0.25D0*S*S9*(1.0D0+T2)
       DERIV(2,5)=0.5D0*(1.0D0-SS)*(-1.0D0+T2)
       DERIV(2,6)=-ST*S1
       DERIV(2,7)=0.5D0*(1.0D0-SS)*(1.0D0+T2)
       DERIV(2,8)=-ST*S9
       DERIV(2,9)=-T2*(1.0D0-SS)
C
       RETURN
      ENDIF
C
      END
