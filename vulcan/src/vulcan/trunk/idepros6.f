      SUBROUTINE IDEPROS6(PROPST,IPLAT,IA1,INUPM)
C***********************************************************************
C
C**** THIS ROUTINE ORDERS THE MICROSTRUCTURAL PROPERTIES OF MODEL
C     NUMBER 6 (IPCMO=6) OF RATE PHASE-CHANGE FORMULATIONS
C
C***********************************************************************
c
c     IPLAT= number of model. 1 is first model of input, 2 second model
c            of input.
c
c
c


      IMPLICIT REAL*8 (A-H,O-Z)
    
c     COUPLING VARIABLES
c     thermal-microestructural
      INCLUDE 'nued_om.f'
      

c     THERMAL VARIABLES
      INCLUDE 'auxl_omt.f'
      INCLUDE 'prob_omt.f'
     
      DIMENSION PROPST(*)
      IA2=IA1+2                 ! 2=ipcfo,ipcmo
      IVERSI=INT(PROPST(IA2+1))

c*****************************************************************
c 
c     VERSION 1 (SOLIDUS-SOLIDUS PHASE-CHANGE MODEL)
c     
c*****************************************************************
      IF(IVERSI.EQ.1) THEN
         TEMSS=PROPST(IA2+2)
         AKCOO=PROPST(IA2+3)
         AKC11=PROPST(IA2+4)
         AKC22=PROPST(IA2+5)
         ANC11=PROPST(IA2+6)
         HENER=PROPST(IA2+7)
c     transfer to VPLAT array
         VPLAT(IPLAT,6)=REAL(IVERSI)
         VPLAT(IPLAT,7)=TEMSS
         VPLAT(IPLAT,8)=AKCOO
         VPLAT(IPLAT,9)=AKC11
         VPLAT(IPLAT,10)=AKC22
         VPLAT(IPLAT,11)=ANC11
         VPLAT(IPLAT,12)=HENER
c     total number of properties of model 6 v1 (number of
c     elements of VPLAT array)     
         IMODE=7
      ENDIF

c*****************************************************************
c 
c     VERSION 2
c     
c*****************************************************************
      IF(IVERSI.EQ.2) THEN
         TEMFF=PROPST(IA2+2)
         TEMPP=PROPST(IA2+3)
         CONUF=PROPST(IA2+4)
         CONUP=PROPST(IA2+5)
         COGRF=PROPST(IA2+6)
         COGRP=PROPST(IA2+7)
         HENER=PROPST(IA2+8)
C     transfer to VPLAT array
         VPLAT(IPLAT,6)=REAL(IVERSI)
         VPLAT(IPLAT,7)=TEMFF
         VPLAT(IPLAT,8)=TEMPP
         VPLAT(IPLAT,9)=CONUF
         VPLAT(IPLAT,10)=CONUP
         VPLAT(IPLAT,11)=COGRF
         VPLAT(IPLAT,12)=COGRP
         VPLAT(IPLAT,13)=HENER
c     total number of properties of model 6 v2 (number of
c     elements of VPLAT array)
         IMODE=8
      ENDIF

c*****************************************************************
c 
c     VERSION 3 (Austempered Ductile Iron model using
c     Avrami equation)
c     
c*****************************************************************      
      IF(IVERSI.EQ.3) THEN
c     number of models that provide initial conditions
         NCOPC=IPLAC(INUPM,1,1)
c     if there isn't initial condition of other model (NCOPC=0), micros6
c     uses the next parameters
         IF(NCOPC.EQ.0) THEN
            C=PROPST(IA2+2)
            SI=PROPST(IA2+3)
            CU=PROPST(IA2+4)
            AMN=PROPST(IA2+5)
            ANI=PROPST(IA2+6)
            AMO=PROPST(IA2+7)
            CR=PROPST(IA2+8)
            FAINI=PROPST(IA2+9)
            TAUS=PROPST(IA2+10)
            TBS=PROPST(IA2+11)
            TBF=PROPST(IA2+12)
            CCFB=PROPST(IA2+13)
            AKAVRA=PROPST(IA2+14)
            AMAVRA=PROPST(IA2+15)
         endif

c     rest of the parameters load from input file
         IF(NCOPC.NE.0) THEN
            TBS=PROPST(IA2+2)
            TBF=PROPST(IA2+3)
            CCFB=PROPST(IA2+4)
            AKAVRA=PROPST(IA2+5)
            AMAVRA=PROPST(IA2+6)
         ENDIF

c     transfer parameters to VPLAT array
         VPLAT(IPLAT,6)=REAL(IVERSI)
         IF(NCOPC.EQ.0) THEN
            VPLAT(IPLAT,7)=C
            VPLAT(IPLAT,8)=SI
            VPLAT(IPLAT,9)=CU
            VPLAT(IPLAT,10)=AMN
            VPLAT(IPLAT,11)=ANI
            VPLAT(IPLAT,12)=AMO
            VPLAT(IPLAT,13)=CR
            VPLAT(IPLAT,14)=FAINI
            VPLAT(IPLAT,15)=TAUS
         endif

c     for couple model (parameter from micros 14)
c     only for two models
         IF(NCOPC.NE.0) THEN
            NMODC=IPLAC(INUPM,2,1)
            VPLAT(IPLAT,7)=VPLAT(NMODC,IPLAC(INUPM,2,13)) !c
            VPLAT(IPLAT,8)=VPLAT(NMODC,IPLAC(INUPM,2,14)) !si
            VPLAT(IPLAT,9)=VPLAT(NMODC,IPLAC(INUPM,2,18)) !cu
            VPLAT(IPLAT,10)=VPLAT(NMODC,IPLAC(INUPM,2,15)) !mn
            VPLAT(IPLAT,11)=VPLAT(NMODC,IPLAC(INUPM,2,17)) !ni
            VPLAT(IPLAT,12)=VPLAT(NMODC,IPLAC(INUPM,2,19)) !mo
            VPLAT(IPLAT,13)=VPLAT(NMODC,IPLAC(INUPM,2,16)) !cr
            VPLAT(IPLAT,14)=0.D0
            VPLAT(IPLAT,15)=0.D0
         endif
         VPLAT(IPLAT,16)=TBS
         VPLAT(IPLAT,17)=TBF
         VPLAT(IPLAT,18)=CCFB
         VPLAT(IPLAT,19)=AKAVRA
         VPLAT(IPLAT,20)=AMAVRA

c     total number of properties of model 6 v3 (number of
c     elements of VPLAT array)
         IF(NCOPC.EQ.0) THEN
            IMODE=18
         ENDIF
         IF(NCOPC.NE.0) THEN
            IMODE=9
         ENDIF
      ENDIF
      IA1=IA2+IMODE

      RETURN
      END
