      SUBROUTINE INMICRO6(ALPHAM,INUPC,IPLAT,INUPM)
C***********************************************************************
C
C**** THIS ROUTINE INITIALITES MICROSTRUCTURAL PARAMETERS FOR MODEL 6
C
C***********************************************************************
C
C     Index of variables
C
C     ALPHAM= array of microstructural (microscopical) variables
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

c--   COUPLING VARIABLES (thermal-microstructural)
      INCLUDE 'nued_om.f'
    
c--   THERMAL VARIABLES
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
      
      DIMENSION ALPHAM(NHISTM)

c--   transfer from VPLAT to microscopical variables
      IVERSI=INT(VPLAT(IPLAT,6))

c*****************************************************************
c 
c     VERSION 1 (SOLIDUS-SOLIDUS PHASE-CHANGE MODEL)
c     
c*****************************************************************
c--   transfer from VPLAT to microscopical variables
      IF(IVERSI.EQ.1) THEN
         TEMSS=VPLAT(IPLAT,7)
         AKCOO=VPLAT(IPLAT,8)
         AKC11=VPLAT(IPLAT,9)
         AKC22=VPLAT(IPLAT,10)
         ANC11=VPLAT(IPLAT,11)
         HENER=VPLAT(IPLAT,12)  ! constant
c     number of VPLAT defined in input data
         IV=12

c--   initiallites microstructural parameters (ALPHAM and others)
         IN=INUPC
         ALPHAM(IN+1)=1.0
         ALPHAM(IN+2)=0.0
         ALPHAM(IN+3)=0.0
         ALPHAM(IN+4)=TEMSS
c   increments ALPHAM index        
         INUPC=INUPC+4
      ENDIF                 

c*****************************************************************
c 
c     VERSION 2
c     
c*****************************************************************
c--   transfer from VPLAT to microscopical variables
      IF(IVERSI.EQ.2) THEN
         TEMFF=VPLAT(IPLAT,7)
         TEMPP=VPLAT(IPLAT,8)
         CONUF=VPLAT(IPLAT,9)
         CONUP=VPLAT(IPLAT,10)
         COGRF=VPLAT(IPLAT,11)
         COGRP=VPLAT(IPLAT,12)
         HENER=VPLAT(IPLAT,13)
c     number of VPLAT defined in input data
         IV=13

c--   initiallites microstructural parameters (ALPHAM and others)         
         IN=INUPC
         ALPHAM(IN+1)=1.0       ! solid-solid phase-change function
         ALPHAM(IN+2)=0.0       ! ferrite fraction
         ALPHAM(IN+3)=0.0       ! grain density of ferrite
         ALPHAM(IN+4)=0.0       ! ferrite radius
         ALPHAM(IN+5)=0.0       ! pearlite fraction
         ALPHAM(IN+6)=0.0       ! grain density of pearlite
         ALPHAM(IN+7)=0.0       ! pearlite radius
c     increments ALPHAM index        
         INUPC=INUPC+7
      ENDIF                 

c*****************************************************************
c 
c     VERSION 3 (Austempered Ductile Iron model using
c     Avrami equation)
c     
c*****************************************************************
      IF(IVERSI.EQ.3) THEN

c--   number of models that provide initial conditions
         NCOPC=IPLAC(INUPM,1,1)

c--   if parameters are load from input file (not coupled model)
         IF(NCOPC.EQ.0) THEN

c--   transfer from VPLAT to microscopical variables
            C=VPLAT(IPLAT,7)
            SI=VPLAT(IPLAT,8)
            CU=VPLAT(IPLAT,9)
            AMN=VPLAT(IPLAT,10)
            ANI=VPLAT(IPLAT,11)
            AMO=VPLAT(IPLAT,12)
            CR=VPLAT(IPLAT,13)
            FAINI=VPLAT(IPLAT,14)
            TAUS=VPLAT(IPLAT,15)
            TBS=VPLAT(IPLAT,16)
            TBF=VPLAT(IPLAT,17)
            CCFB=VPLAT(IPLAT,18)
            AKAVRA=VPLAT(IPLAT,19)
            AMAVRA=VPLAT(IPLAT,20)
c     number of VPLAT defined in input data
            IV=20

c--   initialize microstructural parameters (ALPHAM and others)
            IN=INUPC
            ALPHAM(IN+1)=0.D0   !vfBF
            ALPHAM(IN+2)=FAINI  !vfA
            cA=-0.435D0+3.35D-4*TAUS+
     .           1.61D-6*TAUS*TAUS+0.006D0*AMN-0.11D0*SI-0.07D0*ANI+
     .           0.014D0*CU-0.3D0*AMO
c     cA and cAR initial aren't equal to zero because there's problem when
c     this variables go to mechanical problem. I think it is for smoothing
c     funtion, because there are a big gradient when cA and cAR are
c     initialized to zero and suddenly is equal to cA. Is correct initialized
c     these variables to cA.
            ALPHAM(IN+3)=cA     !cA initial
            ALPHAM(IN+4)=cA     !cAR initial
            vfGr=1.D0-FAINI
            ALPHAM(IN+5)=vfGr   !vfGr initial
            ALPHAM(IN+6)=0.D0   !httime
            INUPC=INUPC+6       !increments ALPHAM index
         ENDIF

c--   coupled model (micros 6 and 14)
         IF(NCOPC.NE.0) THEN
            IN=INUPC
            ALPHAM(IN+1)=0.D0   !vfBF
            ALPHAM(IN+2)=0.D0   !vfA
            ALPHAM(IN+3)=0.D0   !cA
            ALPHAM(IN+4)=0.D0   !cAR
            ALPHAM(IN+5)=0.D0   !vfGr
            ALPHAM(IN+6)=0.D0   !httime
            INUPC=INUPC+6       !increments ALPHAM index
         ENDIF
      ENDIF
      
      RETURN
      END
