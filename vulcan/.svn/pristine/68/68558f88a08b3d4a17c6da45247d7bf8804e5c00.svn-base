      SUBROUTINE INMICRO7(ALPHAM,INUPC,IPLAT)
C***********************************************************************
C
C**** THIS ROUTINE INITIALITES MICROSTRUCTURAL PARAMETERS FOR MODEL 7
C
C***********************************************************************
C
C     Index of variables
C
C     ALPHAM= array of microstructural (microscopical) variables
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES (thermal-microstructural)
C
      INCLUDE 'nued_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION ALPHAM(NHISTM)
C
C**** TRANSFERS "VPLAT" TO MICROSCOPICAL VARIABLES
C
      IVERSI=INT(VPLAT(IPLAT,6))
C
      IF(IVERSI.EQ.1) THEN
      AM0=VPLAT(IPLAT,7)
      AK=VPLAT(IPLAT,8)
      VISC7=VPLAT(IPLAT,9)
      UNIVC=VPLAT(IPLAT,10)
      ENERA=VPLAT(IPLAT,11)
      HENER=VPLAT(IPLAT,12)      ! constant
      CHININF=VPLAT(IPLAT,13)
      CONTA=VPLAT(IPLAT,14)
      CONTB=VPLAT(IPLAT,15)
      IMECHMIC=INT(VPLAT(IPLAT,16))
C
      IV=16                      ! number of VPLAT defined in input data
C
C**** INITIALITES MICROSTRUCTURAL PARAMETERS ("ALPHAM" & OTHERS)
C
      IN=INUPC   
      ALPHAM(IN+1)=0.0           ! hydration degree
      ALPHAM(IN+2)=0.0           ! compression strength
      ALPHAM(IN+3)=0.0           ! -- (Young modulus?)
      ALPHAM(IN+4)=0.0           ! hydration amplitude
C
C**** INCREMENTS "ALPHAM" INDEX
C
      INUPC=INUPC+4
      ENDIF                      ! iversi.eq.1
C
      IF(IVERSI.EQ.2) THEN 
      AKBAR=VPLAT(IPLAT,7)
      ETA00=VPLAT(IPLAT,8)
      ETABAR=VPLAT(IPLAT,9)
      GHINF=VPLAT(IPLAT,10)
      ALPHABAR=VPLAT(IPLAT,11)
      AABAR=VPLAT(IPLAT,12)
      UNIVC=VPLAT(IPLAT,13)
      ENERA=VPLAT(IPLAT,14)
      HENER=VPLAT(IPLAT,15)      ! constant
      CONTA=VPLAT(IPLAT,16)
      CONTB=VPLAT(IPLAT,17)
      ENETE=VPLAT(IPLAT,18)
      FINFIN=VPLAT(IPLAT,19)
      FINFINM=VPLAT(IPLAT,20)
      GSETT=VPLAT(IPLAT,21)
      TMAX=VPLAT(IPLAT,22)
      TREF=VPLAT(IPLAT,23)
      CONAE=VPLAT(IPLAT,24)
      IMECHMIC=INT(VPLAT(IPLAT,25))
C
      IV=25                      ! number of VPLAT defined in input data
C
C**** INITIALITES MICROSTRUCTURAL PARAMETERS ("ALPHAM" & OTHERS)
C
      IN=INUPC
      ALPHAM(IN+1)=0.0           ! hydration degree
      ALPHAM(IN+2)=0.0           ! compression strength
      ALPHAM(IN+3)=0.0           ! tension strength
      ALPHAM(IN+4)=0.0           ! Young modulus
      ALPHAM(IN+5)=0.0           ! rate of hydration degree
      ALPHAM(IN+6)=0.0           ! hydration amplitude
      ALPHAM(IN+7)=0.0           ! chemical afinity
      ALPHAM(IN+8)=0.0           ! chemical dissipation
C
C**** INCREMENTS "ALPHAM" INDEX
C
      INUPC=INUPC+8
      ENDIF                      ! iversi.eq.2
C
      RETURN
      END
