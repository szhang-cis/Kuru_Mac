      SUBROUTINE IDEPROM7(PROPS,IPLAT,IA1)
C***********************************************************************
C
C**** THIS ROUTINE ORDERS THE MICROSTRUCTURAL PROPERTIES OF MODEL
C     NUMBER 7 (IPCMO=7) OF RATE PHASE-CHANGE FORMULATIONS
C
C     Notes:
C
C     IMODE: number of thermal-mechanical-microscopical properties
C            defined in propmic7.f and idepros7.f
C
C     IMECHMIC: it must be the last property defined in idepros7.f
C
C     The following variables are defined in idenpr.f (ideipt.f, etc.)
C     VPLATM(IPLAT,1)=ILSPC
C     VPLATM(IPLAT,2)=ISSPC
C     VPLATM(IPLAT,3)=EXPAN
C     VPLATM(IPLAT,7)=ITEMRO
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nued_om.f'   ! thermal-microstructural
      INCLUDE 'nuee_om.f'   ! mechanical-microstructural
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'prob_om.f'
C
      DIMENSION PROPS(*)
C
      IA2=IA1+2     ! 2=ipcfo,ipcmo
C
      IVERSI=INT(PROPS(IA2+1))
C
      IF(IVERSI.EQ.1) THEN
C
C**** THERMAL-MECHANICAL-MICROSCOPICAL PROPERTIES (see idepros7.f)
C
      AM0=PROPS(IA2+2)
      AK=PROPS(IA2+3)
      VISC7=PROPS(IA2+4)
      UNIVC=PROPS(IA2+5)
      ENERA=PROPS(IA2+6)
      HENER=PROPS(IA2+7)
      CHININF=PROPS(IA2+8)
      CONTA=PROPS(IA2+9)
      CONTB=PROPS(IA2+10)
      IMECHMIC=INT(PROPS(IA2+11))
C
      VPLATM(IPLAT,8)=REAL(IVERSI)
      VPLATM(IPLAT,9)=AM0
      VPLATM(IPLAT,10)=AK
      VPLATM(IPLAT,11)=VISC7
      VPLATM(IPLAT,12)=UNIVC
      VPLATM(IPLAT,13)=ENERA
      VPLATM(IPLAT,14)=HENER
      VPLATM(IPLAT,15)=CHININF
      VPLATM(IPLAT,16)=CONTA
      VPLATM(IPLAT,17)=CONTB
      VPLATM(IPLAT,18)=REAL(IMECHMIC)
      IMODE=11      ! imode=total number of therm.-mic. prop. of model 7
      ENDIF         ! iversi.eq.1
C
      IF(IVERSI.EQ.2) THEN
C
C**** THERMAL-MECHANICAL-MICROSCOPICAL PROPERTIES (see idepros7.f)
C
      AKBAR=PROPS(IA2+2)
      ETA00=PROPS(IA2+3)
      ETABAR=PROPS(IA2+4)
      GHINF=PROPS(IA2+5)
      ALPHABAR=PROPS(IA2+6)
      AABAR=PROPS(IA2+7)
      UNIVC=PROPS(IA2+8)
      ENERA=PROPS(IA2+9)
      HENER=PROPS(IA2+10)
      CONTA=PROPS(IA2+11)
      CONTB=PROPS(IA2+12)
      ENETE=PROPS(IA2+13)
      FINFIN=PROPS(IA2+14)
      FINFINM=PROPS(IA2+15)
      GSETT=PROPS(IA2+16)
      TMAX=PROPS(IA2+17)
      TREF=PROPS(IA2+18)
      CONAE=PROPS(IA2+19)
      IMECHMIC=INT(PROPS(IA2+20))
C
      VPLATM(IPLAT,8)=REAL(IVERSI)
      VPLATM(IPLAT,9)=AKBAR
      VPLATM(IPLAT,10)=ETA00
      VPLATM(IPLAT,11)=ETABAR
      VPLATM(IPLAT,12)=GHINF
      VPLATM(IPLAT,13)=ALPHABAR
      VPLATM(IPLAT,14)=AABAR
      VPLATM(IPLAT,15)=UNIVC
      VPLATM(IPLAT,16)=ENERA
      VPLATM(IPLAT,17)=HENER
      VPLATM(IPLAT,18)=CONTA
      VPLATM(IPLAT,19)=CONTB
      VPLATM(IPLAT,20)=ENETE
      VPLATM(IPLAT,21)=FINFIN
      VPLATM(IPLAT,22)=FINFINM
      VPLATM(IPLAT,23)=GSETT
      VPLATM(IPLAT,24)=TMAX
      VPLATM(IPLAT,25)=TREF
      VPLATM(IPLAT,26)=CONAE
      VPLATM(IPLAT,27)=REAL(IMECHMIC)
      IMODE=20      ! imode=total number of therm.-mic. prop. of model 7
      ENDIF         ! iversi.eq.2 
C
C**** IMECHMIC=0 => mechanical model not influenced by microstructure
C     IMECHMIC=1 => mechanical model influenced by microstructure
C
      IF(IMECHMIC.EQ.1) THEN
C
C**** CORRECTS THE ASSIGMENT OF VPLATM GIVEN IN idenpr.f (ideipt.f, etc)
C
       I4A=IA1+1
       VPLATM(IPLAT,4)=PROPS(I4A)          ! IPCFOM
       I4B=I4A+1
       VPLATM(IPLAT,5)=PROPS(I4B)          ! IPCMOM
      ENDIF
C
      IA1=IA2+IMODE
C
      RETURN
      END
