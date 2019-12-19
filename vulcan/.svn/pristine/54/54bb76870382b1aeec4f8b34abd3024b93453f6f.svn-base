      SUBROUTINE IDEPROS1(PROPST,IPLAT,IA1)
C***********************************************************************
C
C**** THIS ROUTINE ORDERS THE MICROSTRUCTURAL PROPERTIES OF MODEL
C     NUMBER 1 (IPCMO=1) OF RATE PHASE-CHANGE FORMULATIONS
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** THERMAL VARIABLES
C
      INCLUDE 'auxl_omt.f'
      INCLUDE 'prob_omt.f'
C
      DIMENSION PROPST(*)
C
      IA2=IA1+2         ! 2=ipcfo,ipcmo
C
      VPLAT(IPLAT,1)=PROPST(IA2+1)   ! Ts
      VPLAT(IPLAT,2)=PROPST(IA2+2)   ! Tl
      VPLAT(IPLAT,3)=PROPST(IA2+3)   ! L
      VPLAT(IPLAT,6)=PROPST(IA2+4)   ! 0 for this model (see ideprom1.f)
C
      IMODE=4                   ! imode=total number of prop. of model 1
C
      IA1=IA2+IMODE
C
      RETURN
      END
