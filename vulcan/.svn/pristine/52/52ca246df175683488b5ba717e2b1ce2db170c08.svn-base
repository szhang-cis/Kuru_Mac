      SUBROUTINE IDEPROM1(PROPS,IPLAT,IA1)
C***********************************************************************
C
C**** THIS ROUTINE ORDERS THE MICROSTRUCTURAL PROPERTIES OF MODEL
C     NUMBER 1 (IPCMO=1) OF RATE PHASE-CHANGE FORMULATIONS
C
C     Notes:
C
C     IMODE: number of thermal-mechanical-microscopical properties
C            defined in propmic1.f and idepros1.f
C
C     IMECHMIC: it must be the last property defined in idepros1.f
C               (it does not play any role in this model, i.e., =0)
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
C**** THERMAL-MECHANICAL-MICROSCOPICAL PROPERTIES (see idepros1.f)
C
      VPLATM(IPLAT, 8)=PROPS(IA2+1)   ! Ts
      VPLATM(IPLAT, 9)=PROPS(IA2+2)   ! Tl
      VPLATM(IPLAT,10)=PROPS(IA2+3)   ! L
      VPLATM(IPLAT,11)=PROPS(IA2+4)   ! mechanical-microstructural index
C
      IMODE=4       ! imode=total number of therm.-mic. prop. of model 1
C
C**** IMECHMIC=0 => mechanical model not influenced by microstructure
C     IMECHMIC=1 => mechanical model influenced by microstructure
C
      IMECHMIC=INT(VPLATM(IPLAT,11))
      IF(IMECHMIC.EQ.1) THEN               ! not useful in this model
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
