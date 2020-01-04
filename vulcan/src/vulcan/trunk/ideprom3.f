      SUBROUTINE IDEPROM3(PROPS,IPLAT,IA1)
C***********************************************************************
C
C**** THIS ROUTINE ORDERS THE MICROSTRUCTURAL PROPERTIES OF MODEL
C     NUMBER 3 (IPCMO=3) OF RATE PHASE-CHANGE FORMULATIONS
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
      IA2=IA1+2         ! 2=ipcfo,ipcmo
      IA3=3             ! 3=thermal-microscopical properties
C
      IF(IMICRM.GT.0) THEN

c     mechanical-microscopical properties

c      VPLATM(IPLAT,1)=PROPS(IA2+IA3+1)    ! ILSPC
c      VPLATM(IPLAT,2)=PROPS(IA2+IA3+2)    ! ISSPC
c      VPLATM(IPLAT,3)=PROPS(IA2+IA3+3)    ! EXPAN

      ENDIF
C
      IMODE=0             ! imode=total number of mech. prop. of model 1
C
      IA1=IA2+IA3+IMODE
C
      RETURN
      END