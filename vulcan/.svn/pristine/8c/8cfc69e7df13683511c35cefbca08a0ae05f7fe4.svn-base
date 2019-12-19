      SUBROUTINE IDEPROM14(PROPS,IPLAT,IA1)
c*******************************************************************
c     THIS ROUTINE IDENTIFICATE IDENTIFIES THE MICROSTRUCTURAL
c     MATERIAL PROPERTIES FOR MODEL 14 OF RATE PHASE-CHANGE
C     FORMULATION
c*******************************************************************
c
c
c
c
c
c
c*******************************************************************
c     autor: adrian boccardo
c     date: last write 03-26-2014
c*******************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
     
C--   COUPLING VARIABLES
      INCLUDE 'nued_om.f'       ! thermal-microstructural
      INCLUDE 'nuee_om.f'       ! mechanical-microstructural

C--   MECHANICAL VARIABLES
      INCLUDE 'auxl_om.f'
      INCLUDE 'prob_om.f'

      DIMENSION PROPS(*)
      
      IA2=IA1+2                 ! 2=ipcfo,ipcmo

C--   OF PROPS TO PARAMETERS
      nmm=INT(PROPS(IA2+17))
      
ccccccccccccccccccccccccccccccccccccccccccccc
c--   mechanical model 1 (constant expansion)
ccccccccccccccccccccccccccccccccccccccccccccc
      if (nmm.eq.1) then

c--   of props to parameter of mechanical model 1
         svc=PROPS(IA2+18)

c--   of parameter of model 1 to vplatm
         VPLATM(IPLAT,4)=1      !IPCFOM (gt 0)
         VPLATM(IPLAT,5)=14     !IPCMOM (micro model number)
         VPLATM(IPLAT,8)=nmm
         VPLATM(IPLAT,9)=svc

         IMODE=18
c         IMODE=2                !num of parameters of mechanical model 14
      endif

      IA1=IA2+IMODE

c      open(UNIT=11,FILE='salida.txt',STATUS='UNKNOWN')
c      write(11,*) 'ideprom14',svc

      RETURN
      END
