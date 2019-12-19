      SUBROUTINE IDEPROM6(PROPS,IPLAT,IA1,INUPM)
C***********************************************************************
C
C**** THIS ROUTINE ORDERS THE MICROSTRUCTURAL PROPERTIES OF MODEL
C     NUMBER 6 (IPCMO=6) OF RATE PHASE-CHANGE FORMULATIONS
C
C***********************************************************************
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
      IVERSI=INT(PROPS(IA2+1))
      
c*****************************************************************
c     
c     VERSION 1 and 2
c     
c*****************************************************************
      if ((IVERSI.eq.1).or.(IVERSI.eq.2)) then
         IA3=6                  ! 6=thermal-microscopical properties
         IF(IMICRM.GT.0) THEN
c     mechanical-microscopical properties
c            VPLATM(IPLAT,1)=PROPS(IA2+IA3+1) ! ILSPC
c            VPLATM(IPLAT,2)=PROPS(IA2+IA3+2) ! ISSPC
c            VPLATM(IPLAT,3)=PROPS(IA2+IA3+3) ! EXPAN
         ENDIF
         IMODE=0                ! imode=total number of mech. prop. of model 6
         IA1=IA2+IA3+IMODE
      endif

c*****************************************************************
c     
c     VERSION 3 (Austempered Ductile Iron model)
c     
c*****************************************************************
      if (IVERSI.eq.3) then

C--   OF PROPS TO PARAMETERS
         NCOPC=INT(IPLAC(INUPM,1,1))
         if (NCOPC.eq.0) then
            nmm=INT(PROPS(IA2+16))
            TBF=(PROPS(IA2+12))
            alphaf=(PROPS(IA2+17))
            alphaa=(PROPS(IA2+18))
         endif
         if (NCOPC.eq.1) then
            nmm=INT(PROPS(IA2+7))
            TBF=(PROPS(IA2+3))
            alphaf=(PROPS(IA2+8))
            alphaa=(PROPS(IA2+9))
         endif
      endif

c--   OF PARAMETERS TO VPLATM
      if(IVERSI.eq.3) then
         VPLATM(IPLAT,4)=1      !IPCFOM (gt 0)
         VPLATM(IPLAT,5)=6      !IPCMOM (micro model number)
         VPLATM(IPLAT,8)=IVERSI
         VPLATM(IPLAT,9)=nmm

cccccccccccccccccccccccccccccccccccccccccccccc
c--   mechanical model 1 (bhadeshia's model)
cccccccccccccccccccccccccccccccccccccccccccccc
         if (nmm.eq.1) then
            VPLATM(IPLAT,10)=TBF
            VPLATM(IPLAT,11)=alphaf
            VPLATM(IPLAT,12)=alphaa
            if (NCOPC.eq.0) then
               IMODE=18
            endif
            if (NCOPC.eq.1) then
               IMODE=9
            endif
         endif
      endif
c     
      IA1=IA2+IMODE
      
c      open(UNIT=16,FILE='ideprom6.txt',STATUS='UNKNOWN')
c      write(16,*) 'ideprom6:',INUPM,IVERSI,nmm,TBF,alphaf,alphaa
      
      RETURN
      END
