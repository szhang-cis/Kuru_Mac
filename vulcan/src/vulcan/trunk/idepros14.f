      SUBROUTINE IDEPROS14(PROPST,IPLAT,IA1)
c*******************************************************************
c     THIS ROUTINE IDENTIFICATE IDENTIFIES THE MICROSTRUCTURAL
c     MATERIAL PROPERTIES FOR MODEL 14
c*******************************************************************
c
c
c
c
c
c
c*******************************************************************
c     autor: adrian boccardo
c     date: last write 06-01-2014
c*******************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'auxl_omt.f'
      INCLUDE 'prob_omt.f'
      DIMENSION PROPST(*)
      IA2=IA1+2         ! 2=ipcfo,ipcmo
C
C     OF PROPST TO PARAMETERS
C
      cc=PROPST(IA2+1)
      csi=PROPST(IA2+2)
      cmn=PROPST(IA2+3)
      ccr=PROPST(IA2+4)
      cni=PROPST(IA2+5)
      ccu=PROPST(IA2+6)
      cmo=PROPST(IA2+7)
      vfGri=PROPST(IA2+8)
      vfFi=PROPST(IA2+9)
      vfPi=PROPST(IA2+10)
      sip=PROPST(IA2+11)
      grnden=PROPST(IA2+12)
      henerf=PROPST(IA2+13)
      henerp=PROPST(IA2+14)
      sds=PROPST(IA2+15)
      IFPCDT=PROPST(IA2+16)

c     
c     OF PARAMETERS TO VPLAT. begin of VPLAT(IPLAT,6)
c
      VPLAT(IPLAT,6)=cc
      VPLAT(IPLAT,7)=csi
      VPLAT(IPLAT,8)=cmn
      VPLAT(IPLAT,9)=ccr
      VPLAT(IPLAT,10)=cni
      VPLAT(IPLAT,11)=ccu
      VPLAT(IPLAT,12)=cmo
      VPLAT(IPLAT,13)=vfGri
      VPLAT(IPLAT,14)=vfFi
      VPLAT(IPLAT,15)=vfPi
      VPLAT(IPLAT,16)=sip
      VPLAT(IPLAT,17)=grnden
      VPLAT(IPLAT,18)=henerf
      VPLAT(IPLAT,19)=henerp
      VPLAT(IPLAT,20)=sds
      VPLAT(IPLAT,21)=IFPCDT
C     
      IMODE=18                  !num of parameters model 14
c     
      IA1=IA2+IMODE
      RETURN
      END
