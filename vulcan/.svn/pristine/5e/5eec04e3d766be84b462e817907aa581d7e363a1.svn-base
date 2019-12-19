      SUBROUTINE INMICRO14(ALPHAM,INUPC,IPLAT)
C***********************************************************************
C     
C**** THIS ROUTINE INITIALITES MICROSTRUCTURAL PARAMETERS FOR MODEL 14
C     
C***********************************************************************
C     
C     Index of variables
C     
C     ALPHAM= array of microstructural (microscopical) variables
C     
C***********************************************************************
C     autor: adrian boccardo
C     date: last write 06-01-2014
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C     
C--   COUPLING VARIABLES (thermal-microstructural)
C     
      INCLUDE 'nued_om.f'
C     
C--   THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C     
      DIMENSION ALPHAM(NHISTM)
C     
C--   TRANSFERS "VPLAT" TO MICROSCOPICAL VARIABLES
C     
      cc=VPLAT(IPLAT,6)
      csi=VPLAT(IPLAT,7)
      cmn=VPLAT(IPLAT,8)
      ccr=VPLAT(IPLAT,9)
      cni=VPLAT(IPLAT,10)
      ccu=VPLAT(IPLAT,11)
      cmo=VPLAT(IPLAT,12)
      vfGri=VPLAT(IPLAT,13)
      vfFi=VPLAT(IPLAT,14)
      vfPi=VPLAT(IPLAT,15)
      sip=VPLAT(IPLAT,16)
      grnden=VPLAT(IPLAT,17)
      henerf=VPLAT(IPLAT,18)
      henerp=VPLAT(IPLAT,19)
      sds=VPLAT(IPLAT,20)
      IFPCDT=VPLAT(IPLAT,21)
c     
      IV=21                     !number of VPLAT defined in input data
C     
C--   INITIALITES MICROSTRUCTURAL VARIABLES ("ALPHAM")
C     see micros14
      IN=INUPC
      ALPHAM(IN+1)=0.D0         !fphch
      ALPHAM(IN+2)=vfGri        !vfGr
      ALPHAM(IN+3)=vfFi         !vfFs
      ALPHAM(IN+4)=0.D0         !vfAs
      ALPHAM(IN+5)=vfPi         !vfP
c     
      cFCm=0.03D0 
      cCm=6.67D0                !carbon concentration cementite
      TempPo=727.D0+30.07D0*csi-1.98D0*csi**2.D0-10.7D0*ccu-13.7D0*cmn+ ![alloying wi/wt%]
     .     9.3D0*cmo+24.3D0*ccr-12.D0*cni ! Ghergu and Lacaze
      cP=(0.1876D0-4.112D-4*TempPo+2.26D-7*TempPo**2.D0+ !look micros14
     .     0.125D0*csi/100.D0)*100.D0
      vfCmi=vfPi*(cP-cFCm)/(cCm-cFCm)
      ALPHAM(IN+6)=vfCmi        !vfCm
      ALPHAM(IN+7)=vfPi-vfCmi   !vfFm
      ALPHAM(IN+8)=0.D0         !vfAm
      ALPHAM(IN+9)=0.D0         !vfA
      pi=4.D0*atan(1.D0)
      grnden=1.D9*grnden
      RGri=(3.D0*vfGri/(4.D0*pi*grnden))**(1.D0/3.D0) ![m]
      ALPHAM(IN+10)=RGri        !RGr [m]
      ALPHAM(IN+11)=0.D0        !RA [m]
      ALPHAM(IN+12)=0.D0        !cA [m]
      if (vfPi.eq.0.D0) then
         XCmi=0.D0
      else
         XCmi=(sip/2.D0)*(vfCmi/vfPi)
      endif
      ALPHAM(IN+13)=XCmi        !XCm [m]
      ALPHAM(IN+14)=0.D0        !XA [m]
c     
      INUPC=INUPC+14            !less than NBASES; see pointes.f. 
                                !se suma el numero de variables alpham
C
      RETURN
      END
