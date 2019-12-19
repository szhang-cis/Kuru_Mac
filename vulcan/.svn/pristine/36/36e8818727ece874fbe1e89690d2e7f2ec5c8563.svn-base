      SUBROUTINE MICROS14(TGAUAT,TGAUST,TGAUIT,DTEMPT,
     .     BASMM, BASCC, BASKK,
     .     TSOE2, TSOE1, TSOC1,
     .     IPLAT,
     .     ALPHAM,INUPC)  
c*****************************************************************
c 
c     MODEL 14 (Boccardo's model for inverse eutectoid
c     transformation)
c     
c*****************************************************************
c
c     ALPHAM(IN+1)= function of phase change (fphch).
c     ALPHAM(IN+2)= volume fraction of graphite (vfGr).
c     ALPHAM(IN+3)= volume fraction of ferrite (stable transformation) (vfFs).
c     ALPHAM(IN+4)= volume fraction of austenite (stable transformation) (vfAs).
c     ALPHAM(IN+5)= volume fraction of pearlite (vfP).
c     ALPHAM(IN+6)= volume fraction of cementite (vfCm).
c     ALPHAM(IN+7)= volume fraction of ferrite (metastable transformation) (vfFm).
c     ALPHAM(IN+8)= volume fraction of austenite (metastable transformation) (vfAm).
c     ALPHAM(IN+9)= volume fraction of austenite (stable and metastable transformation) (vfA).
c     ALPHAM(IN+10)= graphite radius [m] (RGr).
c     ALPHAM(IN+11)= austenite radius [m] (RA).
c     ALPHAM(IN+12)= average carbon content in austenite [W%] (cA).
c     ALPHAM(IN+13)= position of cementite-austenite interface [m] (XCm).
c     ALPHAM(IN+14)= position of austenite-ferrite interface [m] (XA).
c
c*****************************************************************
c     autor: adrian boccardo.
c     date: last write 06-01-2014 by adrian boccardo.
c*****************************************************************
c
      IMPLICIT REAL*8 (A-H,O-Z)
C     
C--   COUPLING VARIABLES    ! thermal-microstructural
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
      DIMENSION BASKK(NDIMETO)
C     
      DIMENSION TSOE1(5), TSOE2(5), TSOC1(5)
c
C--   TRANSFERS "VPLAT" TO MICROSCOPICAL PARAMETERS
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
      nsds=VPLAT(IPLAT,20)
      IFPCDT=VPLAT(IPLAT,21)
c
      IV=21                     !number of VPLAT defined in input data

C     
C--   TRANSFERS "ALPHAM" TO MICROSCOPICAL VARIABLES
C
      IN=INUPC
      fphch=ALPHAM(IN+1)
      vfGr=ALPHAM(IN+2)
      vfFs=ALPHAM(IN+3)
      vfAs=ALPHAM(IN+4)
      vfP=ALPHAM(IN+5)
      vfCm=ALPHAM(IN+6)
      vfFm=ALPHAM(IN+7)
      vfAm=ALPHAM(IN+8)
      vfA=ALPHAM(IN+9)
      RGr=ALPHAM(IN+10)
      RA=ALPHAM(IN+11)
      cA=ALPHAM(IN+12)
      XCm=ALPHAM(IN+13)
      XA=ALPHAM(IN+14)

C     
C--   SOLVE MODEL 14
C
     
c--   some declarations
c     for linear density of graphite
      TempidGr=20.D0
      TempfdGr=1000.D0
      alphagr=9.26D-6           !alpha graphite
      cGr=100.D0                !carbon concentration graphite
      cCm=6.67D0                !carbon concentration cementite
      denCm=7660.D0             !cementite
      denGret=2200.D0           !graphite to enviroment temperature

c     form nodule/mm3 to nodule/m3
      grnden=1.D9*grnden       

c     index phase change.
c     if ipc=0 don't input to perlite or ferrite phase transformation.
c     if ipc=1 input to ferrite phase transformation.
c     if ipc=2 input to perlite phase transformation.
      ipc=0                    
      vfAslast=vfAs
      vfAmlast=vfAm
      TSOE1(IPLAT)=0.D0         !f_pc at time t+dt
      TSOE2(IPLAT)=0.D0         !f_pc at time t (last step)
      
c--   temperatures
      TempFo=739.D0+31.5D0*csi-7.7D0*ccu-18.7D0*cmn+3.3D0*cmo- ![alloying wi/wt%] Ghergu and Lacaze
     .     10.7D0*ccr-26.D0*cni
      TempF=739.D0+18.4D0*csi+2.D0*csi**2.D0-14.D0*ccu-
     .     45.D0*cmn+2.D0*cmo-24.D0*ccr-27.5D0*cni
      TempP=727.D0+21.6D0*csi+0.023D0*csi**2.D0+8.D0*cmo+13.D0*ccr- ![alloying wi/w%t] Ghergu and Lacaze
     .     21.D0*ccu-25.D0*cmn-33.D0*cni
      TempPo=727.D0+30.07D0*csi-1.98D0*csi**2.D0-10.7D0*ccu-13.7D0*cmn+ ![alloying wi/wt%] Ghergu and Lacaze
     .     9.3D0*cmo+24.3D0*ccr-12.D0*cni

c--   carbon concentration
      CP=(0.1876D0-4.112D-4*TempPo+2.26D-7*TempPo**2.D0+ !cp equal to cAF for Temp 
     .     0.125D0*csi/100.D0)*100.D0 !equal to TempPo, look at the isopleth from
                                      ! Ghergu and Laceze for understand Temp
                                      ! equal to TempPo and not Temp
                                      !equal to TempP.
      cFGr=0.03D0
      cFCm=0.03D0
      cFA=(6.7D-4-5D-10*TGAUST**2.D0-2.8D-7*TGAUST+ !lacaze gerval [alloying wi/wt] [cFA wi/wt%]
     .     (1.2D-8*TGAUST**2.D0-4.5D-3)*csi/100.D0)*100.D0
      cAF=(0.1876D0-4.112D-4*TGAUST+2.26D-7*TGAUST**2.D0+
     .     0.125D0*csi/100.D0)*100.D0 !lacaze gerval [alloying wi/wt] [cAF wi/wt%]
      cAGr=(TGAUST-1154.6D0-6.5D0*csi)*(1.5D0-0.216D0*csi)/
     .     (354.6D0+6.5D0*csi)+2.1D0-0.216D0*csi  !boeri [all wi/wt%]
      cACm=(TGAUST-491.27D0)/365.95D0

c--   density
      denA=8099.79D0-0.506D0*TGAUST+(-118.26D0+0.00739D0*TGAUST)*cAGr
     .     -6.01D0*cmn          !I use it because cFo and cFop aren't stable.
                                !I can change it.
      denF=7875.96D0-0.297D0*TGAUST-5.62D-5*TGAUST**2.D0+
     .     (-206.35D0+0.00778D0*TGAUST+1.472D-6*TGAUST**2.D0)*cFGr-
     .     7.24D0*cmn
      denP=0.9D0*denF+0.1D0*denCm !aprox, where fvFp=0.9 and fvCmp=0.1
      denGrf=denGret/(3.D0*alphagr*(TempfdGr-TempidGr)+1.D0)
      denGr=(TGAUST-TempfdGr)/(TempidGr-TempfdGr)*denGret+
     .     (TGAUST-TempidGr)/(TempfdGr-TempidGr)*denGrf
      denT=vfGr*denGr+vfFs*denF+vfP*denP+vfA*denA

c--   radius
      pi=4.*atan(1.D0)
      Rcel=(3.D0/(4.D0*pi*grnden))**(1.D0/3.D0) ![m]
      RF=(3.D0*(1-vfPi)/(4.D0*pi*grnden))**(1.D0/3.D0) ![m]
      RAi=1.01D0*RGr            !increment ra [m]
      
c--   length
      vfcelP=vfPi               !vfAp+vfP
      XF=sip/2.D0               ![m]
      XAi=1.01D0*XCm            !position of xa [m]

c--   diffusion coeficient of carbon in austenite and ferrite
c      DfA=2.343D-5*exp(-17767.D0/(TGAUST+273.D0)) !Liu and Agren [m2/s]
      DfA=1.67D-6*exp(-1.2D5/(8.314472*(TGAUST+273.D0))) !Kapturkiewicz et al [m2/s]

      Tc=1043.D0-1000.D0*csi/100.D0 ![alloying wi/wt]
c      DfF=2.D-6*exp(-10115.D0/(TGAUST+273.D0))* !Agren [m2/s]
c     .     exp(0.5898D0*(1.D0+2.D0/pi*atan(15629.D0/Tc-
c     .     15309.D0/(TGAUST+273.D0))))
      DfF=7.9D-7*exp(-7.58D4/(8.314472*(TGAUST+273.D0))) !Kapturkiewicz et al [m2/s]
      
cccccccccccccccccccccccccccccccccccccccccccccccccc
cc----------- ferrite homogenization -----------cc
cccccccccccccccccccccccccccccccccccccccccccccccccc
      if (vfAs.eq.0.D0.and.vfFs.gt.0.D0) then 
c--   cbm(cFo)
         denT1=(denGr*vfGr+denF*vfFs)/(1.D0-vfcelP) !for upload in the loop
         denT=vfGr*denGr+vfFs*denF+vfP*denP+vfA*denA !for upload in the loop
         CT1=(cc*denT-CP*denP*vfcelP)/(denT1*(1.D0-vfcelP))
         cbm=(CT1*denT1-vfGr/(1.D0-vfcelP)*cGr*denGr)/
     .        (vfFs/(1.D0-vfcelP)*denF)

c     controls: very important because there are problems if the solver doesn't
c     use it. 1) Is very difficult that ferrite doesn't have carbon. 2) Is very
c     difficult that ferrite has more than two time of carbon concentration
c     (I suppose). 
         if (cbm.lt.0.01D0) cbm=0.01D0
         if (cbm.gt.0.06D0) cbm=0.06D0 

c--   delta graphite
         delRGr=DfF*denF*RF*(cbm-cFGr)/(RGr*(RF-RGr)*
     .        (cGr*denGr-cFGr*denF))*DTIMET

c--   radius graphite at time t+delt
         RGr=RGr+delRGr

c--   volume fraction of graphite, austenite and ferrite at time t+delt
         vfGr=4.D0/3.D0*pi*grnden*RGr**3.D0
         vfFs=4.D0/3.D0*pi*grnden*(RF**3.D0-RGr**3.D0)
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccc
cc----------- phase change ferrite -------------cc
cccccccccccccccccccccccccccccccccccccccccccccccccc
      if ((TGAUST.ge.TempFo.or.vfAs.gt.0.D0).and.vfFs.gt.0.D0) then
c--   ipc=1
         ipc=1
         
c--   RA for frist input
         if (RA.le.RGr) RA=RAi

c--   density
         denT1=(denGr*vfGr+denF*vfFs+denA*vfAs)/(1.D0-vfcelP) !for upload in the loop
         denT=vfGr*denGr+vfFs*denF+vfP*denP+vfA*denA !for upload in the loop
         CT1=(cc*denT-CP*denP*vfcelP)/(denT1*(1.D0-vfcelP))

c--   cbm(cFo)
         cbm=(CT1*denT1-vfGr/(1.D0-vfcelP)*cGr*denGr-(cAGr+cAF)/2.D0*
     .        vfAs/(1.D0-vfcelP)*denA)/(vfFs/(1.D0-vfcelP)*denF)

c     controls: very important because there are problems if the solver doesn't
c     use it. 1) Is very difficult that ferrite doesn't have carbon. 2) Is very
c     difficult that ferrite has more than two time of carbon concentration
c     (I suppose). 
         if (cbm.lt.0.01D0) cbm=0.01D0
         if (cbm.gt.0.06D0) cbm=0.06D0 

c--   delta radius of graphite
         delRGr=DfA*denA*RA*(cAF-cAGr)/(RGr*(RA-RGr)*
     .        (cGr*denGr-cAGr*denA))*DTIMET
         if (delRGr.gt.0.D0) delRGr=0.D0 ! I use it for stabilize the solve

c--   delta radius of austenite
         delRA=(DfF*denF*RF*(cbm-cFA)/(RA*(RF-RA))-DfA*denA*
     .        RGr*(cAF-cAGr)/(RA*(RA-RGr)))/
     .        (cAF*denA-cFA*denF)*DTIMET
         if (delRA.lt.0.D0) delRA=0.D0 ! I use it for stabilize the solve

c--   radius of graphite at time t+delt
         RGr=RGr+delRGr

c--   radius of austenite at time t+delt
         RA=RA+delRA

c--   volumetrix fraction of graphite, austenite and ferrite at time t+delt
         vfGr=4.D0/3.D0*pi*grnden*RGr**3.D0
         vfAs=4.D0/3.D0*pi*grnden*(RA**3.D0-RGr**3.D0)
         vfFs=4.D0/3.D0*pi*grnden*(RF**3.D0-RA**3.D0)

c--   for solve some problems
         if ((RF-RA)/RF.le.0.01D0) then
            vfAs=4.D0/3.D0*pi*grnden*(RF**3.D0-RGr**3.D0)
            vfFs=0.D0
         endif
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccc
cc----------- phase change pearlite ------------cc
cccccccccccccccccccccccccccccccccccccccccccccccccc
      if (TempPo.gt.TempP) then
         TempTP=TempPo
      else
         TempTP=TempP
      endif
      if ((TGAUST.ge.TempTP.or.vfAm.gt.0.D0).and.(vfP.gt.0.D0)) then
c--   ipc=2
         ipc=2

c--   XA for frist input
         if (XA.lt.XAi) XA=XAi

c--   cbm(cFo)
         cbm=2.D0*(CP*denP-cCm*denCm*vfCm/vfcelP-
     .        (cACm+cAF)/2.D0*denA*vfAm/vfcelP-
     .        cFA/2.D0*denF*(1.D0-(vfCm+vfAm)/vfcelP))/
     .        (denF*(1.D0-(vfCm+vfAm)/vfcelP))

c     controls: very important because there are problems if the solver doesn't
c     use it. 1) Is very difficult that ferrite doesn't have carbon. 2) Is very
c     difficult that ferrite has more than two time of carbon concentration
c     (I suppose). 
         if (cbm.lt.0.01D0) cbm=0.01D0
         if (cbm.gt.0.06D0) cbm=0.06D0 
         
c--   delta cementite
         delXCm=DfA*denA*(cAF-cACm)*DTIMET/((XA-XCm)*
     .        (cCm*denCm-cACm*denA))
         if (delXCm.gt.0.D0) delXCm=0.D0 !for help to solve

c--   delta austenite
         delXA=-DfA*denA*(cAF-cACm)*DTIMET/((XA-XCm)*
     .        (cAF*denA-cFA*denF))+DfF*denF*(cbm-cFA)/
     .        ((XF-XA)*(cAF*denA-cFA*denF))*DTIMET
         if (delXA.lt.0.D0) delXA=0.D0 !for stabilize the solve
         
c     this control and change the growth of XCm and XA. It is very importand
c     when perlite transformation begins
         vfPc=1.D0-vfP/vfcelP
         if (vfPc.le.0.2D0) then
            delXCmmax=-0.05D0*XCm !proposed
            if (delXCm.lt.delXCmmax) delXCm=delXCmmax !for stabilize the solve
            delXAmax=0.05D0*(XF-XA) !proposed
            if (delXA.gt.delXAmax) delXA=delXAmax !for stabilize the solve
         endif
         if ((vfPc.gt.0.2D0).and.(vfPc.le.0.4D0)) then
            delXCmmax=-0.1D0*XCm !proposed
            if (delXCm.lt.delXCmmax) delXCm=delXCmmax !for stabilize the solve
            delXAmax=0.1D0*(XF-XA) !proposed
            if (delXA.gt.delXAmax) delXA=delXAmax !for stabilize the solve
         endif
         
c--   length cementite at time t+delt
         XCm=XCm+delXCm
         if (XCm.lt.0.D0) XCm=0.D0

c--   length austenite at time t+delt
         XA=XA+delXA
         if (XA.gt.XF) XA=XF

c--   for help solve
         if ((XF-XA)/XF.le.0.01D0) then
            XCm=0.D0
            XA=XF
         endif

c--   volumetrix fraction of cementite, austenite, ferrite and perlite at
c     time t+delt
         vfP=vfcelP-vfAm
         vfCm=(XCm/XF)*vfcelP
         vfFm=vfP-vfCm
         vfAm=((XA-XCm)/XF)*vfcelP
      endif
      
cccccccccccccccccccccccccccccccccccccccccccccccccc
cc-------- volume fraction of austenite --------cc
cccccccccccccccccccccccccccccccccccccccccccccccccc
      if ((ipc.eq.1).or.(ipc.eq.2)) then
         vfA=vfAs+vfAm
         fphch=vfA/(1.D0-vfGr)  !look
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccc
cc---------- austenite homogenization ----------cc
cccccccccccccccccccccccccccccccccccccccccccccccccc
c--   TempAH is the temperatute for begin to homogenization. The 
c     subroutine has to use it because in upon cooling, when it's
c     couple with austempered subroutine, austenite loses mass
c     carbon if temp is smaller than TempAH. it's not possible 
c     because austenite transfor in perlite and ferrite.
      if (TempP.gt.TempF) then
         TempAH=TempP
      else
         TempAH=TempF
      endif
      if (TGAUST.ge.TempAH) then
         if (vfP.eq.0.D0.and.vfFs.eq.0.D0) then
            
c--   cbm(cAo)
            denT=vfGr*denGr+vfA*denA
            cbm=(cc*denT-vfGr*cGr*denGr)/(vfA*denA)
            
c     controls: very important because there are problems if the solver doesn't
c     use it. 1) Is very difficult that austenite doesn't have carbon. 2) Is very
c     difficult that austenite has more than two (more or less) time of carbon
c     concentration (I suppose). 
            if (cbm.lt.0.3D0) cbm=0.3D0
            if (cbm.gt.1.5D0) cbm=1.5D0 
            
c--   delta graphite
            delRGr=DfA*denA*Rcel*(cbm-cAGr)/(RGr*(Rcel-RGr)*
     .           (cGr*denGr-cAGr*denA))*DTIMET

c--   radius graphite at time t+delt
            RGr=RGr+delRGr

c--   volumetrix fraction of graphite, ferrite and austenite at time t+delt
            vfGr=4.D0/3.D0*pi*grnden*RGr**3.D0
            vfA=1.D0-vfGr

c     cA=(cAGr+cbm)/2.D0       !propose
            cA=cbm              !propose
c            cA=cAGr             !propose
         endif
      endif

c     
C**** OPTIONS IN THE EVALUATION OF df_pc/dT !influye en el jacobiano
C
C     1) standard df_pc/dT = delta f_pc / delta T (Diego's thesis) (con recalecencia da negativo y rope algoritmo)
C     2) standard df_pc/dT = delta f_pc / delta T (always ge 0)
C     3) df_pc/dT "exact" (not implemented)
C     4) df_pc/dT "exact" (always ge 0) (not implemented)
C     5) df_pc/dT = 0
C     
C     Notes:
C     Only one value of ICACOT can be defined, i.e., ICACOT does not
C     depend of any particular phase-change. This can not be useful for
C     more than one microstructural phase-changes.
C
      TSOE1(IPLAT)=vfAm*henerp+vfAs*henerf !L*f_pc at time t+dt
      TSOE2(IPLAT)=vfAmlast*henerp+vfAslast*henerf !L*f_pc at time t

      IF(IFPCDT.NE.3.AND.ICONVT.EQ.1)
     .     CALL RUNENDT('ERROR: IFPCDT NE 3 WHEN ICONVT=1')
C     
      GO TO (1,2,3,4,5) IFPCDT
C     
    1 ICACOT=0
      VELTOT=1.D-10
      VELA1T=DABS(DTEMPT)
      IF(VELA1T.GT.VELTOT) THEN
         IF(IITERT.GT.0)
     .        TSOC1(IPLAT)=(TSOE1(IPLAT)-TSOE2(IPLAT))/(DTEMPT*DTIMET)
      ENDIF
      GO TO 10
C     
    2 ICACOT=1
      VELTOT=1.D-10
      VELA1T=DABS(DTEMPT)
      IF(VELA1T.GT.VELTOT) THEN
         IF(IITERT.GT.0)
     .        TSOC1(IPLAT)=(TSOE1(IPLAT)-TSOE2(IPLAT))/(DTEMPT*DTIMET)
      ENDIF
      GO TO 10
C     
    3 ICACOT=0
      TSOC1(IPLAT)=TSOCGN+TSOCGG+TSOCGA
      TSOC1(IPLAT)=TSOC1(IPLAT)*HENER
      GO TO 10
C     
    4 ICACOT=1
      TSOC1(IPLAT)=TSOCGN+TSOCGG+TSOCGA
      TSOC1(IPLAT)=TSOC1(IPLAT)*HENER
      GO TO 10
C     
    5 ICACOT=1
      TSOC1(IPLAT)=0.0D0
      GO TO 10
C     
 10   CONTINUE

C     
C--   TRANSFERS MICROSCOPICAL VARIABLES TO "ALPHAM"
C
      ALPHAM(IN+1)=fphch
      ALPHAM(IN+2)=vfGr
      ALPHAM(IN+3)=vfFs
      ALPHAM(IN+4)=vfAs
      ALPHAM(IN+5)=vfP
      ALPHAM(IN+6)=vfCm
      ALPHAM(IN+7)=vfFm
      ALPHAM(IN+8)=vfAm
      ALPHAM(IN+9)=vfA
      ALPHAM(IN+10)=RGr
      ALPHAM(IN+11)=RA
      ALPHAM(IN+12)=cA
      ALPHAM(IN+13)=XCm
      ALPHAM(IN+14)=XA
      
      INUPC=INUPC+14            !increments ALPHAM index. less than NBASES; 
                                !see pointes.f
      
      RETURN
      END
