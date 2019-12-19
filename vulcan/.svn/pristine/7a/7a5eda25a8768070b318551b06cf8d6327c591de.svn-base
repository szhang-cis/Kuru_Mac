      SUBROUTINE MICROS6(TGAUAT,TGAUST,TGAUIT,DTEMPT,
     .                    BASMM, BASCC, BASKK,
     .                    TSOE2, TSOE1, TSOC1,
     .                    IPLAT,INUPM,
     .                   ALPHAM,INUPC)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE PHASE-CHANGE FUNCTION, DENSITY,
C     CAPACITY & CONDUCTIVITY ACCORDING WITH THE MICROSTRUCTURAL MODEL
C     NUMBER 6 (IPCMO=6) OF RATE PHASE-CHANGE FORMULATIONS
C
C     S.G. CAST IRON MODEL: AVRAMI-JOHNSON-MEHL EQUATION
C
C***********************************************************************
C
C     Index of variables
C
C     TGAUAT= Temperature at time t
C     TGAUST= Temperature at time t+dt
C     TGAUIT= Initial temperature
C     DTEMPT= Temperature rate
C
C     BASMM = Density
C     BASCC = Capacity coefficient
C     BASKK = Conductivity
C
C     TSOE2 = L*Phase-change function at time t
C     TSOE1 = L*Phase-change function at time t+dt
C
C     ALPHAM= array of microstructural (microscopical) variables
C             First indexes of ALPHAM devoted to f_pc functions
C             (see outmic.f)
C
C     IN=INUPC (INUPC: number of microscopical phase-changes)
C     First indexes of ALPHAM devoted to f_pc functions
C
C     IVERSI=1
C     ALPHAM(IN+1): solid-solid phase-change function (at time t+dt)
C     ALPHAM(IN+2): solid-solid nucleation fraction (at time t+dt)
C     ALPHAM(IN+3): nucleation time
C     ALPHAM(IN+4): Td (temperature when nucleation time=1)
C
C     IVERSI=2
C     ALPHAM(IN+1): solid-solid phase-change function (at time t+dt)
C     ALPHAM(IN+2): ferrite fraction (at time t+dt)
C     ALPHAM(IN+3): density grains of ferrite (at time t+dt)
C     ALPHAM(IN+4): ferrite radius (at time t+dt)
C     ALPHAM(IN+5): pearlite fraction (at time t+dt)
C     ALPHAM(IN+6): density grains of pearlite (at time t+dt)
C     ALPHAM(IN+7): pearlite radius (at time t+dt)
C
C     IVERSI=3 (Boccardo's model. Last write 06-01-2014)
c     IN= position in ALPHAM array that take account the position
c     used for a previous model.
c     with microstructural couple model:
c     ALPHAM(IN+1)= volume fraction of bainitic ferrite (at time t+dt) (vfBF)
c     ALPHAM(IN+2)= volume fraction of residual austenite (at time t+dt) (vfA)
c     ALPHAM(IN+3)= average carbon content in austenite [W%] (cA).
c     ALPHAM(IN+4)= average carbon content in residual austenite [W%] (cAR).
c     ALPHAM(IN+5)= volume fraction of graphite (at time t+dt) (vfGr)
c     ALPHAM(IN+6)= time of bainitic transformation [min] (httime)
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

c     coupling variables (thermal-microstructural)
      INCLUDE 'nued_om.f'

c     thermal variables
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C     
      DIMENSION ALPHAM(NHISTM)
      DIMENSION BASKK(NDIMETO)
      DIMENSION TSOE1(5), TSOE2(5), TSOC1(5)

c     transfer VPLAT to microscopical variables
      IVERSI=INT(VPLAT(IPLAT,6))

c*****************************************************************
c 
c     VERSION 1 (SOLIDUS-SOLIDUS PHASE-CHANGE MODEL)
c     
c*****************************************************************
      IF(IVERSI.EQ.1) THEN
         TEMSS=VPLAT(IPLAT,7)
         AKCOO=VPLAT(IPLAT,8)
         AKC11=VPLAT(IPLAT,9)
         AKC22=VPLAT(IPLAT,10)
         ANC11=VPLAT(IPLAT,11)
         HENER=VPLAT(IPLAT,12)  ! constant
c     number of VPLAT defined in input data
         IV=12                
c     transfer ALPHAM to microscopical variables (load last)
c     converged values
         IN=INUPC
         FPCSO=ALPHAM(IN+1)
         FRACO=ALPHAM(IN+2)
         TTINO=ALPHAM(IN+3)
         TEMPD=ALPHAM(IN+4)
c     solve
         IF(TEMSS.GT.TGAUST) THEN
            DFRAC=(TEMSS-TGAUST)*(TEMSS-TGAUST)/(AKCOO*TGAUST)
            DFRAC=DFRAC*DTIMET
            DTINO=DTIMET
         ELSE
            DFRAC=0.0
            DTINO=0.0
         ENDIF
         FRACT=FRACO+DFRAC
         TTINU=TTINO
         IF(FRACT.LE.1.0) THEN
            TTINU=TTINO+DTINO
            IF(FRACT.GT.0.0) TEMPD=TGAUST
         ELSE
            IF(FRACO.LT.1.0) THEN
               DTIMEX=(FRACT-1.0)/(FRACT-FRACO)*DTIMET
               TTINU=TTINO+DTIMEX
               TEMPD=TGAUAT+DTEMPT*DTIMEX
            ENDIF
         ENDIF
         AUX00=-AKC22/TGAUST
         AUX0D=-AKC22/TEMPD
         AUX11=EXP(AUX00)*(TTINU**ANC11)
         AUX1D=EXP(AUX0D)*(TTINU**ANC11)
         AUX22=-AKC11*AUX11
         AUX2D=-AKC11*AUX1D
         AUX33=0.0
         IF(FRACT.GE.1.0) AUX33=1.0
         FPCSS=(1.0-EXP(AUX22-AUX2D))*AUX33
         FPCSS=1.0-FPCSS
         TSOE1(IPLAT)=FPCSS     ! f_pc at time t+dt
         TSOE2(IPLAT)=FPCSO     ! f_pc at time t
      ENDIF                     ! iversi.eq.1

c*****************************************************************
c 
c     VERSION 2
c     
c*****************************************************************
      IF(IVERSI.EQ.2) THEN
c     interface with microst.f
         DELT=DTIMET
         TEMP=TGAUST
         TP=TGAUAT
         TEMFF=VPLAT(IPLAT,7)
         TEMPP=VPLAT(IPLAT,8)
         CONUF=VPLAT(IPLAT,9)
         CONUP=VPLAT(IPLAT,10)
         COGRF=VPLAT(IPLAT,11)
         COGRP=VPLAT(IPLAT,12)
         HENER=VPLAT(IPLAT,13)
C     number of VPLAT defined in input data
         IV=13        
c     transfer ALPHAM to microscopical variables (load last)
c     converged values
         IN=INUPC
         FPCSO=ALPHAM(IN+1)
         FRAFO=ALPHAM(IN+2)
         DENFO=ALPHAM(IN+3)
         RADFO=ALPHAM(IN+4)
         FRAPO=ALPHAM(IN+5)
         DENPO=ALPHAM(IN+6)
         RADPO=ALPHAM(IN+7)
c     load last converged values
         FPCSS=FPCSO
         FRAFN=FRAFO
         DENFN=DENFO
         RADFN=RADFO
         FRAPN=FRAPO
         DENPN=DENPO
         RADPN=RADPO
c     establish variables of model
         DELFG=0.0
         DELFC=0.0
         TSOCG=0.0
         TSOCC=0.0
         MZ=1
         IF(FPCSS.LE.0.0) THEN
            MZ=100
         ELSE
            IF(TEMP.LT.TEMFF) MZ=2
            IF(TEMP.LT.TEMPP) MZ=3
         ENDIF
c     solve
         IF(MZ.EQ.2.OR.MZ.EQ.3) THEN
            FG=FRAFN
            FGP=FRAFO
            XNGN=DENFN
            RGM=RADFN
            FL=FPCSS
            FC=FRAPN
            FCP=FRAPO
            XNCN=DENPN
            RCM=RADPN
            A1=CONUF
            A2=CONUP
            DG=COGRF
            DC=COGRP
            TEG=TEMFF
            TEC=TEMPP
            CALL GRAFCARB(TEMP,DELT,FG,FGP,DELFG,XNGN,RGM,FL,FA,FC,MZ,
     .           FCP,DELFC,XNCN,RCM,
     .           TSOCG,TSOCC,
     .           A1,A2,DG,DC,TEG,TEC)
            FP=FPCSO
            HEAT=HENER
            FRAFN=FG
            DENFN=XNGN
            RADFN=RGM
            FPCSS=FL
            FRAPN=FC
            DENPN=XNCN
            RADPN=RCM
         ENDIF
         TSOE1(IPLAT)=FL        ! f_pc at time t+dt
         TSOE2(IPLAT)=FP        ! f_pc at time t
         HENER=HEAT
      ENDIF                     ! iversi.eq.2

c*****************************************************************
c 
c     VERSION 3 (Austempered Ductile Iron model using
c     Avrami equation)
c     
c*****************************************************************
      IF(IVERSI.EQ.3) THEN
cc--  number of model that provide initial conditions
         NCOPC=IPLAC(INUPM,1,1)

cc--  transfer VPLAT to microscopical variables
         C=VPLAT(IPLAT,7)
         SI=VPLAT(IPLAT,8)
         CU=VPLAT(IPLAT,9)
         AMN=VPLAT(IPLAT,10)
         ANI=VPLAT(IPLAT,11)
         AMO=VPLAT(IPLAT,12)
         CR=VPLAT(IPLAT,13)
         FAINI=VPLAT(IPLAT,14)
         TAUS=VPLAT(IPLAT,15)
         TBS=VPLAT(IPLAT,16)
         TBF=VPLAT(IPLAT,17)
         CCFB=VPLAT(IPLAT,18)
         AKAVRA=VPLAT(IPLAT,19)
         AMAVRA=VPLAT(IPLAT,20)

c     number of VPLAT defined in input data
         IV=20

cc--  transfers ALPHAM to microscopical variables
         IN=INUPC

c     not couple micro model
         if (NCOPC.eq.0) then
            vfBF=ALPHAM(IN+1)
            vfA=ALPHAM(IN+2)
            cA=ALPHAM(IN+3)
            cAR=ALPHAM(IN+4)
            httime=ALPHAM(IN+6)
         endif

c     couple micro model (with micros14)
         if (NCOPC.ne.0) then
            do ICOPC=1,NCOPC
               IX1=IPLAC(INUPM,ICOPC+1,3) !vfGr
               IX2=IPLAC(INUPM,ICOPC+1,4) !vfFs
               IX5=IPLAC(INUPM,ICOPC+1,7) !vfCm
               IX8=IPLAC(INUPM,ICOPC+1,10) !vfA
               IX10=IPLAC(INUPM,ICOPC+1,11) !cA
            enddo
            vfGr=ALPHAM(IX1)
            vfFs=ALPHAM(IX2)
            vfCm=ALPHAM(IX5)
            vfA=ALPHAM(IX8)
            cA=ALPHAM(IX10)
            vfBF=ALPHAM(IN+1)
            cAR=ALPHAM(IN+4)
            httime=ALPHAM(IN+6)
         endif

c     cAR initial is great equal to cA each time micros6 is called. Using this
c     there isn't problem when this variables go to mechanical problem. I think
c     the problem it is for smoothing funtion, because there are a big gradient
c     when cAR is initialized to zero and suddenly is equal to cA when vfBF is
c     great to zero. Is correct rewrite cAR equal to cA. It is too important only
c     for microstructural coupled model, because for input file parameter is solve
c     in inmicros6.
         if (cAR.le.cA) cAR=cA

cc--  solve
c     parameters if don't use couple micro model
         if (NCOPC.eq.0) then
            vfCm=0.D0
            vfFs=0.D0
            vfGr=1.D0-FAINI
         endif

c     if subrountine doesn't use this temperatute, when temperature
c     is lesser than Taustempered the model continue works. When the
c     temperature is lesser than TBF there isn't more transformation.
         if ((TGAUST.le.TBS).and.(TGAUST.ge.TBF).and.(vfFs.eq.0.D0).and.
     .        (vfCm.eq.0.D0)) then
c     austenite carbon concentration if don't use couple micro model
            if (NCOPC.eq.0) cA=-0.435D0+3.35D-4*TAUS+
     .           1.61D-6*TAUS*TAUS+0.006D0*AMN-0.11D0*SI-0.07D0*ANI+
     .           0.014D0*CU-0.3D0*AMO
            
            cmAR=3.072D0-0.0016D0*TGAUST-0.24D0*SI-0.161D0*AMN
     .           -0.115D0*ANI+0.25D0*CU+0.06D0*AMO+2.69D0*CR
            vfARmax=(cA-CCFB)/(cmAR-CCFB)*(1.D0-vfGr) 
            vfBFmax=(cmAR-cA)/(cmAR-CCFB)*(1.D0-vfGr)
            httime=httime+DTIMET/60.D0 ! in minutes
            apcss2=1.D0-exp(-AKAVRA*httime**AMAVRA)
            vfBF=apcss2*vfBFmax
            vfA=1.D0-(vfGr+vfBF)
            
c     cAR
c     for model 6 v3 with initial data form input file. I use without mass
c     convervation because there are problems for put a initial graphite volumen
c     that its mass couldn't be in equilibrium with matriz austenetine.
            if (NCOPC.eq.0) then
               cAR=cA+apcss2*(cmAR-cA)

c     for model 14 use mass conservation.
            else
               TempidGr=20.D0
               TempfdGr=1000.D0
               alphagr=9.26D-6  !alpha graphite
               cGr=100.D0       !carbon concentration graphite
               denGret=2200.D0  !graphite to enviroment temperature
               denA=8099.79D0-0.506D0*TGAUST+
     .              (-118.26D0+0.00739D0*TGAUST)*
     .              cAR-6.01D0*AMN !I use it because cFo and cFop aren't stable.
                                   !I can change it.
               denF=7875.96D0-0.297D0*TGAUST-5.62D-5*TGAUST**2.D0+
     .              (-206.35D0+0.00778D0*TGAUST+1.472D-6*TGAUST**2.D0)*
     .              CCFB-7.24D0*AMN
               denGrf=denGret/(3.D0*alphagr*(TempfdGr-TempidGr)+1.D0)
               denGr=(TGAUST-TempfdGr)/(TempidGr-TempfdGr)*denGret+
     .              (TGAUST-TempidGr)/(TempfdGr-TempidGr)*denGrf
               denT=vfGr*denGr+vfBF*denF+vfA*denA

               cAR=(C*denT-vfGr*cGr*denGr-vfBF*CCFB*denF)/(vfA*denA) !look
               if(cAR.lt.cA) cAR=cA !control
            endif
      endif

         TSOE1(IPLAT)=0.0D0
         TSOE2(IPLAT)=0.0D0
         HENER=0.0D0       
      ENDIF
         
c****************************************************************
c     DEFINES THE "MACROSCOPICAL" PROPERTIES (DENSITY, CAPACITY,
c     CONDUCTIVITY, PHASE-CHANGE FUNCTION AND LATENT HEAT)
c****************************************************************
      TSOE2(IPLAT)=TSOE2(IPLAT)*HENER ! L*f_pc at time t
      TSOE1(IPLAT)=TSOE1(IPLAT)*HENER ! L*f_pc at time t+dt

c     OPTIONS IN THE EVALUATION OF df_pc/dT
c     
c     1) standard df_pc/dT = delta f_pc / delta T (Diego's thesis)
c     2) standard df_pc/dT = delta f_pc / delta T (always ge 0)
c     3) df_pc/dT "exact"
c     4) df_pc/dT "exact" (always ge 0)
c     5) df_pc/dT = 0
c     
c     Notes:
c     L is constant
c     One possible option to obtain recalescense is to consider a
c     graphite/cementite-dependent latent heat value (L increases with
c     FG/FC)
c     
      IFPCDT=2                  ! better as input
      IF(IFPCDT.NE.3.AND.ICONVT.EQ.1)
     .     CALL RUNENDT('ERROR: IFPCDT NE 3 WHEN ICONVT=1')
      GO TO (1,2,3,4,5) IFPCDT
c     
    1 ICACOT=0
      VELTOT=1.D-10
      VELA1T=DABS(DTEMPT)
      IF(VELA1T.GT.VELTOT) THEN
         IF(IITERT.GT.0)
     .        TSOC1(IPLAT)=(TSOE1(IPLAT)-TSOE2(IPLAT))/(DTEMPT*DTIMET)
      ENDIF
      GO TO 10
c     
    2 ICACOT=1
      VELTOT=1.D-10
      VELA1T=DABS(DTEMPT)
      IF(VELA1T.GT.VELTOT) THEN
         IF(IITERT.GT.0)
     .        TSOC1(IPLAT)=(TSOE1(IPLAT)-TSOE2(IPLAT))/(DTEMPT*DTIMET)
      ENDIF
      GO TO 10
c     
    3 ICACOT=0
      TSOC1(IPLAT)=TSOCA+TSOCG+TSOCC
      TSOC1(IPLAT)=TSOC1(IPLAT)*HENER
      GO TO 10
c     
    4 ICACOT=1
      TSOC1(IPLAT)=TSOCA+TSOCG+TSOCC
      TSOC1(IPLAT)=TSOC1(IPLAT)*HENER
      GO TO 10
c     
    5 ICACOT=1
      TSOC1(IPLAT)=0.0D0
      GO TO 10
c     
 10   CONTINUE
      
      
c     transfer microscopical variables to ALPHAM and increments
c     ALPHAM index
c     version 1
      IF(IVERSI.EQ.1) THEN
         ALPHAM(IN+1)=FPCSS
         ALPHAM(IN+2)=FRACT
         ALPHAM(IN+3)=TTINU
         ALPHAM(IN+4)=TEMPD
         INUPC=INUPC+4
      ENDIF
c     version 2
      IF(IVERSI.EQ.2) THEN
       ALPHAM(IN+1)=FPCSS
       ALPHAM(IN+2)=FRAFN
       ALPHAM(IN+3)=DENFN
       ALPHAM(IN+4)=RADFN
       ALPHAM(IN+5)=FRAPN
       ALPHAM(IN+6)=DENPN
       ALPHAM(IN+7)=RADPN
       INUPC=INUPC+7
      ENDIF

c     version 3: austempered ductile iron model
      IF(IVERSI.EQ.3) THEN
c     not couple micro model
         if (NCOPC.eq.0) then
            ALPHAM(IN+1)=vfBF
            ALPHAM(IN+2)=vfA
            ALPHAM(IN+3)=cA
            ALPHAM(IN+4)=cAR
            ALPHAM(IN+5)=vfGr
            ALPHAM(IN+6)=httime
            INUPC=INUPC+6       !increments ALPHAM index. less than NBASES;
                                !see pointes.f
         endif
c     couple micro model (with micros14)
         if (NCOPC.ne.0) then
            ALPHAM(IX1)=vfGr
            ALPHAM(IX8)=vfA
            ALPHAM(IN+1)=vfBF
            ALPHAM(IN+2)=vfA
            ALPHAM(IN+3)=cA
            ALPHAM(IN+4)=cAR
            ALPHAM(IN+5)=vfGr
            ALPHAM(IN+6)=httime
            INUPC=INUPC+6       !increments ALPHAM index. less than NBASES;
                                !see pointes.f
         endif
      ENDIF
      
      RETURN
      END

