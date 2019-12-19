      SUBROUTINE MICROS7(TGAUAT,TGAUST,TGAUIT,DTEMPT,
     .                    BASMM, BASCC, BASKK,
     .                    TSOE2, TSOE1, TSOC1,
     .                    IPLAT,
     .                   ALPHAM,INUPC)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE PHASE-CHANGE FUNCTION, DENSITY,
C     CAPACITY & CONDUCTIVITY ACCORDING WITH THE MICROSTRUCTURAL MODEL
C     NUMBER 7 (IPCMO=7) OF RATE PHASE-CHANGE FORMULATIONS
C
C     CONCRETE HYDRATION
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
C             Second indexes (the same for any version) of ALPHAM
C             devoted to NNUPO variables (see pointes.f)
C             (variables to be transferred to the mechanical problem)
C
C     IN=INUPC (=0 in this case)
C
C     IVERSI=1
C     ALPHAM(IN+1): hydration degree (at time t+dt)
C     ALPHAM(IN+2): compression strength (at time t+dt)
C     ALPHAM(IN+3): -- (Young modulus?)
C     ALPHAM(IN+4): hydration amplitude (at time t+dt)
C 
C     IVERSI=2
C     ALPHAM(IN+1): hydration degree (at time t+dt)
C     ALPHAM(IN+2): compression strength (at time t+dt)
C     ALPHAM(IN+3): tension strength (at time t+dt)
C     ALPHAM(IN+4): Young modulus (at time t+dt)
C     ALPHAM(IN+5): rate of hydration degree (at time t+dt)
C     ALPHAM(IN+6): hydration amplitude (at time t+dt)
C     ALPHAM(IN+7): chemical afinity (at time t+dt)
C     ALPHAM(IN+8): chemical dissipation (at time t+dt)
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES (thermal-microstructural)
C
      INCLUDE 'nued_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION ALPHAM(NHISTM)
      DIMENSION BASKK(NDIMETO)
C
      DIMENSION TSOE1(5), TSOE2(5), TSOC1(5)
C
C**** TRANSFERS "VPLAT" TO MICROSCOPICAL VARIABLES
C
      IVERSI=INT(VPLAT(IPLAT,6))
C
      IF(IVERSI.EQ.1) THEN
      AM0=VPLAT(IPLAT,7)
      AK=VPLAT(IPLAT,8)
      VISC7=VPLAT(IPLAT,9)
      UNIVC=VPLAT(IPLAT,10)
      ENERA=VPLAT(IPLAT,11)
      HENER=VPLAT(IPLAT,12)      ! constant
      CHININF=VPLAT(IPLAT,13)
      CONTA=VPLAT(IPLAT,14)
      CONTB=VPLAT(IPLAT,15)
      IMECHMIC=INT(VPLAT(IPLAT,16))
C
      IV=16                      ! number of VPLAT defined in input data
C
C**** TRANSFERS "ALPHAM" TO MICROSCOPICAL VARIABLES (LOAD LAST CONVERGED
C     VALUES)
C
      IN=INUPC
      CHIO=ALPHAM(IN+1)          ! hydration degree
      FCO=ALPHAM(IN+2)           ! compression strength
      YOO=ALPHAM(IN+3)           ! -- (Young modulus?)
      RRO=ALPHAM(IN+4)           ! hydration amplitude
C
C**** SOLVE MODEL 7 (INTEGRATION OF HYDRATION EQUATION)
C
      CHIN=CHIO
      FCN=FCO
      YON=YOO
      RRN=RRO
c
c     Cambio para hacer visc7 variable segun la amplitud de hidratacion
c
      VISC77=VISC7+5.5E-05*CHIN
c      
      DCHI=0.0
      TOLCHI=1.0E-03
      BAUX=EXP(ENERA/(UNIVC*TGAUST))
      RECHI=(CHIN-CHIO)/DTIMET-(AM0-AK*CHIN)/(VISC77*BAUX)
C
      MAXCHI=100
      DO I=1,MAXCHI
       AJCHI=-1.0/DTIMET-(AK/(VISC77*BAUX))
       DCHI=RECHI/AJCHI
       CHIN=CHIN+DCHI
       RECHI=(CHIN-CHIO)/DTIMET-(AM0-AK*CHIN)/(VISC77*BAUX)
       IF(RECHI.LT.TOLCHI) GO TO 200
       IF(I.EQ.100)
     .  CALL RUNENDT('MODEL 7 DOES NOT CONVERGE') 
      ENDDO
  200 CONTINUE
C
C**** COMPUTES MECHANICAL PARAMETERS (as a function of CHIN)
C
      RRN=CHIN*CHININF
      FCN=CONTA*CHIN*CHIN+CONTB*CHIN
C
      TSOE1(IPLAT)=-CHIN                ! f_pc at time t+dt
      TSOE2(IPLAT)=-CHIO                ! f_pc at time t
C
C**** DEFINES THE "MACROSCOPICAL" PROPERTIES (DENSITY, CAPACITY,
C     CONDUCTIVITY, PHASE-CHANGE FUNCTION AND LATENT HEAT)
C

C
      TSOE2(IPLAT)=TSOE2(IPLAT)*HENER/BASMM   ! L*f_pc at time t
      TSOE1(IPLAT)=TSOE1(IPLAT)*HENER/BASMM   ! L*f_pc at time t+dt
C
C**** OPTIONS IN THE EVALUATION OF df_pc/dT
C
C     1) standard df_pc/dT = delta f_pc / delta T (Diego's thesis)
C     2) standard df_pc/dT = delta f_pc / delta T (always ge 0)
C     3) df_pc/dT "exact"
C     4) df_pc/dT "exact" (always ge 0)
C     5) df_pc/dT = 0
C
      IFPCDT=5                ! better as input
C
      IF(IFPCDT.NE.3.AND.ICONVT.EQ.1)
     . CALL RUNENDT('ERROR: IFPCDT NE 3 WHEN ICONVT=1')
C
      GO TO (1,2,3,4,5) IFPCDT
C
    1 ICACOT=0
      VELTOT=1.D-10
      VELA1T=DABS(DTEMPT)
      IF(VELA1T.GT.VELTOT) THEN
       IF(IITERT.GT.0)
     .  TSOC1(IPLAT)=(TSOE1(IPLAT)-TSOE2(IPLAT))/(DTEMPT*DTIMET)
       GO TO 10
      ENDIF
C
    2 ICACOT=1
      VELTOT=1.D-10
      VELA1T=DABS(DTEMPT)
      IF(VELA1T.GT.VELTOT) THEN
       IF(IITERT.GT.0)
     .  TSOC1(IPLAT)=(TSOE1(IPLAT)-TSOE2(IPLAT))/(DTEMPT*DTIMET)
       GO TO 10
      ENDIF
C
    3 ICACOT=0
      TSOC1(IPLAT)=DTIMET*(AM0-AK*CHIN)*ENERA/
     .             (VISC77*BAUX*UNIVC*TGAUST*TGAUST)
      TSOC1(IPLAT)=-TSOC1(IPLAT)*HENER/BASMM
      GO TO 10
C
    4 ICACOT=1                                        
      TSOC1(IPLAT)=DTIMET*(AM0-AK*CHIN)*ENERA/        
     .             (VISC77*BAUX*UNIVC*TGAUST*TGAUST)   
      TSOC1(IPLAT)=-TSOC1(IPLAT)*HENER/BASMM
      GO TO 10
C      
    5 ICACOT=1                                        
      TSOC1(IPLAT)=0.0
      GO TO 10                                        
C
   10 CONTINUE
C
C**** TRANSFER MICROSCOPICAL VARIABLES TO "ALPHAM"
C
      ALPHAM(IN+1)=CHIN
      ALPHAM(IN+2)=FCN
      ALPHAM(IN+3)=YON
      ALPHAM(IN+4)=RRN
C
C**** INCREMENTS "ALPHAM" INDEX
C
      INUPC=INUPC+4
      ENDIF                       ! iversi.eq.1
C
      IF(IVERSI.EQ.2) THEN
      AKBAR=VPLAT(IPLAT,7)
      ETA00=VPLAT(IPLAT,8)
      ETABAR=VPLAT(IPLAT,9)
      GHINF=VPLAT(IPLAT,10)
      ALPHABAR=VPLAT(IPLAT,11)
      AABAR=VPLAT(IPLAT,12)
      UNIVC=VPLAT(IPLAT,13)
      ENERA=VPLAT(IPLAT,14)
      HENER=VPLAT(IPLAT,15)      ! constant
      CONTA=VPLAT(IPLAT,16)
      CONTB=VPLAT(IPLAT,17)
      ENETE=VPLAT(IPLAT,18)
      FINFIN=VPLAT(IPLAT,19)
      FINFINM=VPLAT(IPLAT,20)
      GSETT=VPLAT(IPLAT,21)
      TMAX=VPLAT(IPLAT,22)
      TREF=VPLAT(IPLAT,23)
      CONAE=VPLAT(IPLAT,24)
      IMECHMIC=INT(VPLAT(IPLAT,25))
C
      IV=25                      ! number of VPLAT defined in input data
C
C**** TRANSFERS "ALPHAM" TO MICROSCOPICAL VARIABLES (LOAD LAST CONVERGED
C     VALUES)
C
      IN=INUPC
      GHO=ALPHAM(IN+1)           ! hydration degree
      FCO=ALPHAM(IN+2)           ! compression strength
      FCMO=ALPHAM(IN+3)          ! tension strength
      YOUNGO=ALPHAM(IN+4)        ! Young modulus
      RGHO=ALPHAM(IN+5)          ! rate of hydration degree
      CHIO=ALPHAM(IN+6)          ! hydration amplitude
      CHEMAO=ALPHAM(IN+7)        ! chemical afinity
      DISSO=ALPHAM(IN+8)         ! chemical dissipation
C
C**** SOLVE MODEL 7 (INTEGRATION OF HYDRATION EQUATION)
C
      GHN=GHO
      FCN=FCO
      FCMN=FCMO
      YOUNGN=YOUNGO
      RGHN=RGHO
      CHIN=CHIO
      CHEMAN=CHEMAO
      DISSN=DISSO
C
      IITER7=10                          ! better as input
      IF(IITERT.GT.IITER7) TGAUST=TGAUAT ! tuning to improve convergence
C
      DGH=0.0D0                          ! initialization
      REGH=0.0D0
      AJGH=0.0D0
      TOLGH=1.0E-05
      BAUX=EXP(ENERA/(UNIVC*TGAUST))
      ETA=ETA00
      AKBAR2=AKBAR/ETA
      CHEMAN=(AABAR+GHN)*(GHINF-GHN)-ALPHABAR*(TGAUST-TGAUIT)/TGAUIT
c     CHEMAN=AABAR*(GHINF-GHN)-ALPHABAR*(TGAUST-TGAUIT)/TGAUIT  ! simpl.
      IF(CHEMAN.GT.0.0) THEN
c      ETA=ETA00+ETABAR*GHN                ! linear
       ETA=ETA00*EXP(ETABAR*GHN/GHINF)     ! exponential
       AKBAR2=AKBAR/ETA
       REGH=-(GHN-GHO)/DTIMET+CHEMAN*AKBAR2/BAUX
C
       MAXGH=100
       DO I=1,MAXGH
        AJGH=-AKBAR2/BAUX*(GHINF-AABAR-2.0*GHN-
     .                     ETABAR/GHINF*CHEMAN)
c       AJGH=1.0/DTIMET+AKBAR2/BAUX*AABAR+
c    .       AKBAR*ETABAR/(ETA**2.0)*AABAR*(GHINF-GHN)/BAUX     ! simpl.
        IF(AJGH.LT.0.0) THEN
         AJGH=1.0D0/DTIMET
        ELSE
         AJGH=AJGH+1.0D0/DTIMET
        ENDIF
        DGH=REGH/AJGH
        GHN=GHN+DGH
C
        CHEMAN=(AABAR+GHN)*(GHINF-GHN)-ALPHABAR*(TGAUST-TGAUIT)/TGAUIT
c       CHEMAN=AABAR*(GHINF-GHN)-ALPHABAR*(TGAUST-TGAUIT)/TGAUIT !simpl.
        IF(CHEMAN.GT.0.0) THEN
c        ETA=ETA00+ETABAR*GHN                ! linear
         ETA=ETA00*EXP(ETABAR*GHN/GHINF)     ! exponential
         AKBAR2=AKBAR/ETA
         REGH=-(GHN-GHO)/DTIMET+CHEMAN*AKBAR2/BAUX
         IF(DABS(REGH).LT.TOLGH) GO TO 201
         IF(I.EQ.100)
     .    CALL RUNENDT('MODEL 7 DOES NOT CONVERGE') 
        ELSE
         CHEMAN=0.0
         GHN=GHN-DGH
         CALL RUNMENT('WARNING: CHEMAN < 0 DURING ITERATION')
         GO TO 201
        ENDIF
       ENDDO
  201  CONTINUE
C
      ELSE
       CHEMAN=0.0D0
      ENDIF
C
C**** COMPUTES MECHANICAL PARAMETERS (as a function of GHN)
C
      RGHN=(GHN-GHO)/DTIMET
      IF(GHN.LT.GSETT) THEN
       AMFCN=CONTA*GHN
      ELSE
       AMFCN=CONTA*GHN+CONTB
      ENDIF
      AMAFCN=AMFCN*((TMAX-TGAUST)/(TMAX-TREF))**(ENETE)
      IF(AMAFCN.GT.0.0D0) THEN
       RFCN=AMAFCN*RGHN*FINFIN
       FCN=RFCN*DTIMET+FCO                        ! compression strength
       YOUNGN=CONAE*DSQRT(FCN/FINFIN)             ! Young modulus
       AKAPPA=FCN/FINFIN
       FCMN=AKAPPA**(2.0D0/3.0D0)*FINFINM         ! tension strength
      ENDIF
C
c     CHIN=CHIO                        ! not used yet!!!!
c     DISSN=DISSO                      ! not used yet!!!!
c     FCN=CONTA*GHN*GHN+CONTB*GHN      ! old form
c     YOUNGN=CONAE*DSQRT(FCN/FINFIN)   ! old form
C
      TSOE1(IPLAT)=-GHN                ! f_pc at time t+dt
      TSOE2(IPLAT)=-GHO                ! f_pc at time t
C
C**** DEFINES THE "MACROSCOPICAL" PROPERTIES (DENSITY, CAPACITY,
C     CONDUCTIVITY, PHASE-CHANGE FUNCTION AND LATENT HEAT)
C

C
      TSOE2(IPLAT)=TSOE2(IPLAT)*HENER/BASMM   ! L*f_pc at time t
      TSOE1(IPLAT)=TSOE1(IPLAT)*HENER/BASMM   ! L*f_pc at time t+dt
C
C**** OPTIONS IN THE EVALUATION OF df_pc/dT
C
C     1) standard df_pc/dT = delta f_pc / delta T (Diego's thesis)
C     2) standard df_pc/dT = delta f_pc / delta T (always ge 0)
C     3) df_pc/dT "exact"
C     4) df_pc/dT "exact" (always ge 0)
C     5) df_pc/dT = 0
C
      IFPCDT=5                ! better as input
C
      IF(IFPCDT.NE.3.AND.ICONVT.EQ.1)
     . CALL RUNENDT('ERROR: IFPCDT NE 3 WHEN ICONVT=1')
C
      GO TO (11,12,13,14,15) IFPCDT
C
   11 ICACOT=0
      VELTOT=1.D-10
      VELA1T=DABS(DTEMPT)
      IF(VELA1T.GT.VELTOT) THEN
       IF(IITERT.GT.0)
     .  TSOC1(IPLAT)=(TSOE1(IPLAT)-TSOE2(IPLAT))/(DTEMPT*DTIMET)
       GO TO 20
      ENDIF
C
   12 ICACOT=1
      VELTOT=1.D-10
      VELA1T=DABS(DTEMPT)
      IF(VELA1T.GT.VELTOT) THEN
       IF(IITERT.GT.0)
     .  TSOC1(IPLAT)=(TSOE1(IPLAT)-TSOE2(IPLAT))/(DTEMPT*DTIMET)
       GO TO 20
      ENDIF
C
   13 ICACOT=0
      CORCHE=1.0D0                                   ! option 1
c     CORCHE=1.0D0/(DTIMET*AJGH)                     ! option 2
C
      TSOC1(IPLAT)=(GHN-GHO)*ENERA/(UNIVC*TGAUST*TGAUST)*CORCHE
      TSOC1(IPLAT)=TSOC1(IPLAT)*HENER/BASMM
      GO TO 20
C
   14 ICACOT=1                                        
      CORCHE=1.0D0                                   ! option 1
c     CORCHE=1.0D0/(DTIMET*AJGH)                     ! option 2
C
      TSOC1(IPLAT)=(GHN-GHO)*ENERA/(UNIVC*TGAUST*TGAUST)*CORCHE
      TSOC1(IPLAT)=TSOC1(IPLAT)*HENER/BASMM
      GO TO 20
C
   15 ICACOT=1                                        
      TSOC1(IPLAT)=0.0
      GO TO 20                                        
C
   20 CONTINUE
C
C**** TRANSFER MICROSCOPICAL VARIABLES TO "ALPHAM"
C
      ALPHAM(IN+1)=GHN
      ALPHAM(IN+2)=FCN
      ALPHAM(IN+3)=FCMN
      ALPHAM(IN+4)=YOUNGN
      ALPHAM(IN+5)=RGHN
      ALPHAM(IN+6)=CHIN
      ALPHAM(IN+7)=CHEMAN
      ALPHAM(IN+8)=DISSN
C
C**** INCREMENTS "ALPHAM" INDEX
C
      INUPC=INUPC+8
      ENDIF                       ! iversi.eq.2
C
      RETURN
      END
