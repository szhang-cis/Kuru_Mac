      SUBROUTINE MICROS11(TGAUAT,TGAUST,TGAUIT,DTEMPT,
     .     BASMM, BASCC, BASKK,
     .     TSOE2, TSOE1, TSOC1,
     .     IPLAT, INUPM,
     .     ALPHAM,INUPC)  
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE PHASE-CHANGE FUNCTION ACCORDING TO THE
C     MICROSTRUCTURAL MODEL NUMBER 11 (IPCMO=11)
C
C     GRAY CAST IRON MICROSTRUCTURAL MODEL FOR EUTECTOID TRANSFORMATION
C     (Alejandro Urrutia's (AU) master thesis, july 2013)
C
C***********************************************************************
C     
C     Index of variables
C     
C     Input:
C     TGAUAT= Temperature at time t
C     TGAUST= Temperature at time t+dt
C     TGAUIT= Initial temperature
C     DTEMPT= Temperature rate
C     
C     Output:
C     TSOE2 = each L*phase-change function at time t
C     TSOE1 = each L*phase-change function at time t+dt
C     
C     ALPHAM= array of microstructural (microscopical) variables
C     
C     ALPHAM(IN+1)=ferrite fraction
C     ALPHAM(IN+2)=ferrite grain density
C     ALPHAM(IN+3)=ferrite grain radius
C     ALPHAM(IN+4)=pearlite fraction
C     ALPHAM(IN+5)=pearlite grain density
C     ALPHAM(IN+6)=pearlite grain radius
C     ALPHAM(IN+7)=solidus temperature
C     
C     Microstructural variables from solidification (model 2)
C     ALPHAM(IPLAC(INUPM,ICOPC+1,3))=graphite eutectic fraction
C     ALPHAM(IPLAC(INUPM,ICOPC+1,4))=austenite fraction
C     ALPHAM(IPLAC(INUPM,ICOPC+1,5))=graphite fraction
C     ALPHAM(IPLAC(INUPM,ICOPC+1,6))=lamellar spacing
C     ALPHAM(IPLAC(INUPM,ICOPC+1,7))=lamellae thickness
C     ALPHAM(IPLAC(INUPM,ICOPC+1,8))=liquid fraction
C     ALPHAM(IPLAC(INUPM,ICOPC+1,9))=white eutectic fraction
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C     
C**** COUPLING VARIABLES    ! thermal-microstructural
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
      TWOPIT=6.283185307179586D0
C
C**** TRANSFERS "VPLAT" TO MICROSCOPICAL VARIABLES
C 
      IX=0
      NCOPC=IPLAC(INUPM,1,1)
      IF(NCOPC.EQ.0) IX=6
C
      TEFERI=VPLAT(IPLAT, 6+IX)
      HENERF=VPLAT(IPLAT, 7+IX)
      IFENU=INT(VPLAT(IPLAT, 8+IX))
      OMEGAF=VPLAT(IPLAT, 9+IX)
      IFEGR=INT(VPLAT(IPLAT,10+IX))
      RFER0=VPLAT(IPLAT,11+IX)
C
      TEPERI=VPLAT(IPLAT,12+IX)
      HENERP=VPLAT(IPLAT,13+IX)
      IPENU=INT(VPLAT(IPLAT,14+IX))
      OMEGAP=VPLAT(IPLAT,15+IX)
      IPEGR=INT(VPLAT(IPLAT,16+IX))
      DELTAP=VPLAT(IPLAT,17+IX)
      DIFFBC=VPLAT(IPLAT,18+IX)
      SIGMAP=VPLAT(IPLAT,19+IX)
C
      IFPCDT=INT(VPLAT(IPLAT,20+IX))
      IAFLOJ=INT(VPLAT(IPLAT,21+IX))
C
      IV=21+IX                   ! number of VPLAT defined in input data
C
      CA=VPLAT(IPLAT,IV+1)
      SI=VPLAT(IPLAT,IV+2)
      QM=VPLAT(IPLAT,IV+3)
C
      FGREUT=VPLAT(IPLAT,IV+4)   
C
      TEFER=VPLAT(IPLAT,IV+5)
      TEPER=VPLAT(IPLAT,IV+6)
      TALFA=VPLAT(IPLAT,IV+7)
      TCURIE=VPLAT(IPLAT,IV+8)
      CEUD=VPLAT(IPLAT,IV+9)
      CAMAX=VPLAT(IPLAT,IV+10)
C
C**** TRANSFERS "ALPHAM" TO MICROSCOPICAL VARIABLES
C
      IN=INUPC
C
      FSOLF=ALPHAM(IN+1)        ! ferrite fraction
      ANUCF=ALPHAM(IN+2)        ! ferrite grain density
      ARADF=ALPHAM(IN+3)        ! ferrite grain radius
      FSOLP=ALPHAM(IN+4)        ! pearlite fraction
      ANUCP=ALPHAM(IN+5)        ! pearlite grain density
      ARADP=ALPHAM(IN+6)        ! pearlite grain radius
      TSOLD=ALPHAM(IN+7)        ! solidus temperature
C
C**** DEALS WITH COUPLING BETWEEN MICROSTRUCTURAL MODELS (ALPHAM
C     POINTERS OF MODEL THAT PROVIDES INITIAL CONDITION; see miccoup.f)
C
      IPLX2=0
      NCOPC=IPLAC(INUPM,1,1)    ! only for NCOPC=1
      IF(NCOPC.EQ.0) THEN       ! initial cond. (ic) from this model
       CALL RUNENDT('ERROR= INIT. COND. FOR MICRO MODEL 11 NOT IMPLEM.')
      ELSE                      ! initial cond. (ic) from other model
       DO ICOPC=1,NCOPC
        IPLX1=IPLAC(INUPM,ICOPC+1,1)   ! ph-ch that provides ic
        IPLX2=IPLAC(INUPM,ICOPC+1,2)   ! model that provides ic
C
        IF(IPLX2.EQ.2) THEN     ! IPCMO=2 (Gray cast iron solid. model)
         IPLX3 =IPLAC(INUPM,ICOPC+1,3)
         IPLX4 =IPLAC(INUPM,ICOPC+1,4)
         IPLX5 =IPLAC(INUPM,ICOPC+1,5)
         IPLX6 =IPLAC(INUPM,ICOPC+1,6)
         IPLX7 =IPLAC(INUPM,ICOPC+1,7)
         IPLX8 =IPLAC(INUPM,ICOPC+1,8)
         IPLX9 =IPLAC(INUPM,ICOPC+1,9)
         FSOLG=ALPHAM(IPLX3)           ! graphite eutectic fraction
         FAUST=ALPHAM(IPLX4)           ! austenite fraction
         FGRAP=ALPHAM(IPLX5)           ! graphite fraction
         AGLAS=ALPHAM(IPLX6)           ! lamellar spacing
         AGLAT=ALPHAM(IPLX7)           ! lamellae thickness
         FLIQD=ALPHAM(IPLX8)           ! liquid fraction
         FSOLC=ALPHAM(IPLX9)           ! white eutectic fraction
        ENDIF
       ENDDO                    ! 1,ncopc
      ENDIF                     ! ncopc=0
C
C**** CHECKS LIQUID FRACTION AND OBTAINS SOLIDUS TEMPERATURE
C
      IF(FLIQD.EQ.0.0D0) THEN
       IF(TSOLD.EQ.0.0D0) THEN
        TSOLD=TGAUST            ! solidus temperature
       ENDIF
      ELSE
       INUPC=INUPC+7            ! less than NBASES; see pointes.f
       RETURN
      ENDIF
C
C**** LOAD LAST CONVERGED VALUES
C
      FSOLFO=FSOLF
      FSOLPO=FSOLP
C
C**** ESTABLISHES VARIABLES OF MODEL 11
C
      FAUST0=(1.0D0-FGREUT)*FSOLG      ! initial austenite fraction
      FGRAP0=FGREUT*FSOLG              ! initial graphite fraction
C
      AGLAT0=FGREUT*AGLAS              ! initial lamellae thickness
C
C**** SOLVES MODEL 11
C
      TSOCG= 0.0D0              ! only useful to thermal jacobian matrix
      TSOCFN=0.0D0
      TSOCFR=0.0D0
      TSOCPN=0.0D0
      TSOCPR=0.0D0
C
C**** EUTECTIC-EUTECTOID CARBON DIFFUSION
C
      IF(TGAUST.LT.TSOLD.AND.TGAUST.GE.TALFA) THEN
       FGREE=(TGAUST-TSOLD)/(TALFA-TSOLD)*(2.2D0-CEUD)*0.01D0*
     .        7300.0D0/2300.0D0
       FGRAP=FGRAP0+FGREE
       FAUST=FAUST0-FGREE
       TSOCG=1.0D0/(TALFA-TSOLD)*(2.2D0-CEUD)*0.01D0*
     .        7300.0D0/2300.0D0
      ENDIF
C
      FGRAP1=0.0D0
      FAUST1=0.0D0
      FEUTEC=0.0D0
      IF(TGAUST.LE.TALFA.OR.TGAUST.LE.TEPER) THEN
       FGRAP1=FGRAP0+ ! initial graphite fraction for eutectoid transf.
     .  (2.2D0-CEUD)*0.01D0*7300.0D0/2300.0D0    ! FGREE at TALFA
C
       FAUST1=FAUST0- ! initial austenite fraction for eutectoid transf.
     .  (2.2D0-CEUD)*0.01D0*7300.0D0/2300.0D0    ! FGREE at TALFA
C
       FEUTEC=FAUST1-FAUST              ! eutectoid fraction
       IF(FEUTEC.LT.0.0D0) THEN         ! control
        FEUTEC=0.0D0
       ENDIF
      ENDIF
C
C**** FERRITE EUTECTOID PHASE CHANGE
C
C     GUGF=TEFER-TGAUST     ! eutectoid stable (ferrite) undercooling
      GUGF=TALFA-TGAUST     ! eutectoid stable (ferrite) undercooling-AU
      ICORREC=0             ! not used
C
      IF(GUGF.GT.0.0D0.AND.FSOLF.LT.FAUST1.AND.FEUTEC.LT.FAUST1) THEN
       ANUCF=OMEGAF*FGRAP1/AGLAS              ! instantaneous nucleation
       TSOCFN=0.0D0
C
       FSOLF=1.0D0/3.0D0*TWOPIT*ANUCF*RFER0**3.0D0             ! control
       FEUTECX=FSOLF+FSOLP+FGRAP-FGRAP1
       IF(FEUTECX.GT.FAUST1)
     .  CALL RUNENDS('ERROR=NUCLEATION VALUE NOT COMPATIBLE WITH RFER0')
C
       CGAAL=18.76D0-4.112D-2*TGAUST+2.26D-5*TGAUST*TGAUST+0.125D0*SI
       CALGA=6.7D-2-5.0D-8*TGAUST*TGAUST-2.8D-5*TGAUST+
     .       (1.2D-8*TGAUST*TGAUST-4.5D-3)*SI
       CALGR=CALGA-(TALFA-TGAUST)*(2.9D-4-2.8D-5*SI)       ! TALFA-AU
C      CALGR=CALGA-(TEFER-TGAUST)*(2.9D-4-2.8D-5*SI)       ! TALFA=TEFER
       DALC=2.0D-6*DEXP(-10115.0D0/(TGAUST+273.15D0))*
     .      DEXP(0.5898*(1.0D0+4.0D0/TWOPIT*
     .      DATAN(15629.D0/TCURIE-15309.0D0/(TGAUST+273.15D0))))
       DGAC=(1.46D-5-0.36D-5*CA+0.509D-5*SI-0.315D-5*QM)*
     .      DEXP(-(144.3D0-15.0D0*CA+0.37D0*CA*CA+4.0507D0*SI-
     .               4.3663D0*QM)/(8.314D-3*(TGAUST+273.15D0)))
       XXXX=1.0D0/(ARADF*(CGAAL-CALGA))*
     .      (DALC*AGLAT*(CALGA-CALGR)/(2.0D0*ARADF-AGLAT)+
     .       DGAC*(CGAAL-CEUD)*(1.0D0-2.0D0*ARADF/AGLAS))
       ARADF=ARADF+XXXX*DTIMET
       TSOCFR=0.0D0         ! to be improved
C
       FSOLF=1.0D0/3.0D0*TWOPIT*ANUCF*ARADF**3.0D0
       FEUTECX=FSOLF+FSOLP+FGRAP-FGRAP1                        ! control
       IF(FEUTECX.GT.FAUST1) THEN
        FSOLF=FAUST1-FSOLP-FGRAP+FGRAP1
        ARADF=(3.0D0*FSOLF/(TWOPIT*ANUCF))**(1.0D0/3.0D0)
       ENDIF
C
       XXXX=4.0D0*DALC*(CALGA-CALGR)*ARADF*7510.D0/
     .      (AGLAT*(100.0D0-CALGR)*(2.0D0*ARADF-AGLAT)*2300.0D0)
       AGLAT=AGLAT+XXXX*DTIMET
C
       FGRAP=AGLAT/AGLAT0*FGRAP1
       TSOCG=0.0D0          ! to be improved
C
       FEUTECX=FSOLF+FSOLP+FGRAP-FGRAP1                        ! control
       IF(FEUTECX.GT.FAUST1) THEN
        FGRAP=FAUST1-FSOLF-FSOLP+FGRAP1
        AGLAT=FGRAP/FGRAP1*AGLAT0
       ENDIF
C
       ICORREC=1
      ENDIF                            ! gugf.gt.0.0.and....
C
      IF(ICORREC.EQ.1) THEN            ! correction
       IF(IAFLOJ.EQ.1) THEN            ! strict conditions
        FEUTEC=FSOLF+FSOLP+FGRAP-FGRAP1
       ENDIF                           ! iafloj.eq.1
      ENDIF                            ! icorrec.eq.1
C
      TSOCFN=TSOCFN*TWOPIT/3.0D0*ARADF**3.0D0
      TSOCFR=TSOCFR*TWOPIT*ANUCF*ARADF**2.0D0
C
C**** PEARLITE EUTECTOID PHASE CHANGE
C
      GUGP=TEPER-TGAUST     ! eutectoid metastable (pearlite) underc.
      ICORREC=0             ! not used
C
      IF(GUGP.GT.0.0D0.AND.FSOLP.LT.FAUST1.AND.FEUTEC.LT.FAUST1) THEN
       FXX=FAUST1           ! aproximate austenite fraction at TEPER
       ANUCP=OMEGAP*2.0D0*FXX/(TWOPIT*AGLAS)  ! instantaneous nucleation
       TSOCPN=0.0D0
C
       CGATH=(TGAUST-491.27D0)/365.95D0
       CTH=6.67D0
       ENTAL=-4.4D3*TGAUST*TGAUST+3.8D6*TGAUST+1.7D8
       SSC=2.0D0*SIGMAP*TEPER/(GUGP*ENTAL)
       SSP=2.0D0*SSC
       SSAL=7.0D0/8.0D0*SSP
       SSTH=1.0D0/8.0D0*SSP
       XXXX=(CGAAL-CGATH)/(CTH-CAMAX)*
     .      (2.0D0*DGAC+12.0D0*DELTAP*DIFFBC/SSP)*SSP/(SSAL*SSTH)*
     .      (1.0D0-SSC/SSP)*(1.0D0-FSOLP)
       ARADP=ARADP+XXXX*DTIMET
       TSOCPR=0.0D0         ! to be improved
C
       IF(CGAAL-CGATH.LE.0) CALL RUNENDS('CGAAL-CGATH<=0/')    ! control
       IF(CTH-CAMAX.LE.0) CALL RUNENDS('CTH-CAMAX<=0/')        ! AU
       IF(SSP.LE.0) CALL RUNENDS('SSP<=0/')
       IF(ENTAL.LE.0) CALL RUNENDS('ENTAL<=0')
       IF(DGAC.LE.0) CALL RUNENDS('DGAC<=0')
C
       FSOLP=1.0D0/3.0D0*TWOPIT*ANUCP*ARADP**3.0D0
       FEUTECX=FSOLF+FSOLP+FGRAP-FGRAP1                        ! control
       IF(FEUTECX.GT.FAUST1) THEN
        FSOLP=FAUST1-FSOLF-FGRAP+FGRAP1
        ARADP=(3.0D0*FSOLP/(TWOPIT*ANUCP))**(1.0D0/3.0D0)
       ENDIF
C
       ICORREC=1
      ENDIF                            ! gugp.gt.0.0.and....
C
      IF(ICORREC.EQ.1) THEN            ! correction
       IF(IAFLOJ.EQ.1) THEN            ! strict conditions
        FEUTEC=FSOLF+FSOLP+FGRAP-FGRAP1
       ENDIF                           ! iafloj.eq.1
      ENDIF                            ! icorrec.eq.1
C
      TSOCPN=TSOCPN*TWOPIT/3.0D0*ARADP**3.0D0
      TSOCPR=TSOCPR*TWOPIT*ANUCP*ARADP**2.0D0
C
C**** AUSTENITE FRACTION
C
      IF(TGAUST.LT.TALFA.OR.TGAUST.LT.TEPER) THEN
       FEUTEC=FSOLF+FSOLP+FGRAP-FGRAP1
       FAUST=FAUST1-FEUTEC
      ENDIF
C
C**** ESTABLISHES LATENT HEAT * PHASE-CHANGE FUNCTION (at time t & t+dt)
C
      TSOE2(IPLAT)=-HENERF*FSOLFO-HENERP*FSOLPO         ! L*f_pc at t
      TSOE1(IPLAT)=-HENERF*FSOLF-HENERP*FSOLP           ! L*f_pc at t+dt
C     
C**** OPTIONS IN THE EVALUATION OF df_pc/dT
C     
C     1) standard df_pc/dT = delta f_pc / delta T (Diego's thesis)
C     2) standard df_pc/dT = delta f_pc / delta T (always ge 0)
C     3) df_pc/dT "exact"
C     4) df_pc/dT "exact" (always ge 0)
C     5) df_pc/dT = 0
C     
C     Notes:
C     Only one value of ICACOT can be defined, i.e., ICACOT does not
C     depend of any particular phase-change. This can not be useful for
C     more than one microstructural phase-changes.
C     
      IF(IFPCDT.NE.3.AND.ICONVT.EQ.1)
     .     CALL RUNENDT('ERROR: IFPCDT NE 3 WHEN ICONVT=1')
C     
      GO TO (1,2,3,4,5) IFPCDT
C     
 1    ICACOT=0
      VELTOT=1.D-10
      VELA1T=DABS(DTEMPT)
      IF(VELA1T.GT.VELTOT) THEN
         IF(IITERT.GT.0)
     .        TSOC1(IPLAT)=(TSOE1(IPLAT)-TSOE2(IPLAT))/(DTEMPT*DTIMET)
      ENDIF
      GO TO 10
C     
 2    ICACOT=1
      VELTOT=1.D-10
      VELA1T=DABS(DTEMPT)
      IF(VELA1T.GT.VELTOT) THEN
         IF(IITERT.GT.0)
     .        TSOC1(IPLAT)=(TSOE1(IPLAT)-TSOE2(IPLAT))/(DTEMPT*DTIMET)
      ENDIF
      GO TO 10
C     
 3    ICACOT=0
      TSOC1(IPLAT)=TSOCG+TSOCFN+TSOCFR+TSOCPN+TSOCPR
      TSOC1(IPLAT)=TSOC1(IPLAT)*HENER
      GO TO 10
C     
 4    ICACOT=1
      TSOC1(IPLAT)=TSOCG+TSOCFN+TSOCFR+TSOCPN+TSOCPR
      TSOC1(IPLAT)=TSOC1(IPLAT)*HENER
      GO TO 10
C     
 5    ICACOT=1
      TSOC1(IPLAT)=0.0D0
      GO TO 10
C     
 10   CONTINUE
C
C**** TRANSFER MICROSCOPICAL VARIABLES TO "ALPHAM"
C
      ALPHAM(IN+1)=FSOLF        ! ferrite fraction
      ALPHAM(IN+2)=ANUCF        ! ferrite grain density
      ALPHAM(IN+3)=ARADF        ! ferrite grain radius
      ALPHAM(IN+4)=FSOLP        ! pearlite fraction
      ALPHAM(IN+5)=ANUCP        ! pearlite grain density
      ALPHAM(IN+6)=ARADP        ! pearlite grain radius
      ALPHAM(IN+7)=TSOLD        ! solidus temperature
C
      IF(IPLX2.EQ.2) THEN       ! IPCMO=2 (Gray cast iron solid. model)
       ALPHAM(IPLX4)=FAUST      ! only those that change
       ALPHAM(IPLX5)=FGRAP
       ALPHAM(IPLX7)=AGLAT
      ENDIF
C
C**** INCREMENTS "ALPHAM" INDEX
C
      INUPC=INUPC+7             ! less than NBASES; see pointes.f
C
      RETURN
      END
