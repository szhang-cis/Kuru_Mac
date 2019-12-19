      SUBROUTINE MICROS12(TGAUAT,TGAUST,TGAUIT,DTEMPT,
     .     BASMM, BASCC, BASKK,
     .     TSOE2, TSOE1, TSOC1,
     .     IPLAT, INUPM,
     .     ALPHAM,INUPC)  
C***********************************************************************
C     
C**** THIS ROUTINE EVALUATES THE PHASE-CHANGE FUNCTION ACCORDING TO THE
C     MICROSTRUCTURAL MODEL NUMBER 12 (IPCMO=12)
C     
C     S.G. CAST IRON MICROSTRUCTURAL MODEL FOR EUTECTOID TRANSFORMATION
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
C     New microstructural variables (eutectoid)
C     
C     ALPHAM(IN+1)= ferrite fraction
C     ALPHAM(IN+2)= pearlite fraction
C     
C     NNUM4T= number of different grain densities and radii (limit number)
C     ALPHAM(IN+3:IN+3+(NNUM4T-1))= ferrite grain radius
C     ALPHAM(IN+4+(NNUM4T-1):IN+4+2*(NNUM4T-1))= pearlite grain radius
C     ALPHAM(IN+5+2*(NNUM4T-1):IN+5+3*(NNUM4T-1))= pearlite grain densities
C     
C     Auxiliar microstructural variables (not printed; see pointes.f):
C     ALPHAM(IN+6+3*(NNUM4T-1)) = initial carbon content in austenite
C     ALPHAM(IN+7+3*(NNUM4T-1))= initial austenite fraction
C     ALPHAM(IN+8+3*(NNUM4T-1))= initial relation ferrite/austenite
C     initial pearlites colonies densities (instantaneus nucleation) 
C     ALPHAM(IN+9+3*(NNUM4T-1))= 
C     ALPHAM(IN+10+3*(NNUM4T-1)) = pearlite nucleation index
C     index: nucleation up to Tmin
C     ALPHAM(IN+11+3*(NNUM4T-1))= recalescence or T_min
C     
C     Old microstructural variables (solidification)
C     
C     NNUM4T=number of different grain densities and radii
C     
C     Model 4 (Boeri)
C     ALPHAM(IPLAC(INUPM,ICOPC+1, 3))= austenite fraction
C     ALPHAM(IPLAC(INUPM,ICOPC+1, 5))= graphite fraction
C     ALPHAM(IPLAC(INUPM,ICOPC+1, 6):IPLAC(INUPM,ICOPC+1, 6)+(NNUM4T-1))
C     = graphite grain density
C     ALPHAM(IPLAC(INUPM,ICOPC+1, 7):IPLAC(INUPM,ICOPC+1, 7)+(NNUM4T-1))
C     = graphite grain radius
C     ALPHAM(IPLAC(INUPM,ICOPC+1,11))= graphite nucleation index (jo)
C     Model 5 (Su)
C     ALPHAM(IPLAC(INUPM,ICOPC+1, 3))= austenite fraction
C     ALPHAM(IPLAC(INUPM,ICOPC+1, 4):IPLAC(INUPM,ICOPC+1, 4)+(NNUM4T-1))
C     = first austenite gr. R
C     ALPHAM(IPLAC(INUPM,ICOPC+1, 5))= graphite fraction
C     ALPHAM(IPLAC(INUPM,ICOPC+1, 6):IPLAC(INUPM,ICOPC+1, 6)+(NNUM4T-1))
C     = graphite grain density
C     ALPHAM(IPLAC(INUPM,ICOPC+1, 7):IPLAC(INUPM,ICOPC+1, 7)+(NNUM4T-1))
C     = graphite grain radius
C     ALPHAM(IPLAC(INUPM,ICOPC+1,11))= graphite nucleation index (jo)
C     s
C     Model 9 (Thesis)
C     ALPHAM(IPLAC(INUPM,ICOPC+1, 3))= austenite fraction
C     ALPHAM(IPLAC(INUPM,ICOPC+1, 4))= total gr. radius
C     ALPHAM(IPLAC(INUPM,ICOPC+1, 5))= graphite fraction
C     ALPHAM(IPLAC(INUPM,ICOPC+1, 6):IPLAC(INUPM,ICOPC+1, 6)+(NNUM4T-1))
C     = graphite grain density
C     ALPHAM(IPLAC(INUPM,ICOPC+1,10):IPLAC(INUPM,ICOPC+1,10)+(NNUM4T-1))
C     = graphite grain radius
C     ALPHAM(IPLAC(INUPM,ICOPC+1,11))= graphite nucleation index (jo)
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
      TWOPIT= 8.0D0*DATAN(1.0D0)
C     
C     
C***  TRANSFERS "VPLAT" TO MICROSCOPICAL VARIABLES
      TEA1F = VPLAT(IPLAT, 6)   ! stable temperature (ferrite)
      TEA1P = VPLAT(IPLAT, 7)   ! metastable temperature (pearlite)
      HENERF= VPLAT(IPLAT, 8)   ! latent heat of ferrite (eutectoid)
      HENERP= VPLAT(IPLAT, 9)   ! latent heat of pearlite (eutectoid)
      GNF   = VPLAT(IPLAT,10)   ! kinetics params of ferrite (growth)
      RINIFE= VPLAT(IPLAT,11)   ! diffusion coef. C in ferrite
C     
C***  ferrite growth values from input data file *-s.dat
      IGRMFX= INT(VPLAT(IPLAT,12)) ! ferrite's growth model
      IF(IGRMFX.EQ.1.OR.IGRMFX.EQ.2.OR.IGRMFX.EQ.3) THEN
         DENSA = VPLAT(IPLAT,13) ! austenite density
         DENSG = VPLAT(IPLAT,14) ! ferrite density
c$$$  DIFCA = VPLAT(IPLAT,15) ! diffusion coef. C in austenite
         DIFCFF= VPLAT(IPLAT,15) ! diffusion coef. C in ferrite
         HENERF= VPLAT(IPLAT,16) ! latent heat of ferrite (eutectoid)
      ENDIF
C     
C***  pearlite nucleation values from input data file *-s.dat
      INUCMPX= INT(VPLAT(IPLAT,17)) ! pearlite's nucleation model index
      IF(INUCMPX.EQ.1.OR.INUCMPX.EQ.2.OR.INUCMPX.EQ.3) THEN
         GNUCCP = VPLAT(IPLAT,18) ! nucleation coefficient of pearlite
         GNUCEP = VPLAT(IPLAT,19) ! nucleation exponent of pearlite
         IGNUCOP= INT(VPLAT(IPLAT,20)) ! pearlite law of nucleation
      ENDIF
C     nucleation arrest criterion index for pearlite
      INUCAX= INT(VPLAT(IPLAT,21)) 
C     
C***  pearlite growth  values from input data file *-s.dat
C     
      IGROMPX= INT(VPLAT(IPLAT,22)) ! pearlite's growth model index
      IF(IGROMPX.EQ.1.OR.IGROMPX.EQ.2.OR.IGROMPX.EQ.3) THEN
         GROCP = VPLAT(IPLAT,23) ! growth coefficient of pearlite's law
         GROEP = VPLAT(IPLAT,24) ! growth exponent of pearlite's law
         HENERP= VPLAT(IPLAT,25) ! latent heat of pearlite
         IGROMP= INT(VPLAT(IPLAT,26)) ! pearlite's growth model
      ENDIF
C     
      IKMICX= INT(VPLAT(IPLAT,27)) ! index for micro.-dep. conductivity
      IKAUX = 0
      IF(IKMICX.EQ.1.OR.IKMICX.EQ.2.OR.IKMICX.EQ.3) THEN
         IKAUX=3
         BASKS= VPLAT(IPLAT,28) ! solid conductivity
         BASKM= VPLAT(IPLAT,29) ! mushy conductivity
         BASKL= VPLAT(IPLAT,30) ! liquid conductivity
      ENDIF
C     
      IFPCDT= INT(VPLAT(IPLAT,28+IKAUX)) ! index for temp. derivative
      IAFLOJ= INT(VPLAT(IPLAT,29+IKAUX)) ! index for fraction correction
C     
      IV= 29+IKAUX              ! number of VPLAT defined in input data
C     
C**** DEALS WITH COUPLING BETWEEN MICROSTRUCTURAL MODELS
C     
      NCOPC= IPLAC(INUPM,1,1)
      IF(NCOPC.NE.0) THEN       ! initial condition from other model
         IPLX2= IPLAC(INUPM,2,2) ! model that provides initial condition
C     
C**** ALPHAM POINTERS OF MODEL THAT PROVIDE INITIAL CONDITION
C     (see miccoup.f)
         IPLX3 = IPLAC(INUPM,2,3) ! austenite fraction
         IPLX4 = IPLAC(INUPM,2,4) ! depend of solid. model
         IPLX5 = IPLAC(INUPM,2,5) ! graphite fraction
         IPLX6 = IPLAC(INUPM,2,6) ! su & boeri - graphite densities
         IPLX7 = IPLAC(INUPM,2,7) ! su & boeri - graphite radius
         IPLX8 = IPLAC(INUPM,2,8) ! dgc - first graphite gr. Nz3 
         IPLX9 = IPLAC(INUPM,2,9) ! dgc - graphite densities
         IPLX10= IPLAC(INUPM,2,10) ! dgc - graphite radius
         IPLX11= IPLAC(INUPM,2,11) ! jo nucleation index
         IPLX12= IPLAC(INUPM,2,12) ! eutectic cell radius (RADC)
         IPLX13= IPLAC(INUPM,2,13) ! silicon content
C     secondary dendrite arm spacing -- SDAS 
         IPLX14= IPLAC(INUPM,2,14) ! SDAS -- secondary dendrite arm spacing
         IPLX15= IPLAC(INUPM,2,15) ! copper content
         IPLX16= IPLAC(INUPM,2,16) ! manganese content
         IPLX17= IPLAC(INUPM,2,17) ! magnesium content
         IPLX18= IPLAC(INUPM,2,18) ! niobium content
         IPLX19= IPLAC(INUPM,2,19) ! tin content
         IPLX20= IPLAC(INUPM,2,20) ! chromium content
         IPLX21= IPLAC(INUPM,2,21) ! molybdenum content
         IPLX22= IPLAC(INUPM,2,22) ! nickel content
         IPLX23= IPLAC(INUPM,2,23) ! cpro
      ENDIF                     ! ncopc.ne.0
C     index of graphite/ferrite nucleation & pearlite colonies
      JO= INT(ALPHAM(IPLX11))   ! graphite nucl. index (numb. of families)
C     
C**** TRANSFERS "ALPHAM" TO MICROSCOPICAL VARIABLES
      IN= INUPC
c     
      FSOLF= ALPHAM(IN+1)       ! ferrite fraction
      FSOLP= ALPHAM(IN+2)       ! pearlite fraction
      CAG  = ALPHAM(IN+3)       ! carbon content in austenite in contact with graphite
      CPROS= ALPHAM(IN+4)       ! carbon content in austenite (mass equilibrium)
C     initialize ferrite grain radius
      IF(JO.GT.0) THEN
         DO INU=1,JO
            DELTAG(INU)= ALPHAM(IN+5+INU-1) ! layer of graphite nodules
            DELTAF(INU)= ALPHAM(IN+6+NNUM4T-1+INU-1) ! layer of ferrite grains
            RAFGGN(INU)= ALPHAM(IN+7+2*(NNUM4T-1)+INU-1) ! area of nodules cubiertos
            FNU(INU)   = ALPHAM(IN+8+3*(NNUM4T-1)+INU-1) ! ferrite grains radius
         ENDDO
      ENDIF
      JP= INT(ALPHAM(IN+15+6*(NNUM4T-1))) ! pearlite nucl. index
C     initialize pearlite grain densities and pearlite nodules radius
      IF(JP.GT.0) THEN
         DO JNU=1,JP
            DPERC(JNU,1)= ALPHAM(IN+ 9+4*(NNUM4T-1)+JNU-1) ! pear. gr. dens.
            PNUC (JNU)  = ALPHAM(IN+10+5*(NNUM4T-1)+JNU-1) ! pear. gr. rad. (Rp)
         ENDDO
      ENDIF
C     initial carbon content in austenite
C     carbon content corresponding to S' point in Fe-C-Si diagram
      CSP     = ALPHAM(IN+11+6*(NNUM4T-1))
C     initial fraction of austenite
      FSOLAINI= ALPHAM(IN+12+6*(NNUM4T-1))
C     relation of austenite fraction transformed
      RFA     = ALPHAM(IN+13+6*(NNUM4T-1))
C     initial pearlites colonies densities (for instantaneus nucleation)
CCC   DPERI   = ALPHAM(IN+14+6*(NNUM4T-1))
C     initial pearlites colonies densities (for instantaneus nucleation)
      CFGO    = ALPHAM(IN+14+6*(NNUM4T-1))
C     index to set TSP and CAFO
      ICEp    = INT(ALPHAM(IN+16+6*(NNUM4T-1))) 
C     temperature corresponding to S prime point in FE-C-Si phase diagram
      TSp     = ALPHAM(IN+18+6*(NNUM4T-1)) ! temperature for final solidification
C     content corresponding to S prime point in FE-C-Si phase diagram
      CAFO    = ALPHAM(IN+17+6*(NNUM4T-1)) 
C     graphite volumetric fraction
c$$$  UVVGR   = ALPHAM(IN+19+6*(NNUM4T-1)) 
C     index to set RFA for eutectoid transformation
      ICEpp   = INT(ALPHAM(IN+19+6*(NNUM4T-1)))
C     index to set CFG (if it is lower than one)
      ICFG    = INT(ALPHAM(IN+20+6*(NNUM4T-1)))
C     
      IF(INUCAX.EQ.1)
     .     INDEXG= INT(ALPHAM(IN+21+6*(NNUM4T-1))) ! lt NBASES, see pointes.f
      IF(INUCAX.EQ.2)
     .     TINDEXG= ALPHAM(IN+21+6*(NNUM4T-1)) ! lt NBASES, see pointes.f
C     
      IF(IPLX2.EQ.4) THEN       ! IPCMO=4 (Boeri's model)
         FSOLA = ALPHAM(IPLX3)
         FSOLG = ALPHAM(IPLX5)
         RADT  = ALPHAM(IPLX4)     
C     
         IF(JO.GT.0) THEN
            IF(FSOLF.EQ.0.0D0) THEN ! first time
               DO INU=1,JO
                  FNU(INU)= 1.0D-5 ! Rf= 1 micrometro
                  RFA=0.0D0
               ENDDO
            ENDIF
            DO INU=1,JO
               DNU(INU,1)= ALPHAM(IPLX6+INU-1)
               RNU(INU)  = ALPHAM(IPLX7+INU-1)
            ENDDO
         ENDIF
      ENDIF
C     
      IF(IPLX2.EQ.5) THEN       ! IPCMO=5  (SU' model)
         FSOLA= ALPHAM(IPLX3)
         FSOLG= ALPHAM(IPLX5)
         RADT = ALPHAM(IPLX4)
C     
         IF(JO.GT.0) THEN
            IF(FSOLF.EQ.0.0D0) THEN ! first time
               DO INU=1,JO
                  FNU(INU)= 1.01D0*ALPHAM(IPLX10+INU-1) ! Rf= 1 micrometro
                  RFA=0.0D0
               ENDDO
            ENDIF
            DO INU=1,JO
               DNU(INU,1)= ALPHAM(IPLX6+INU-1)
               RNU(INU)  = ALPHAM(IPLX7+INU-1)
            ENDDO
         ENDIF
      ENDIF
C     
      IF(IPLX2.EQ.9) THEN       ! IPCMO=9  (tesis' model)
         FSOLA= ALPHAM(IPLX3)   ! austenita fraction
         FSOLG= ALPHAM(IPLX5)   ! graphite fraction
         RADT = ALPHAM(IPLX4)   ! eutectic grain radius
         SICON= ALPHAM(IPLX13)  ! silicon content
         CUCON= ALPHAM(IPLX15)  ! copper content
         QNCON= ALPHAM(IPLX16)  ! manganses content
         QOCON= ALPHAM(IPLX21)  ! molybdenum content
         CRCON= ALPHAM(IPLX20)  ! chromium content
         QICON= ALPHAM(IPLX22)  ! nickel content
         CPRO = ALPHAM(IPLX23)  ! cpro
C     
         IF(JO.GT.0) THEN
            IF(FSOLF.EQ.0.0D0) THEN ! first time
               DO INU=1,JO
                  FNU(INU)= RINIFE
                  RFA     = 0.0D0
               ENDDO
            ENDIF
            DO INU=1,JO
               DNU(INU,1)= ALPHAM(IPLX6+INU-1)
               RNU(INU)  = ALPHAM(IPLX10+INU-1)
            ENDDO
         ENDIF
      ENDIF
C     
C**** LOAD LAST CONVERGED VALUES
C     
      FSOLGO= FSOLG
      FSOLAO= FSOLA
      FSOLFO= FSOLF
      FSOLPO= FSOLP
      CPROSO= CPROS
      CAGO  = CAG
      RFAO  = RFA
      RFA   = 0.0D0
C**   initial porcent of Si - manganese an copper 
C     (for PhaseDiagram calculation)
      SI= VPLAT((IPLAT-1),11)   ! second-phase change-1 (first phase change)
      PH= VPLAT((IPLAT-1),16)
      CU= VPLAT((IPLAT-1),21)   ! second-phase change-1 (first phase change)
      QN= VPLAT((IPLAT-1),26)   ! second-phase change-1 (first phase change)
      CR= VPLAT((IPLAT-1),46)   ! second-phase change-1 (first phase change)
      QO= VPLAT((IPLAT-1),51)   ! second-phase change-1 (first phase change)
      QI= VPLAT((IPLAT-1),56)   ! second-phase change-1 (first phase change)
C     
C**** COMPUTES EQUILIBRIUM EUTECTOID TEMPERATURE & OTHER PARAMETERS
C     
C     Nomenclature:
C     TEaTF: upper limit of intercritic stable eutectoid.  
C     TEA1F: lower limit of intercritic stable eutectoid.  
C     TEaTP: upper limit of intercritic metastable eutectoid.  
C     TEA1P: lower limit of intercritic metastable eutectoid.  
C     
C***  paper Valerie GERVAL and Jacques LACAZE -- ISIJ INternational, 
C     Vol. 40 (2000), No. 4, pp. 386-392
C     stable eutectoid temperature 
      TEaTF= 739.0D0+ 31.5D0*SI- 7.7D0*CU- 18.7D0*QN
     .     + 3.3D0*QO- 10.7D0*CR- 26.0D0*QI
      TEA1F= 739.0D0+ 18.4D0*SI+ 2.0D0*SI**2.0D0
     .     - 14.0D0*CU- 45.0D0*QN+ 2.0D0*QO- 24.0*CR- 27.5*QI
C     metastable eutectoid temperature
      TEaTP= 727.0D0+ 30.07D0*SICON- 1.98D0*SICON**2.0D0
     .     - 10.7D0*CUCON- 13.7D0*QNCON+ 9.3D0*QOCON+ 24.3D0*CRCON
     .     - 12.0D0*QICON
      TEA1P= 727.0D0+ 21.6D0*SICON+ 2.3D-2*SICON**2.0D0
     .     - 21.0D0*CUCON- 25.0D0*QNCON+ 8.0D0*QOCON+ 13.0D0*CRCON
     .     - 33.0D0*QICON
C     
C**** COMPUTES EQUILIBRIUM CARBON CONCENTRATION & OTHER PARAMETERS
C     
C     Nomenclature:
C     CFA: ferrite %C at the ferrite-austenite interface
C     CAF: austenite %C at the ferrite-austenite interface
C     CFG: ferrite %C at the ferrite-graphite interface
C     
C***  paper lacaze -- gerval - ISIJ  International, Vol. No 38 (1998), No. 7
C     
      CAF    = (0.1876D0-((4.112D-4)*TGAUST)+((2.26D-7)
     .     *TGAUST**2.0D0)+((0.125D0)*(SI/100.0D0)))*100.0D0
      CFA    = (6.7D-4-(5.0D-10*TGAUST**2.0D0)-
     .     2.8D-7*TGAUST+(1.2D-8*TGAUST**2.0D0-
     .     4.5D-3)*SI/100.0D0)*100.0D0
      SOBR   = TEA1F-TGAUST     ! undercooling
      DELTAWC= (SOBR*(2.9D-6-(2.8D-5*(SI/100.0D0))))*100.0D0
      CFG    = CFA-DELTAWC
C**   control on cfg value
      IF(CFG.LE.0.0D0) THEN
         IF(ICFG.EQ.0) THEN
C**   set cfg index and calculate CFA and CFG with temperature 
C**   corresponding to time t
            ICFG= 1
            CFA    = (6.7D-4-(5.0D-10*TGAUAT**2.0D0)-
     .           2.8D-7*TGAUAT+(1.2D-8*TGAUAT**2.0D0-
     .           4.5D-3)*SI/100.0D0)*100.0D0
            SOBR   = TEA1F-TGAUAT ! undercooling
            DELTAWC= (SOBR*(2.9D-6-(2.8D-5*(SI/100.0D0))))*100.0D0
            CFG    = CFA-DELTAWC
            CFGO= CFG
c$$$  CF= CFG
         ELSE
            CFG= CFGO
c$$$  CF= CFG
         ENDIF
      ENDIF
C     
C***  eutectoid metaestable equilibrium carbon concentration by 
C     Kapturkiewicz W. et al. 
      CMAC = (TGAUST  - 491.27D+0) / 365.95D+0 ! not used yet
      CMAF = (976.3D+0-TGAUST)/372.3D+0
      CMFA = (991.73D+0-TGAUST)/15498.0D+0
      CMCEM= 6.67D+0            ! carbon concentration in cementite
C     
C***  C**** difussion coeficient of carbon in austenite -- Lacaze-Gerval --  
C     ISIJ international, Vol. 38 (1998), No.7.
C***  paper lacaze -- gerval - ISIJ  International, Vol. No 38 (1998), No. 7
C***  carbon difussion coefficiente in ferrite
      TCURIE= 1043.0D0-1000.0D0*(SI/100.0D0)
      DIFCF = 2.0D-6*DEXP(-10115.0D0/(TGAUST+274.15D0))*
     .     DEXP(0.5898D0*(1.0D0+2.0D0/(TWOPIT/2.0D0)*
     .     DATAN((15629.0D0/TCURIE)-(15309.0D0/(TGAUST+274.15D0)))))
c$$$  DIFCF= 8.0D-11
C***  carbon difussion coefficiente in austenite
C     esta dividido por diez en base a la propuesta hecha por 
C     Lacaze and Gerval en el
C***  paper lacaze -- gerval - ISIJ  International, Vol. No 38 (1998), No. 7
      DIFCA = (2.343D-5*DEXP(-17767.0D0/(TGAUST+274.15D+0)))/10.0D0
C***  computes stable and metastable undercooling
C     
      GUGE   = TEaTF-TGAUST     ! eutectoid stable (ferrite) undercooling
c$$$  GUGM   = TEaTP-TGAUST     ! eutectoid metastable (pearlite) undercooling
C     metastable eutextoide temperature defined in base of 
C     Kapturkiewicz W. et al. 
      TEpER  = 731.699
      GUGM   = TEpER-TGAUST     ! eutectoid metastable (pearlite) undercooling
      ICORREC= 0                ! not used      
C     
      IF(FSOLF.EQ.0.0D0) CSP= CAF ! initial carbon in austenite
C     
C**** ESTABLISHES VARIABLES OF MODEL 12
C     
      FEUTEC = 1.0D0-FSOLA      ! eutectoid fraction
      FEUTECO= 1.0D0-FSOLAO
C     
C**** SOLVES MODEL 12
C     
      TSOCGF= 0.0D0             ! only useful to thermal jacobian matrix
      TSOCGP= 0.0D0
      TSOCGG= 0.0D0
C***  calculate volume of carbon in graphite, ferrite and pearlite 
C     respectively
      UVVGRO= UVVGR
      UVVFO = UVVFR
      UVVPO = UVVPR
C     
C**** AUSTENITE-GRAPHITE PHASE CHANGE 
C     CORRESPONDING BETWEEN EUTECTIC AND STABLE EUTECTOID TEMPERATURE
C     temperature is major of TEaTF
      IF(GUGE.LT.0.0D0.AND.FSOLA.GT.0.0D0.AND.FSOLF.EQ.0.0D+0) THEN 
C***  control to verify that cast has solidified ¿is fliq == 0.0D0?
         FCONTROL= 1.0D0-(FSOLA+FSOLG+FSOLF+FSOLP)
C***  
         IF(FCONTROL.LE.1.0D-5.AND.ICEp.EQ.0) THEN
C     
C**** Nomenclature
C     
C     CEP: max %C in aust.
C     CAG: austenite %C at the austenite-graphite interface
C***  set control
            ICEp= 1.0D+0
C***  ====================================================================
C***  calculamos CAGR entre la temperatura eutectica y la 
C     eutectoide estable
C     
C***  asigno a la variable TSP el valor correspondiente a la temperatura 
C     del limite eutectoide estable superior
            TSp= TEaTF
C     
C***  asigno a la variable CAGRO el valor de la concentracion de 
C     equilibrio de carbono correspondiente al limite eutectoide 
C     estable superior correspondiente a la temperatura TEaTF calculado 
C     con la concentracion de Si correspondiente a las primeras zonas en 
C     solidificar (a las cuales las considero iguales a las inciciales) 
            CAFO= (0.1876D0-((4.112D-4)*TSp)+((2.26D-7)
     .           *TSp**2.0D0)+((0.125D0)*(SI/100.0D0)))*100.0D0
C     
C**** for impingement and latent heat purposes
c$$$  FSOLAINI= FSOLA
C***  calculate porcentage of austenite transformed
c$$$  RFA     = FSOLA/FSOLAINI
            CPROS   = CPRO               
c$$$  write(89,*) 'cpros_1= ', cpros
            CPROSO  = CPRO
C     
C***  ====================================================================
C     
         ENDIF
C***  we calculate graphite radius growth only if cpro is greater than cag
C     and fliq == 0.0D0
         IF(FCONTROL.LE.1.0D-5) THEN
c$$$  IF(FCONTROL.EQ.0.0D0) THEN
C     
C***  asigno a la variable TEP el valor de la temperatura eutectica
            TEp= 1154.0D0 +4.0D0*SI -2.0*QN -30.0*PH
C     
C***  concentracion correspondiente al punto E prima
            CEp= 21.0D-1 -2.0D-1*SI +1.112D-2*QN +1.68D-1*PH
C     
C***  concentracion de equilibrio de carbono en la austenita en 
C     contacto con el grafito
            CAG  = ((CAFO-CEp)/(TSp-TEp))*TGAUST- ((CAFO-CEp)/(TSp-TEp)) 
     .           *TEp+ CEp
C***  thermodynamics parameters uses in graphite growth
            AKG= DIFCA*(CPROS-CAG)/ ! for graphite nodule growth
     .           ((DENSG/DENSA)*100.0D0-CAG)
C***  inicializo los valores de las fracciones y de las variables
C     empleados para el calculo del jacobiano
            FSOLGX= 0.0D0
            TSOCGX= 0.0D0
C     
            IF(JO.GT.0) THEN
               DO IJO=1,JO
C***  graphite radius growth
                  RINIGR    = RNU(IJO)
c$$$  RNU(IJO)  = RNU(IJO)+AKG*DTIMET/RNU(IJO)*RFA ! uncomment
                  RNU(IJO)  = RNU(IJO)+AKG*DTIMET/RNU(IJO)
                  RFINGR    = RNU(IJO)
C***  calculate new graphite phase fraction
                  FSOLGX  = FSOLGX+2.0D0/3.0D0*TWOPIT*DNU(IJO,1)*
     .                 RNU(IJO)**3.0D0
C***  calculo el area delos nodulos cubiertas por los granos de ferrita
C     
CCCC  write(87,*) GNF, RINIFE
C     
                  RAFGGN(IJO)= (GNF*(FNU(IJO)**2.0D0))/
     .                 (4.0D0*(RNU(IJO)**2.0D0))
               ENDDO            ! ijo=1,jo
            ENDIF               ! jo.gt.0
C     calculate new values of graphite and austenite fractions
            FSOLG= FSOLGX
C     FSOLA= 1.0D0-FSOLG
            FSOLA= 1.0D0-FSOLG-FSOLF-FSOLP
C***  calculate porcentage of austenite transformed
c$$$  RFA  = FSOLA/FSOLAINI     
C     
            UVVGRO= FSOLGO
            UVVGR = FSOLG
C     
            IF(UVVGRO.NE.0.0D0) THEN
               CPROSP= ((CPROS*(1.0D0-UVVGRO)*DENSA) + 
     .              (100.0D0*(UVVGRO-UVVGR))*DENSG)/
     .              ((1.0D0-UVVGR)*DENSA) 
C     
c$$$  write(89,*) 'cpros_2= ', cpros
            ELSE
               CPROSP= CPROS
            ENDIF
C     
CCC   write(79,*) tgaust, fsolg
CCC   write(79,*) tgaust, fsola
CCC   write(79,*) tgaust, fsolf
CCC   write(79,*) tgaust, fsolp
C     
         ELSE                   ! if fliqd .ne. 0.0D0
            FSOLG= FSOLGO
            FSOLA= FSOLAO
            UVVGR= UVVGRO
            CPROP= CPROP
            RFA  = RFAO
         ENDIF
C     
C***  ====================================================================C     
C**** STABLE PHASE CHANGE === GRAPHITE+AUSTENITE --> GRAPHITE+FERRITE
C     
      ELSEIF(GUGE.GT.0.0D0.AND.FSOLA.GT.0.0D0) THEN
C***  
         IF(ICEpp.EQ.0) THEN
            ICEpp= 1.0D+0
C**** for impingement and latent heat purposes
            FSOLAINI= FSOLA
         ENDIF
C     write(79,*) fsolg
C     write(79,*) fsola
C     write(79,*) fsolf
C     write(79,*) fsolp
C     
C***  ferrite and graphite growth during stable eutectoid intercritic
         FSOLGX= 0.0D0
         FSOLFX= 0.0D0
         TSOCGX= 0.0D0
         TSOCFX= 0.0D0
C     
c$$$  write(79,*) rfa, fsola, fsolaini
         RFA   = FSOLA/FSOLAINI ! austenite fraction percent transformed
C     
CCC   write(79,*) tgaust, rfa
C     
C***  calcula CAG for graphite growth
C     
C***  equilibrium carbon concentration in austenita in contact 
C     graphitefrom TEaTF to TEA1TFC
C     
C***  asigno a la variable TEP el valor de la temperatura eutectica
         TEp= 1154.0D0 +4.0D0*SI -2.0*QN -30.0*PH
C     
C***  concentracion correspondiente al punto E prima
         CEp  = 21.0D-1 -2.0D-1*SI +1.112D-2*QN +1.68D-1*PH
C     
C***  concentracion de equilibrio de carbono en la austenita en 
C     contacto con el grafito
         CAG  = ((CAFO-CEp)/(TSp-TEp))*TGAUST- ((CAFO-CEp)/(TSp-TEp)) 
     .        *TEp+ CEp
C***  control on cag
         IF(CAG.LE.0.0D0) CAG= CAGO
C     
C***  temperatures nomenclature
C     stables eutectoid temperature
C     TEaTF (upper limit) - TEA1F (upper limit), 
C     metastables eutectoid temperature
C     TEaTP (upper limit)  - TEA1P (lower limit)         
c$$$  IF(TGAUST.GE.TEA1F.AND.TGAUST.LE.TEaTF) THEN
CCC   
         IF(TGAUST.GE.TEA1F.AND.TGAUST.LE.TEaTF.AND.CPROS.LE.6.67D0) 
     .        THEN
C     
c$$$  write(79,*) tgaust
C     
C***  ///// FERRITE GRAINS GROWTH RATE /////
C     
C***  considero el crecimiento de los granos de ferrita si y solo si 
C     CAF is greater than CPROS
C     
            IF(CAF.GT.CPROS) THEN
C     
C***  thermodynamics parameters uses in ferrite growth
               AKF= DIFCA*(CAF-CPROS)/ ! for ferrite nodule growth
     .              (CAF-CFA)
C     
c$$$  write(83,*) AKF, rfa, caf, cpros
C     
               IF(JO.GT.0) THEN
                  DO IJO=1,JO
                     TSOCGGX= 0.0D0 ! temperature derivative of -RNU
                     TSOCGFX= 0.0D0 ! temperature derivative of -FNU
C     
C***  relation between radius (part of carbon gradient 2)
CCC   write(79,*) deltaf(ijo)
                     RRFG    = (RNU(IJO)+FNU(IJO)+DELTAF(IJO))
     .                    /((RNU(IJO)+FNU(IJO))*DELTAF(IJO))
C     
c$$$  write(80,*) RRFG
C     
                     RINIF   = FNU(IJO)
c$$$  FNU(IJO)= FNU(IJO)+RRFG*AKF*DTIMET*RFA
                     FNU(IJO)= FNU(IJO)+RRFG*AKF*DTIMET
                     RFINF   = FNU(IJO)
C***  porcentage of graphite nodule ocuped by ferrite grains
C***  relation area ferrite grains graphite nodules (RAFGGN)
                     RAFGGN(IJO)= (GNF*(FNU(IJO)**2.0D0))/
     .                    (4.0D0*(RNU(IJO)**2.0D0))
C     
C***  calculo el espesor de la capa limite delante del frente de 
C     transformacion (ferrita austenita)
C     para evitar problemas numericos debidos al por cero
                     DELTAF(IJO)= 2.0D0*DIFCA/((RFINF-RINIF)/DTIMET)
C     
C***  calculate new ferrite phase fraction
                     FSOLFX  = FSOLFX+1.0D0/3.0D0*TWOPIT*DNU(IJO,1)*
     .                    GNF*FNU(IJO)**3.0D0
C     
C***  only usefull for thermal jacobian
                     TSOCFX= TSOCFX+TSOCGFX*2.0D0*TWOPIT*DNU(IJO,1)*
     .                    FNU(IJO)**2.0D0
C     
                  ENDDO
c$$$  write(80,*) fsolfx
               ENDIF
            ELSE                ! caf.gt.cproso     
               FSOLFX= FSOLFO
c$$$  write(80,*) fsolfx
            ENDIF
C     
C***  ///// GRAPHITE NODULES GROWTH RATE /////
C     
            IF(CPROS.GT.CAG) THEN
C***  thermodynamics parameters uses in graphite growth
               AKG= DIFCA*(CPROS-CAG)/ ! for graphite nodule growth
     .              ((DENSG/DENSA)*100.0D0-CAG)
C     
               IF(JO.GT.0) THEN
                  DO IJO=1,JO
C     
C***  graphite nodule growth only if it is in contact with austenite
                     IF(RAFGGN(IJO).LT.1.0D0) THEN
                        RINIGR    = RNU(IJO)
                        RNU(IJO)  = RNU(IJO)+((AKG*DTIMET*
     .                       (1.0D0-RAFGGN(IJO)))/RNU(IJO))*RFA ! uncomment
c$$$  RNU(IJO)  = RNU(IJO)+((AKG*DTIMET*
c$$$  .                       (1.0D0-RAFGGN(IJO)))/RNU(IJO)) ! uncomment
                        RFINGR    = RNU(IJO)
C***  phase fractions and jacobian
C     
C***  calculate new graphite phase fraction
                        FSOLGX  = FSOLGX+2.0D0/3.0D0*TWOPIT*
     .                       DNU(IJO,1)*RNU(IJO)**3.0D0
C***  porcentage of graphite nodule ocuped by ferrite grains
C***  relation area of ferrite grains and graphite nodules (RAFGGN)
                        RAFGGN(IJO)= (GNF*(FNU(IJO)**2.0D0))/
     .                       (4.0D0*(RNU(IJO)**2.0D0))
C     
C***  only for thermal jacobian (not used)
                        TSOCGX= TSOCGX+TSOCGGX*2.0D0*TWOPIT*
     .                       DNU(IJO,1)*RNU(IJO)**2.0D0
                     ENDIF
                  ENDDO         ! ijo=1,jo
               ENDIF            ! jo.gt.0
            ELSE
               FSOLGX= FSOLGO
            ENDIF               ! cpros.gt.cag
C     
C***  calculate new values of ferrite and austenite fractions, and
C     calculate for thermal jacobian (not used)
            TSOCGG= TSOCGX
            TSOCGF= TSOCFX
C     
C***  calculate new graphite, ferrite and austenite volumetric 
C     fraction
C***  calculate ferrite volumetric fraction 
            FSOLF= FSOLFX
C**   control on ferrite phase fraction
            IF(FSOLF.LT.FSOLFO) FSOLF= FSOLFO
C***  calculate grphite volumetric fraction 
            FSOLG= FSOLGX
C**   control on graphite phase fraction
            IF(FSOLG.LT.FSOLGO) FSOLG= FSOLGO
C***  calculate new eutectoid fraction (previous control)
            IF((FSOLG+FSOLF+FSOLP).GT.1.0D0) THEN ! assumption valid for small dt
               REAGX= (1.0D0-FSOLP)/(FSOLG+FSOLF) ! excess factor
               FSOLG= FSOLG*REAGX ! no control on FSOLP, RNU & FNU
               FSOLF= FSOLF*REAGX
            ELSE
               FSOLEX= FSOLG+FSOLF+FSOLP
            ENDIF
C     
C     FSOLA= 1.0D0-FSOLG-FSOLF
            FSOLA= 1.0D0-FSOLG-FSOLF-FSOLP
C***  calculate porcentage of austenite transformed
            RFA  = FSOLA/FSOLAINI     
C***  calculate carbon concentration in austenite after 
C     ferrite and graphite growth 
            UVVGRO= FSOLGO
            UVVFO = FSOLFO
            UVVGR = FSOLG
            UVVF  = FSOLF
C     
            CF= CFG
C     
            CPROSP= ((CPROS*(1.0D0-UVVGRO-UVVFO)*DENSA) + 
     .           ((100.0D0*(UVVGRO-UVVGR))*DENSG)+
     .           ((CF*(UVVFO-UVVF))*DENSA))/
     .           ((1.0D0-UVVGR-UVVF)*DENSA) 
C     
C     CPROSP= ((CPROS*(1.0D0-UVVGRO-UVVFO)*DENSA) + 
C     .           ((100.0D0*(UVVGRO-UVVGR))*DENSG))/
C     .           ((1.0D0-UVVGR-UVVF)*DENSA) 
CCC   write(86,*) cprosp
C     
CCC   write(89,*) 'cpros_3= ', cprosp
C     write(82,*) fsolg
C     write(82,*) akg1
C     write(82,*) akg2
C     write(82,*) akfg
C     write(82,*) akfa
C     write(82,*) fsola
c$$$  write(82,*) fsolf
C     write(82,*) fsolp
C     write(82,'(/)')
C     
         ELSE
            CPROSP= CPROS
         ENDIF
C***  end ferrite and graphite growth for temperatures in the 
C     stable eutectoid intercritic growth 
C     
C***  start ferrite and graphite growth for temperatures lower then lower 
C     stable eutectoid intercritic temperatura 
C     
c$$$  IF(TEA1F.GT.TGAUST) THEN
         IF(TEA1F.GT.TGAUST.AND.CPROS.LE.6.67D0) THEN
C     
c$$$  write(84,*) tgaust
C     
CCC   IF(CPROS.LE.0.0D0) CPROS= CPROSO
C     
C***  calculate thermodinamics values for ferrite and graphite growth
C***  @ ferrite growth
C     
CCC   write(87,*) difcff
C     
            AKFG= DIFCF*(CFA-CFG)/(CAF-CFA)
C     
CCCC  write(82,*) 'cpros_u= ', cpros
            IF(CAF.GT.CPROS) THEN
C     
CCC   write(81,'(/)') 
CCC   write(81,*) '***********************************' 
CCC   write(81,*) fsolg
c$$$  write(81,*) fsolf
CCC   write(81,*) cpros
CCC   write(81,*) caf
CCC   write(81,*) '***********************************' 
CCC   write(81,'(/)') 
C     
               AKFA= DIFCA*(CPROS-CAF)/(CAF-CFA)
            ELSE
               AKFA= 0.0D0
            ENDIF               ! caf.ge.cproso
C***  @ graphite growth
            AKG1= DIFCF*DENSA*(CFA-CFG)/(DENSG*(100.0D0-CFG))
C     
            IF(CPROS.GT.CAG) THEN
               AKG2= DIFCA*(CPROS-CAG)/ ! for graphite nodule growth
     .              ((DENSG/DENSA)*100.0D0-CAG)
            ELSE
               AKG2= 0.0D0
            ENDIF
C     
C     write(83,*) fsolg
C     write(83,*) akg1
C     write(83,*) akg2
C     write(83,*) akfg
C     write(83,*) akfa
C     write(83,*) fsola
CCCC  write(83,*) fsolf
C     write(83,*) fsolp
C     write(83,'(/)')
C     
            IF(JO.GT.0) THEN
               DO IJO=1,JO
C     
C***  ///// FERRITE GRAINS AND GRAPHITE NODULES GROWTH RATE /////
C     
C***  relation between radius 
CCCC  DELTAF(IJO)= 5.0D-5
C     RRFGG= RNU(IJO)/((RNU(IJO)+FNU(IJO))*FNU(IJO))
C     val = (RNU(IJO)+FNU(IJO))*FNU(IJO)
C     val1= FNU(IJO)
C     val2= RNU(IJO)
C     
CCCC  write(79,*) deltaf(ijo)
C     
                  RRFGG= FNU(IJO)/((RNU(IJO)+FNU(IJO))*FNU(IJO))
                  RRFGA= ((RNU(IJO)+FNU(IJO))+DELTAF(IJO))/
     .                 (DELTAF(IJO)*(RNU(IJO)+FNU(IJO)))
CCCC  write(83,*) rrfgg
c$$$  write(83,*) rrfga
c$$$  write(83,*) DELTAF(IJO)
c$$$  write(83,'(/)') 
                  AKF  = AKFG*RRFGG - AKFA*RRFGA
C     
                  RAFGGN(IJO)= (GNF*(FNU(IJO)**2.0D0))/
     .                 (4.0D0*(RNU(IJO)**2.0D0))
C     
C***  ///// GRAPHITE NODULES GROWTH RATE /////
C     nodules are completly sourronded by ferrite
                  IF(RAFGGN(IJO).GE.1.0D+0) THEN     
                     RRGNG= (RNU(IJO)+FNU(IJO))/
     .                    (RNU(IJO)*FNU(IJO))
                     AKG  = AKG1*RRGNG
                  ELSE
C***  relation between radius 
C     
C     carbon flux from ferrite to austenite bulk
                     RRAG1= (RNU(IJO)+FNU(IJO))/
     .                    (RNU(IJO)*FNU(IJO))
C     carbon flux from austenite to graphite
                     RRAG2= 1.0D0/RNU(IJO)
                     AKG  = AKG1*RRAG1*RAFGGN(IJO) + 
     .                    AKG2*RRAG2*(1.0D0-RAFGGN(IJO))
                  ENDIF
C     
C     write(83,*) '**************************************'
C     write(83,*) fsolg
C     write(83,*) akg1
C     write(83,*) akg2
C     write(83,*) akfg
C     write(83,*) akfa
C     write(83,*) fsola
C     write(83,*) fsolf
C     write(83,*) rrfgg
C     write(83,*) val
C     write(83,*) val1
C     write(83,*) val2
C     write(83,*) rrfga
C     write(83,'(/)')
C     
c***  new graphite nodules radius
                  RINIGR  = RNU(IJO)
                  RNU(IJO)= RNU(IJO)+ AKG*DTIMET*RFA ! uncomment
c$$$  RNU(IJO)= RNU(IJO)+ AKG*DTIMET ! uncomment
                  RFINGR  = RNU(IJO)
C***  new ferrrite grain radius
                  RINIF   = FNU(IJO)
c$$$  FNU(IJO)= FNU(IJO)+ AKF*DTIMET*RFA
                  FNU(IJO)= FNU(IJO)+ AKF*DTIMET
                  RFINF   = FNU(IJO)
C     
C***  calculate new graphite phase fraction
                  FSOLGX  = FSOLGX+2.0D0/3.0D0*TWOPIT*DNU(IJO,1)*
     .                 RNU(IJO)**3.0D0
                  FSOLFX  = FSOLFX+1.0D0/3.0D0*TWOPIT*DNU(IJO,1)*
     .                 GNF*FNU(IJO)**3.0D0
C     
C***  calculo el espesor de la capa limite delante del frente de 
C     la ferrita
C     deltaf(ijo)= deltafini
                  DELTAF(IJO)= 2.0D0*DIFCA/((RFINF-RINIF)/DTIMET)
C     
c$$$  val88      = (RFINF-RINIF)/DTIMET 
c$$$  write(83,*) val88
c$$$  write(83,'(/)') 
C     
C***  porcentage of graphite nodule ocuped by ferrite grains
C***  relation area of ferrite grains and graphite nodules (RAFGGN)
                  RAFGGN(IJO)= (GNF*(FNU(IJO)**2.0D0))/
     .                 (4.0D0*(RNU(IJO)**2.0D0)) 
C     
C***  only for thermal jacobian (not used)
                  TSOCGX= TSOCGX+TSOCGGX*2.0D0*TWOPIT*DNU(IJO,1)*
     .                 RNU(IJO)**2.0D0
                  TSOCFX= TSOCFX+TSOCGFX*1.0D0*TWOPIT*4.0D0*
     .                 FNU(IJO)**2.0D0
               ENDDO
            ENDIF
C     
c$$$  WRITE(89,*) "FSOLF_10= ", FSOLFX
C     
C***  calculate new values of ferrite and austenite fractions, and
C     calculate for thermal jacobian (not used)
            TSOCGG= TSOCGX
            TSOCGF= TSOCFX
C***  calculate new graphite fraction (previus control)
            IF(FSOLGX.LT.FSOLGO) THEN
               FSOLG= FSOLGO
            ELSE
               FSOLG= FSOLGX
            ENDIF
C***  calculate new ferrite fraction (previus control)
            IF(FSOLFX.LT.FSOLFO) THEN
               FSOLF= FSOLFO
            ELSE
               FSOLF= FSOLFX
            ENDIF
C***  calculate new eutectoid fraction (previous control)
            IF((FSOLG+FSOLF+FSOLP).GT.1.0D0) THEN ! assumption valid for small dt
               REAGX= (1.0D0-FSOLP)/(FSOLG+FSOLF) ! excess factor
               FSOLG= FSOLG*REAGX ! no control on FSOLP, RNU & FNU
               FSOLF= FSOLF*REAGX
            ELSE
               FSOLEX= FSOLG+FSOLF+FSOLP
            ENDIF
C***  calculate new graphite, ferrita and austenite fraction
C     FSOLA= 1.0D0-FSOLG-FSOLF
            FSOLA= 1.0D0-FSOLG-FSOLF-FSOLP
C***  calculate porcentage of austenite transformed
            RFA  = FSOLA/FSOLAINI     
C     
C***  calculate carbon concentration in austenite after 
C     ferrite and graphite growth 
            UVVGRO= FSOLGO
            UVVFO = FSOLFO
            UVVGR = FSOLG
            UVVF  = FSOLF
C     
            val1= 1.0d0-uvvgro-uvvfo
            val2= uvvgro-uvvgr
            val3= uvvfo-uvvf
            val4= 1.0d0-uvvgr-uvvf
C     
CCCC  write(82,*) 'cpros_m= ', cproso
            CF= CFG
            CPROSP= ((CPROS*(1.0D0-UVVGRO-UVVFO)*DENSA) + 
     .           (100.0D0*(UVVGRO-UVVGR)*DENSG)+
     .           (CF*(UVVFO-UVVF)*DENSA))/
     .           ((1.0D0-UVVGR-UVVF)*DENSA)
C     CPROSP= ((CPROS*(1.0D0-UVVGRO-UVVFO)*DENSA) + 
C     .           (100.0D0*(UVVGRO-UVVGR)*DENSG))/+
C     .           ((1.0D0-UVVGR-UVVF)*DENSA)
C     
CCC   write(82,*) 'cpros_l= ', cprosp
CCC   write(82,'(/)')
C     
CCC   if(cprosp.lt.0.0D0) then
CCC   write(89,*) 'cprosp'
CCC   write(89,*) cprosp
CCC   write(89,*) val1
CCC   write(89,*) val2
CCC   write(89,*) val3    
CCC   write(89,*) 'cproso = ', cproso    
CCC   write(89,*) 'cpros_4= ', cpros    
CCC   write(89,*) cf    
CCC   write(89,*) tgaust
C     
CCC   write(89,*) fsolg
c$$$  write(89,*) fsolf
CCC   write(89,*) fsolp
CCC   write(89,*) 'uvvgro= ', uvvgro
CCC   write(89,*) 'uvvgr = ', uvvgr
C     
CCC   write(89,*) 'uvvfo= ', uvvfo
CCC   write(89,*) 'uvvf = ', uvvf
C     
CCC   write(89,'(/)')
CCC   endif
C     
C     write(81,*) fsolg
C     write(81,*) akg1
C     write(81,*) akg2
C     write(81,*) akfg
C     write(81,*) akfa
C     write(81,*) fsola
C     write(81,*) fsolf
C     write(81,*) fsolp
C     write(81,'(/)')
C     AKF   = AKFG*RRFGG - AKFA*RRFGA
         ELSE
            CPROSP= CPROS
         ENDIF                  ! tea1f.ge.tgaust
C     
C**** METASTABLE PHASE CHANGE === GRAPHITE+AUSTENITE --> GRAPHITE+PEARLITE
C     
         IF(GUGM.GT.0.0D0.AND.FSOLA.GT.0.0D0) THEN
CCC   IF(CAF.GT.CMAC.AND.FSOLA.GT.0.0D0) THEN
CCC   write(87,*) tgaust
            TSOCPNN= 0.0D0      ! used for thermal jacobian
            RFA    = FSOLA/FSOLAINI
C     
C***  relation eutectic grain surface area and pearlite colonies surface area
            RSEGEC  = 0.0D+0
C     
            IF(JP.GT.0) THEN
C***  numbers of eutectic cells per unit area (ncpa)
               NCPA= 19.0D+0
C***  determine the value of eutectic cell radius
               RADC  = DSQRT(1.0D-6/((TWOPIT/2.0D0)*NCPA))
C     
C***  determine the number of eutectic cells in RVE or in eutectic grain
C***  numbers of eutectic cells per unit volume (ncpv)
               NCPV  = (RADT/RADC)**3.0D0
C     
C***  total area corresponding to eutectic cell in all eutectic grain
               SCELLS= (TWOPIT/2.0D0)*RADC**2.0D0*NCPV
C     
C***  area corresponding to all pearlite colonies
               SPEARCOL= 0.0D+0
               DO IJP=1,JP
                  SPEARCOL= SPEARCOL + (TWOPIT/2.0D0)*DPERC(IJP,1)*
     .                 PNUC(IJP)**2.0D0
               ENDDO
               RSEGEC= SPEARCOL/SCELLS
            ENDIF
C     
C**** recalescence criterio de stop nucleacion
C     
            IF(INUCAX.EQ.1) THEN ! criterio stop nucleation recalescence
C     calculate pending average
               DTEMPT= (TGAUST-TGAUAT)/DTIMET 
               IF(JP.GT.0.AND.DTEMPT.GT.0.0D0)
     .              INDEXG=1    ! stop nucleation
               IF(INDEXG.EQ.0) THEN ! nucleation
C**** instantaneus nucleation laws
                  IF(INUCMPX.EQ.1) THEN ! ley nuestra dgc
                     JP= JP+1   ! increment pearlite colonies families
C     only for control purposses -- initial proofs --
                     IF(JP.GT.NNUM4T) ! see pointes.f
     .                    CALL RUNENDT(' ERROR IN MICROS12: JP GT 
     .                    NNUM4T ')
C     calculo pendiente promedio
                     DERTEMP    = (TGAUST-TGAUAT)/DTIMET 
                     DPERC(JP,1)= 1.0E+13*(-DERTEMP) 
                     INDEXG     = 1
C**** ley de nucleacion continua lacaze - gerval
C**** ley de nucleacion continua lacaze - gerval
                  ELSEIF(INUCMPX.EQ.2) THEN !  nucleation - lacaze-gerval -ISIJ
C     
C     IF(JP.GT.0) THEN
                     IF(DTEMPT.GT.0.0D0.OR.RSEGEC.GE.9.5D-1)
     .                    INDEXG=1
                     IF(INDEXG.EQ.0) THEN
                        JP= JP+1 ! increment pearlite colonies families
C     only for control purposses -- initial proofs --
                        IF(JP.GT.NNUM4T) ! see pointes.f
     .                       CALL RUNENDT(' ERROR IN MICROS12: JP GT 
     .                       NNUM4T ')
C     
                        DPERC(JP,1)= GNUCEP*GNUCCP*GUGM**
     .                       (GNUCEP-1.0D+0)*RFA*DTIMET
C     
c$$$  write(79,*) tgaust
                     ENDIF
C     
C**** ley de nucleacion M. R. Varma et. al.
                  ELSEIF(INUCMPX.EQ.3) THEN ! M. R. Varma nucleation
                     IF(JP.EQ.0.AND.DPERI.EQ.0) THEN
C     calculate limit number of colonies pearlites based  
C     on instantaneous nucleation laws     
                        DERTEMP= (TGAUST-TGAUAT)/DTIMET 
                        DPERI  = 3.0D+13*(-DERTEMP) 
                     ENDIF
                     JP= JP+1   ! increment pearlite colonies families
C     only for control purposses -- initial proofs --
                     IF(JP.GT.NNUM4T) ! see pointes.f
     .                    CALL RUNENDT(' ERROR IN MICROS12: JP GT 
     .                    NNUM4T ')
                     DPERC(JP,1)= 9.0D16*9.0D-5*DEXP(9.0D-5*GUGM)*
     .                    (-1.0D0*DTEMPT)*DTIMET*RFA
C     .                    (-1.0D0*DTEMPT)*DTIMET
                     TSOCPNN    = GNUCCP*GUGM**(GNUCEP-1.0D0)
     .                    *DTIMET*RFA ! temperature derivative of - DPERC
C     .                    *DTIMET ! temperature derivative of - DPERC
                  ELSE
                     CALL RUNENDT('ERROR IN MICROS12: INUCMPX 
     .                    DOESNT EXIST')
                  ENDIF         ! inucmpx .eq. 1 .or. eq. 2 .or. eq. 3
               ENDIF            ! indexg0 .eq.
            ELSE                ! inucax doesnt exist
               CALL RUNENDT('ERROR IN MICROS12: INUCAX 
     .              DOESNT EXIST')
            ENDIF               ! inucax.eq.1 or .eq. 2 or error
C     
C**** growth pearlite colonies -- metastable diagram
C     
            FSOLPX = 0.0D0
            TSOCPX = 0.0D0      ! only used for thermal jacobian
C     growth of pearlite colonies
            IF(JP.GT.0)THEN
               DO IJP=1,JP
C     
                  IF(IGROMPX.EQ.1) THEN
C     lacaze-IJCMR, 1999, 11, 431-436. -- boundary carbon difussion
                     PNUC(IJP)= PNUC(IJP)+ GROCP*GUGM**GROEP
     .                    *DTIMET 
                  ELSEIF(IGROMPX.EQ.2) THEN
C     M.R.Varma et. al., Bull. Mater. Sci., Vol. 24, No 3, 2001, pp. 305-312 
C     (formula 1) -- volumen carbon difussion
                     DIFCP    = DEXP(-125000.0D0/(8.31451D0* 
     .                    (TGAUST+274.15D0))) 
                     PNUC(IJP)= PNUC(IJP) + 2.5D-6*DIFCP*DTIMET*
     .                    GUGM**2*RFA
C     M.R.Varma et. al., Bull. Mater. Sci., Vol. 24, No 3, 2001, pp. 305-312 
C     (formula 2) boundary carbon difussion
                     DIFCP    = DEXP(-110000.0D0/(8.31451D0* 
     .                    (TGAUST+274.15D0))) 
                     PNUC(IJP)= PNUC(IJP)+ 1.6D-8*DIFCP*DTIMET*
     .                    GUGM**3*RFA
C     
                  ELSEIF(IGROMPX.EQ.3) THEN
C     
C     IF(CMAF.GT.CMAC.AND.TGAUST.LT.730.81342715) THEN
c$$$  IF(CMAF.GT.CMAC) THEN
CCC   IF(CAF.GT.CMAC) THEN
c$$$  GUGM= TGAUST-TGAUAT
C     
C**   value of interfacial energy proposed by Pandtit and Bhadeshia
CCC   SIGMAALPHATHETA = 8.0D-1 !((J/m^2))
C**   value of interfacial energy based on maximum growth criterion proposed
C     in:
C     Neural Network interlamellar spacing by C. CAPDEVILA, 
C     F G. CABALLERO and C. GARCIA DE ANDRES (ISIJ, 2005, Vol. 45, No. 2, 
C     pp. 229-237)
                     SIGMAALPHATHETA = 9.4D-1 !((J/m^2))
C**   value of interfacial energy based on maximum rate of entropy production 
C     criterion:
C     Neural Network interlamellar spacing by C. CAPDEVILA, 
C     F G. CABALLERO and C. GARCIA DE ANDRES (ISIJ, 2005, Vol. 45, No. 2, 
C     pp. 229-237)
CCC   SIGMAALPHATHETA = 6.3D-1 !((J/m^2))
C     
C**   valor del cambio de entalpia propuesto por C. Zener (PhD C.C.Montes)
                     DELTAeNTALpERLITA= 4.19D+6 ! (J/m^3)
C**   valor del cambio de entalpia proposed in: "Neural Network 
C     interlamellar spacing" by C. CAPDEVILA, F G. CABALLERO and 
C     C. GARCIA DE ANDRES (ISIJ, 2005, Vol. 45, No. 2, pp. 229-237)
C     whos cited Kramer et al.
c$$$  DELTAeNTALpERLITA= 6.09D+8 ! (J/m^3)
                     PRO              = (7.6D-8*DEXP(-10350.0/
     .                    (8.31451D0*(TGAUST+274.15))))*1.0D-6 ! K_x*DIFCB*DELTAP
C     (formula 2) boundary carbon difussion
C***  critical interlamellar spacing
CCCC  HLAMBDAC= 2.0D0*SIGMAALPHATHETA*TEA1P/
                     HLAMBDAC= (2.0D0*SIGMAALPHATHETA*TEpER)/
     .                    (DELTAeNTALpERLITA*GUGM)
CCCC  HLAMBDAC= (2.0D0*SIGMAALPHATHETA*730.813427159)/
CCC   .                       (DELTAENTALPERLITA*GUGM)
C**   value of interlamellar spacing based on maximum growth rate
                     HLAMBDA= 2.0D0*HLAMBDAC
C**   value of interlamellar spacing based on maximum entropy genreation
CCC   HLAMBDA= 3.0D0*HLAMBDAC          
C     
CCC   write(87,*) tgaust
CCC   write(85,*) TEA1F
CCC   write(85,*) TEaTP
CCC   write(85,*) TEA1P
C     
C**   value of boundary difussion coefficient (Pandtit and Bhadeshia)
                     DIFCB  = 8.51D-5*DEXP(-96851.0D0/
     .                    (8.31451D0*(TGAUST+274.15))) 
C***  values of interlamellar spacing 
                     HLAMBDAPER= HLAMBDA
                     HLAMBDACEM= HLAMBDAPER/8.0D+0
                     HLAMBDAFER= 7.0D0*HLAMBDACEM
C***  
                     PNUC(IJP)= PNUC(IJP)+ (((CMAF-CMAC)/
     .                    (CMCEM-CMFA))*
     .                    ((2.0D0*DIFCA/7.0D-1)+
     .                    (12.0D0*PRO/HLAMBDAPER))*
     .                    (HLAMBDAPER/
     .                    (HLAMBDAFER*HLAMBDACEM))*
     .                    (1.0D0-(HLAMBDAC/HLAMBDAPER)))*DTIMET*RFA
c$$$  .                    (1.0D0-(HLAMBDAC/HLAMBDAPER)))*DTIMET*RFA
C     
c$$$  PNUC(IJP)= PNUC(IJP)+ (((CAF-CMAC)/
c$$$  .                       (CMCEM-CFA))*
c$$$  .                       ((2.0D0*DIFCA/7.0D-1)+
c$$$  .                       (12.0D0*DIFCB*1.0D-4/HLAMBDAPER))
c$$$  .                       (12.0D0*PRO/HLAMBDAPER))
c$$$  .                       *(HLAMBDAPER/
c$$$  .                       (HLAMBDAFER*HLAMBDACEM))*
c$$$  .                       (1.0D0-(HLAMBDAC/HLAMBDAPER)))*DTIMET !*RFA
C     
c$$$  PNUC(IJP)= PNUC(IJP)+ (((CMAF-CMAC)/
c$$$  .                       (CMCEM-CMFA))*
c$$$  .                       ((2.0D0*DIFCA)+
c$$$  .                       (12.0D0*DIFCB*1.0D-4/HLAMBDAPER))
c$$$  .                       *(HLAMBDAPER/
c$$$  .                       (HLAMBDAFER*HLAMBDACEM))*
c$$$  .                       (1.0D0-(HLAMBDAC/HLAMBDAPER)))*DTIMET!*RFA
C     
                     val33= CMAF  - CMAC
                     val22= CMCEM - CMFA
c$$$  write(83,*) ijp, pnuc(ijp), caf, cmaf, cmac, cfa
CCC   write(83,*) val22
CCC   write(83,*) val33
CCC   write(83,*) cmaf, cmac, cmcem, cmfa, tgaust
CCC   write(83,'(/)')                 
c$$$  ELSE
c$$$  PNUC(IJP)= PNUC(IJP)
c$$$  ENDIF
                  ELSE
                     CALL RUNENDT('ERROR IN MICROS12: IGROMPX
     .                    DOESNT EXIST')
                  ENDIF         ! igrompx
C     calculate pearlite fraction
                  TSOCGPX= 0.0D0
                  FSOLPX = FSOLPX+1.0D0/3.0D0*TWOPIT*DPERC(IJP,1)*
     .                 PNUC(IJP)**3.0D0 ! pearlite fraction
C     
C                 if(fsolpx.lt.0.0D0)  write(84,*) fsolpx
C     
C                if(fsolpx.gt.1.0D0) then
C                    write(83,*) fsolpx
C                    write(83,*) ijp, pnuc(ijp)
C                    write(83,*) ijp, pnuc(ijp)
C                    write(83,*) CMFA
C                    write(83,*) CMAF
C                    write(83,*) CMAC
C                    write(83,*) CMCEM
C                    write(83,*) DIFCA
C                    write(83,*) PRO
C                    write(83,*) HLAMBDAPER
C                    write(83,*) HLAMBDAFER
C                    write(83,*) HLAMBDACEM
C                    write(83,*) HLAMBDAC
C                    write(83,'(/)')
C                 endif
C     
c$$$  write(79,*) HLAMBDAPER
c$$$  write(83,*) ijp, pnuc(ijp)
C     
CCC   write(84,*) rfa
CCC   write(84,'(/)')
C     
                  IF(IJP.EQ.JP)
     .                 TSOCPX= TSOCPX+TSOCGPX*2.0D0*TWOPIT*
     .                 DPERC(JP,1)*PNUC(IJP)**2.0D0     
C     print only for control purposses
               ENDDO            ! ijp=1,jp
            ENDIF               ! jp .gt. 0
C     calculate the new pearlite fraction (previous control)
            IF (FSOLPX.LT.0.0D0) THEN
               FSOLP = -FSOLPX
               TSOCPX= -TSOCPX
            ELSE
               FSOLP = FSOLPX
               TSOCGP= TSOCPX
            ENDIF
C***  calculate carbon concentration in austenite after 
C     ferrite and graphite growth 
            UVVGRO= FSOLGO
            UVVFO = FSOLFO
            UVVPO = FSOLPO
            UVVGR = FSOLG
            UVVF  = FSOLF
            UVVP  = FSOLP
C     
            CF= (CMFA+CMAF)/2.0D0
C     
            CF= CFG
C     
c$$$  C     
            CPROSP= ((CPROS*(1.0D0-UVVGRO-UVVFO-UVVPO)*DENSA) + 
     .           ((100.0D0*(UVVGRO-UVVGR))*DENSG)+
     .           ((CF*(UVVFO-UVVF))*DENSA)+
     .           ((CF+CMCEM)*(UVVPO-UVVP)*(8.0D0/7.0D0)*DENSA))/
     .           ((1.0D0-UVVGR-UVVF-UVVP)*DENSA)
C     
            if(cprosp.le.0.0D0) cprosp= cpros
c$$$  write(88,*) cpros, cprosp, uvvpo, uvvp, fsolf
c$$$  write(88,'(/)')
C     
         ENDIF                  ! gugm.gt.0.0.and....
C     
C**** calculamos la nueva fraccion de perlita (previo control)
         IF((FSOLG+FSOLF+FSOLP).GT.1.0D0) THEN ! assumption valid for small dt
c$$$  write(87,*) fsola
c$$$  write(87,*) fsolg
c$$$  write(87,*) fsolf
c$$$  write(87,*) fsolp
c$$$  write(87,'(/)')
            FSOLP= 1.0D0-FSOLF-FSOLG ! only FSOLP is corrected
         ENDIF
C     
         FEUTEC= FSOLG+FSOLF+FSOLP ! always
C     FSOLA = 1.0D0-FEUTEC   ! austenite fraction
         FSOLA = 1.0D0-FEUTEC   ! austenite fraction
C     control on austenite fraction
         IF(FSOLA.LE.1.0D-5) THEN
            FSOLA= 0.0D0
         ENDIF
         RFA= FSOLA/FSOLAINI    ! recalculate ferrite/austenite fraction
C     control on RFA value
c$$$  IF(RFA.LT.0.0D0) THEN
c$$$  RFA = 0.0D0
c$$$  RFAO= 0.0D0
c$$$  ENDIF
C     
C        if(fsolp.lt.0.0D0)  write(85,*) fsolp
      ENDIF                     ! guge.gt.0.0.and....
CCC   IF(FSOLF.GT.0.0D0) THEN
c$$$  esto lo hago para que la primera vez que se cumpla que rfa.ne.0.0D0 
c$$$  y consiguientemente rfao sea igual a cero, no de un calor 
c$$$  latente positivo y la temperatura en lugar de aumentar desecienda
      if(rfa.ne.0.0D0.and.rfao.eq.0.0D0) rfao= rfa
c$$$  write(86,*) rfa, rfao, fsola, fsolaini, fsolg, fsolf, fsolp
      tsoe1(iplat)= -rfa        ! f_pc at time t+dt
      tsoe2(iplat)= -rfao       ! f_pc at time t
      tsoe2(iplat)= -tsoe2(iplat)*henerf ! l*f_pc at time t
      tsoe1(iplat)= -tsoe1(iplat)*henerf ! l*f_pc at time t+dt
CCCC  ENDIF
C**** ESTABLISHES LATENT HEAT * PHASE-CHANGE FUNCTION (at time t & t+dt)
C     
c$$$  TSOE1(IPLAT)=-(FSOLF *HENERF+ FSOLP *HENERP) ! f_pc at time t+dt
c$$$  TSOE2(IPLAT)=-(FSOLFO*HENERF+ FSOLPO*HENERP) ! f_pc at time t
C     
c$$$  ENDIF                     ! guge.gt.0.0.and....
C     
c$$$  write(85,*) cprosp
C     
C**** DEFINES MICROSTRUCTURAL-DEPENDENT "MACROSCOPICAL" PROPERTIES
C     
      IF(IKMICX.EQ.1.OR.IKMICX.EQ.2.OR.IKMICX.EQ.3) THEN
         BASKK(1)= BASKS        ! solid conductivity
      ENDIF
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
 1    ICACOT= 0
      VELTOT= 1.D-10
      VELA1T=DABS(DTEMPT)
      IF(VELA1T.GT.VELTOT) THEN
         IF(IITERT.GT.0)
     .        TSOC1(IPLAT)= (TSOE1(IPLAT)-TSOE2(IPLAT))/
     .        (DTEMPT*DTIMET)
      ENDIF
      GO TO 10
C     
 2    ICACOT= 1
      VELTOT= 1.D-10
      VELA1T=DABS(DTEMPT)
      IF(VELA1T.GT.VELTOT) THEN
         IF(IITERT.GT.0)
     .        TSOC1(IPLAT)= (TSOE1(IPLAT)-TSOE2(IPLAT))/
     .        (DTEMPT*DTIMET)
      ENDIF
      GO TO 10
C     
 3    ICACOT=0
      TSOC1(IPLAT)= TSOCGF+TSOCGP+TSOCGG
      TSOC1(IPLAT)= TSOC1(IPLAT)*HENERF
      GO TO 10
C     
 4    ICACOT=1
      TSOC1(IPLAT)= TSOCGF+TSOCGP+TSOCGG
      TSOC1(IPLAT)= TSOC1(IPLAT)*HENERF
      GO TO 10
C     
 5    ICACOT=1
      TSOC1(IPLAT)= 0.0D0
      GO TO 10
C     
 10   CONTINUE
C     
C**** TRANSFER MICROSCOPICAL VARIABLES TO "ALPHAM"
C     
      ALPHAM(IN+1)= FSOLF       ! ferrite fraction
      ALPHAM(IN+2)= FSOLP       ! pearlite fraction
C     carbon content in austenite in contact with graphite
      ALPHAM(IN+3)= CAG        
C     carbon content in austenite (mass equilibrium) 
      ALPHAM(IN+4)= CPROSP      
C     
      IF(JO.GT.0) THEN
         DO INU=1,JO            ! NNUM4T different ferrite grain radius
C     boundary layer at graphtie/austenita interface
            ALPHAM(IN+5+INU-1)             = DELTAG(INU) 
            ALPHAM(IN+6+NNUM4T-1+INU-1)    = DELTAF(INU) 
            ALPHAM(IN+7+2*(NNUM4T-1)+INU-1)= RAFGGN(INU)
            ALPHAM(IN+8+3*(NNUM4T-1)+INU-1)= FNU(INU) 
c$$$  WRITE(80,*) FNU(INU)
         ENDDO        
      ENDIF
C     
      IF(JP.GT.0) THEN
         DO JNU=1,JP
            ALPHAM(IN+9+4*(NNUM4T-1)+JNU-1) = DPERC(JNU,1) ! pear. gr dens.
            ALPHAM(IN+10+5*(NNUM4T-1)+JNU-1)= PNUC(JNU) ! pear. gr. rad. (Rp)
c$$$  WRITE(80,*) PNUC(JNU)
c$$$  WRITE(81,*) jnu, DPERC(JNU,1)
         ENDDO
      ENDIF
C     
C**** Boeri's model
      IF(IPLX2.EQ.4) THEN       ! IPCMO=4 (Boeri's model)
         ALPHAM(IPLX3)= FSOLA
         ALPHAM(IPLX5)= FSOLG
C     
         IF(JO.GT.0) THEN
            DO INU=1,JO
               ALPHAM(IPLX6+INU-1)= DNU(INU,1)
               ALPHAM(IPLX7+INU-1)= RNU(INU)
            ENDDO
         ENDIF
      ENDIF
C     
C**** SU' model
      IF(IPLX2.EQ.5) THEN       ! IPCMO=5  (SU' model)
         ALPHAM(IPLX3) = FSOLA
         ALPHAM(IPLX5) = FSOLG        
         IF(JO.GT.0) THEN
            DO INU=1,JO
               ALPHAM(IPLX6+INU-1)= DNU(INU,1)
               ALPHAM(IPLX7+INU-1)= RNU(INU)
            ENDDO
         ENDIF
      ENDIF
C     
C**** dgc's model
      IF(IPLX2.EQ.9) THEN       ! IPCMO=9  (dgc's model)
         ALPHAM(IPLX3 )= FSOLA
         ALPHAM(IPLX5 )= FSOLG
CCCC  ALPHAM(IPLX23)= CPROSP
CCCC  ALPHAM(IPLX15)= CUCON  ! copper content
CCCC  ALPHAM(IPLX16)= QNCON  ! manganese content
CCCC  ALPHAM(IPLX4) = RADT     
C     
         IF(JO.GT.0) THEN
            DO INU=1,JO
               ALPHAM(IPLX6+INU-1) = DNU(INU,1)
               ALPHAM(IPLX10+INU-1)= RNU(INU)
c$$$  WRITE(79,*) RNU(INU)
            ENDDO
         ENDIF
      ENDIF
C     
cccc  ALPHAM(IPLX13)            = SICONX ! silicon content
      ALPHAM(IN+11+6*(NNUM4T-1))= CSP ! initial carbon content in austenite
      ALPHAM(IN+12+6*(NNUM4T-1))= FSOLAINI ! initial fraction of austenite
      ALPHAM(IN+13+6*(NNUM4T-1))= RFA ! initial relation ferrite/austenite
C     initial pearlites colonies densities (instantaneus nucleation)
CCC   ALPHAM(IN+14+6*(NNUM4T-1))= DPERI     
      ALPHAM(IN+14+6*(NNUM4T-1))= CFGO     
      ALPHAM(IN+15+6*(NNUM4T-1))= FLOAT(JP) ! pearlite nucleation index
      ALPHAM(IN+16+6*(NNUM4T-1))= FLOAT(ICEp) ! index to set caotic
      ALPHAM(IN+17+6*(NNUM4T-1))= CAFO ! C content for final sol.
      ALPHAM(IN+18+6*(NNUM4T-1))= TSp ! temperature for final solidification
C     ALPHAM(IN+19+6*(NNUM4T-1))= UVVGR ! graphite growth rate
      ALPHAM(IN+19+6*(NNUM4T-1))= ICEpp
      ALPHAM(IN+20+6*(NNUM4T-1))= FLOAT(ICFG) ! index to set CFG
C     
      IF(INUCAX.EQ.1)
     .     ALPHAM(IN+21+6*(NNUM4T-1))= 
     .     FLOAT(INDEXG)        ! index: nucleation up to Trecal
      IF(INUCAX.EQ.2)
     .     ALPHAM(IN+21+6*(NNUM4T-1))= 
     .     TINDEXG              ! index: nucleation up to Tmin
C     
C**** INCREMENTS "ALPHAM" INDEX
C     
      INUPC= INUPC+4+6*NNUM4T+11 ! less than NBASES; see pointes.f
C     
      RETURN
      END
