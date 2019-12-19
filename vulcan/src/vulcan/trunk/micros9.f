      SUBROUTINE MICROS9(TGAUAT,TGAUST,TGAUIT,DTEMPT,
     .     BASMM,BASCC, BASKK,
     .     TSOE2,TSOE1,TSOC1,
     .     IPLAT,
     .     ALPHAM,INUPC)  
C***********************************************************************
C     
C**** THIS ROUTINE EVALUATES THE PHASE-CHANGE FUNCTION ACCORDING TO THE
C     MICROSTRUCTURAL MODEL NUMBER 9 (IPCMO=9)
C     
C     S.G. CAST IRON MICROSTRUCTURAL MODEL: DARDATI'S MODEL
C     nuclea grafito antes y con las formulas de crecimiento de Omar y 
C     con diferenciacion de radios de grafito (is it old? who is Omar?)
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
C     ALPHAM(IN+1) = liquid fraction - FLIQD
C     ALPHAM(IN+2) = austenite fraction (solid)- FSOLA
C     ALPHAM(IN+3) = graphite fraction (solid)- FSOLG
c     ALPHAM(IN+4) = silicon content (current)-SICON
C     ALPHAM(IN+5) = total grain radius - RADT
C     ALPHAM(IN+6) = radio Rg - RADG
C     ALPHAM(IN+7) = Cl/a resultante (liquido interdend.)- CLA
C     ALPHAM(IN+8) = veloc. de crec. de la punta de la dendrita - VGRA
C     ALPHAM(IN+9) = Cpro resultante (liquido intergran.)-CPRO
C     ALPHAM(IN+10)= RN - RADN
C     
C     ALPHAM(IN+11)= RC   - RADC -- cellular radius (see Rivera's Thesis)
C     ALPHAM(IN+12)= ISDASGR - i model of sdas 
C     ALPHAM(IN+13)= total time of solidification
C     ALPHAM(IN+14)= SDAS - secondary dendrite arm spacing
C     ALPHAM(IN+15)= phosporus content (current) - PHCON
C     ALPHAM(IN+16)= cooper content (current) - CUCON
C     ALPHAM(IN+17)= manganese content (current) - QNCON
C     ALPHAM(IN+18)= magnesium (magnesio) content (current) - HGCON
C     ALPHAM(IN+19)= niobium content (current) - SBCON
C     ALPHAM(IN+20)= tin (estano) content (current) - SNCON
C     ALPHAM(IN+21)= chromium  (cromo) content (current) - CRCON
C     ALPHAM(IN+22)= molybdenum (molibdeno) content (current) - QOCON
C     ALPHAM(IN+23)= nickel (nikel) content (current) - QICON
C     ALPHAM(IN+24)= not used yet
C     ALPHAM(IN+25)= not used yet
C     
C     NNUM4T= number of different grain densities and radii - JO maximo
C     
C     ALPHAM(IN+26:IN+26+(NNUM4T-1))= graphite grain density (ggd) zone 1
C     ALPHAM(IN+27+  (NNUM4T-1):IN+27+2*(NNUM4T-1))= ggd zone 2 - DNU
C     ALPHAM(IN+28+2*(NNUM4T-1):IN+28+3*(NNUM4T-1))= ggd zone 3 -DNU
C     ALPHAM(IN+29+3*(NNUM4T-1):IN+29+4*(NNUM4T-1))= gg radius - RNU
C     ALPHAM(IN+30+4*(NNUM4T-1):IN+30+5*(NNUM4T-1))= gg radius - RNUZ1
C     
C     Auxiliar microstructural variables (not printed; see pointes.f):
C     ALPHAM(IN+31+5*(NNUM4T-1))= graphite nucleation index (gni) -jo
C     ALPHAM(IN+32+5*(NNUM4T-1))= for recalescence -INDEXG
C     ALPHAM(IN+33+5*(NNUM4T-1))= for T_min -SINDEXG
C     ALPHAM(IN+34+5*(NNUM4T-1))= growth index. If RADG=RADT INDEXR=1
C     ALPHAM(IN+35+5*(NNUM4T-1))= 0.0D0 ! local final time of solidification
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
      PI    = 4.0D0*DATAN(1.0D0)
C     
      IN   = INUPC
c     
cccc  write(83,*) ' inupc_l= ', inupc
cccc  write(83,*) ' in_l   = ', in
c     
      FLIQD= ALPHAM(IN+1)
C     
C**** CHECKS LIQUID FRACTION
C     
      IF(FLIQD.LE.0.0D0) THEN   ! gets out if liquid fraction is zero
         INUPC= INUPC+26+5*NNUM4T+5
         RETURN
      ENDIF
C     
      FRGZ1 = 0.0D0             ! fraccion grafito zona 1
      FRGZ2 = 0.0D0             ! fraccion grafito zona 2
      FRGZ3 = 0.0D0             ! fraccion grafito zona 3
      FRAZ1 = 0.0D0             ! fraccion austenita zona 1
      FRAZ2 = 0.0D0             ! fraccion austenita zona 2
      FRAZ3 = 0.0D0             ! fraccion austenita zona 3
      RADN  = 0.0D0             !
      FSOLID= 0.0D0             ! fraccion solida
C     
      DTIMETO= DTIMET           !
C     
C**** TRANSFERS "VPLAT" TO MICROSCOPICAL VARIABLES
C     
      TEINF= VPLAT(IPLAT,1)     ! not used
      TESUP= VPLAT(IPLAT,2)     ! not used
      HENER= VPLAT(IPLAT,3)     ! Latent heat
      DELEE= TESUP-TEINF
C**   carbon
      CA     = VPLAT(IPLAT, 6)
      AKCAL  = VPLAT(IPLAT, 7)
      AKCAS  = VPLAT(IPLAT, 8)
      ICALSMX= INT(VPLAT(IPLAT, 9))
      ICASSMX= INT(VPLAT(IPLAT,10))
C**   silicon
      SI     = VPLAT(IPLAT,11)
      AKSIL  = VPLAT(IPLAT,12)
      AKSIS  = VPLAT(IPLAT,13)
      ISILSMX= INT(VPLAT(IPLAT,14))
      ISISSMX= INT(VPLAT(IPLAT,15))
C**   phosporus
      PH     = VPLAT(IPLAT,16)
      AKPHL  = VPLAT(IPLAT,17)
      AKPHS  = VPLAT(IPLAT,18)
      IPHLSMX= INT(VPLAT(IPLAT,19))
      IPHSSMX= INT(VPLAT(IPLAT,20))
C**   copper
      CU     = VPLAT(IPLAT,21)
      AKCUL  = VPLAT(IPLAT,22)
      AKCUS  = VPLAT(IPLAT,23)
      ICULSMX= INT(VPLAT(IPLAT,24))
      ICUSSMX= INT(VPLAT(IPLAT,25))
C**   manganese
      QN     = VPLAT(IPLAT,26)
      AKQNL  = VPLAT(IPLAT,27)
      AKQNS  = VPLAT(IPLAT,28)
      IQNLSMX= INT(VPLAT(IPLAT,29))
      IQNSSMX= INT(VPLAT(IPLAT,30))
C**   magnesium
      QG     = VPLAT(IPLAT,31)
      AKQGL  = VPLAT(IPLAT,32)
      AKQGS  = VPLAT(IPLAT,33)
      IQGLSMX= INT(VPLAT(IPLAT,34))
      IQGSSMX= INT(VPLAT(IPLAT,35))
C**   niobium
      QB     = VPLAT(IPLAT,36)
      AKQBL  = VPLAT(IPLAT,37)
      AKQBS  = VPLAT(IPLAT,38)
      IQBLSMX= INT(VPLAT(IPLAT,39))
      IQBSSMX= INT(VPLAT(IPLAT,40))
C**   tin
      SN     = VPLAT(IPLAT,41)
      AKSNL  = VPLAT(IPLAT,42)
      AKSNS  = VPLAT(IPLAT,43)
      ISNLSMX= INT(VPLAT(IPLAT,44))
      ISNSSMX= INT(VPLAT(IPLAT,45))
C**   chromium
      CR     = VPLAT(IPLAT,46)
      AKCRL  = VPLAT(IPLAT,47)
      AKCRS  = VPLAT(IPLAT,48)
      ICRLSMX= INT(VPLAT(IPLAT,49))
      ICRSSMX= INT(VPLAT(IPLAT,50))
C**   molybdenum
      QO     = VPLAT(IPLAT,51)
      AKQOL  = VPLAT(IPLAT,52)
      AKQOS  = VPLAT(IPLAT,53)
      IQOLSMX= INT(VPLAT(IPLAT,54))
      IQOSSMX= INT(VPLAT(IPLAT,55))
C**   nickel
      QI     = VPLAT(IPLAT,56)
      AKQIL  = VPLAT(IPLAT,57)
      AKQIS  = VPLAT(IPLAT,58)
      IQILSMX= INT(VPLAT(IPLAT,59))
      IQISSMX= INT(VPLAT(IPLAT,60))
C     graphite nucleation model index
C     Su, Boeri, Rappaz or Stefanescu nucl. models
      INUCMX= INT(VPLAT(IPLAT,61))
      ANUCA = VPLAT(IPLAT,62)   ! austenite nucleation coeff. A
      ANUCB = VPLAT(IPLAT,63)   ! graphite nucleation coeff. b
      ANUCC = VPLAT(IPLAT,64)   ! graphite nucleation coeff. c
C     nucleation arrest criterion index
      INUCAX= INT(VPLAT(IPLAT,65))
C     graphite growth model index -- Dardati & Zener
      IGRGMX= INT(VPLAT(IPLAT,66))
      DIFCL = VPLAT(IPLAT,67)   ! diffusion coef. of C in liquid
      RNODA = VPLAT(IPLAT,68)   ! initial radius of austenite
      RNODO = VPLAT(IPLAT,69)   ! initial radius of graphite nodule
      AUSGR = VPLAT(IPLAT,70)   ! auste. shell radius/graphite radius
      DENSA = VPLAT(IPLAT,71)   ! austenite density
      DENSG = VPLAT(IPLAT,72)   ! graphite density
C     index for primary dendrite tip grpwth -- Rappaz & Stefanescu
      IGRAMX= INT(VPLAT(IPLAT,73))
      DIFCA = VPLAT(IPLAT,74)   ! diffusion coef. of C in austenite
      AMLQ  = VPLAT(IPLAT,75)   ! liquidus pending line
      AKCA  = VPLAT(IPLAT,76)   ! carbon partition coefficient
      COGTH = VPLAT(IPLAT,77)   ! Gibb-Thompson coeficcient 
C     
      IDELRN = INT(VPLAT(IPLAT,78)) ! index for growth zone's 2 radius
      ISDASGR= INT(VPLAT(IPLAT,79)) ! index for SDAS growth 
CCCC  write(79,*) ' isdasgr= ', isdasgr
      IMICOUP= INT(VPLAT(IPLAT,80)) ! model to couple micro and macro
C     
      IKMICX = INT(VPLAT(IPLAT,81)) ! index for micro.-dep. conductivity
      IKAUX= 0    
      IF(IKMICX.EQ.1) THEN
         IF(IMICOUP.EQ.1.OR.IMICOUP.EQ.2.OR.IMICOUP.EQ.3) THEN
            IKAUX= 3
            VPLAT(IPLAT,82)= BASKS ! solid conductivity
            VPLAT(IPLAT,83)= BASKM ! mushy conductivity
            VPLAT(IPLAT,84)= BASKL ! liquid conductivity
         ENDIF
      ENDIF
C     
      IFPCDT= INT(VPLAT(IPLAT,82+IKAUX)) ! index for temp. derivative
      IAFLOJ= INT(VPLAT(IPLAT,83+IKAUX)) ! index for fraction correction
C     
      IV= 83+IKAUX              ! number of VPLAT defined in input data
C     
C**** TRANSFERS (REST OF) "ALPHAM" TO MICROSCOPICAL VARIABLES
C     
      FSOLA = ALPHAM(IN+2)
      FSOLG = ALPHAM(IN+3)
      FSOLID= FSOLA+FSOLG
C     
      SICON = ALPHAM(IN+4)
C     
cccc  write(98,*) 'si_l= ', sicon
      RADT  = ALPHAM(IN+5)
      RADG  = ALPHAM(IN+6)
      CLA   = ALPHAM(IN+7)
      VGRA  = ALPHAM(IN+8)
      CPRO  = ALPHAM(IN+9)
      RADN  = ALPHAM(IN+10)
C***  new variables
      RADC     = ALPHAM(IN+11)  ! RADC -- cellular radius (see Rivera's Thesis)
      PDAS     = ALPHAM(IN+12)  ! PDAS - - primary dendrite arm spacing
      totalTime= ALPHAM(IN+13)  ! total time from pouring
      SDAS     = ALPHAM(IN+14)  ! SDAS - secondary dendrite arm spacing
C     SDAS - secondary dendrite arm spacing -- initial for every step
      SDASINI  = ALPHAM(IN+15)    
      SDASINI  = sdas
      PHCON    = ALPHAM(IN+16)  ! phosporus content (current) - PHCON
      CUCON    = ALPHAM(IN+17)  ! cooper content (current) - CUCON
      QNCON    = ALPHAM(IN+18)  ! manganese content (current) - QNCON
      QGCON    = ALPHAM(IN+19)  ! magnesium (magnesio) content (current)- HGCON
      QBCON    = ALPHAM(IN+20)  ! niobium content (current) - NBCON
      SNCON    = ALPHAM(IN+21)  ! tin (estano) content (current)- SNCON
      CRCON    = ALPHAM(IN+22)  ! chromium  (cromo) content (current)- CRCON
      QOCON    = ALPHAM(IN+23)  ! molybdenum (molibdeno) content (current)- QOCON
      QICON    = ALPHAM(IN+24)  ! nickel (nikel) content (current)- QICON
      NUY1     = ALPHAM(IN+25)  ! not used yet 1
      NUY2     = ALPHAM(IN+26)  ! not used yet 2
C     
      JO= INT(ALPHAM(IN+32+5*(NNUM4T-1))) ! graphite nucl. index - jo
      IF(JO.GT.0) THEN
         DO INU=1,JO
            DNU(INU,1)= ALPHAM(IN+27+             INU-1 )
            DNU(INU,2)= ALPHAM(IN+28+   NNUM4T-1 +INU-1 )
            DNU(INU,3)= ALPHAM(IN+29+2*(NNUM4T-1)+INU-1 )
            RNU(INU)  = ALPHAM(IN+30+3*(NNUM4T-1)+INU-1)
            RNUZ1(INU)= ALPHAM(IN+31+4*(NNUM4T-1)+INU-1 )
         ENDDO
      ENDIF
C**** nucleation index -- recalescence or T_min --
      INDEXG = INT(ALPHAM(IN+33+5*(NNUM4T-1))) ! nucl. solo hasta recal.
      SINDEXG= ALPHAM(IN+34+5*(NNUM4T-1)) ! nucl. solo hasta T_min.
      INDEXR = INT(ALPHAM(IN+35+5*(NNUM4T-1))) ! grap. Rg growth index INDEXR=1
C***  new variable
      TF     = ALPHAM(IN+36+5*(NNUM4T-1)) ! local final time of solidification
C     
C**** LOAD LAST CONVERGED VALUES
C     
      FLIQDO= FLIQD             ! last liquid's fraction 
      FSOLAO= FSOLA             ! last austenite's fraction 
      FSOLGO= FSOLG             ! last graphite's fraction 
      RADGO = RADG              ! last graphite's radius 
      RADNO = RADN              ! last austenite's radius 
      CLAO  = CLA               ! last carbon concentration in liquid 
      VGRAO = VGRA              ! last austenite's tip growth 
C     
C**** COMPUTES EQUILIBRIUM EUTECTIC TEMPERATURE & OTHER PARAMETERS
C     
C     Nomenclature:
C     TEG:  equilibrium graphite eutectic temperature
C     CLA:  liquid %C at the liquid-austenite interface
C     CAL:  asutenite %C at the liquid-austenite interface
C     CLG:  liquid %C at the liquid-graphite interface
C     AKCA: coeficiente de particion del carbono (k)
C     
CCCCC TEG = 1154.6D0+6.5D0*SICON
c     
cccc  write(94,*) ' sicon_a*= ', sicon
cccc  sicon= si
C     
C***  TEG segun paper C03
C     modeling stable and metastable eutectic transformation of SGI
C     ISIJ INternational Vol. 41 (2001), No. 9, pp. 986-991
      TEG= 1154.0D0 +4.0D0*SICON -2.0*QNCON - 30.0*PHCON
cccc  write(94,*) ' teg*= ', teg
C     
      CLA= (1569.0D0-TGAUST- 24.32D0* SICON)/ (-1.0D0*AMLQ)
cccc  write(94,*) ' cla*   = ', cla
cccc  write(94,*) ' sicon_d*= ', si
C     C      write(94,*) ' qncon*   = ', qncon
C     C      write(94,*) ' phcon*   = ', phcon
      CAL= (1528.4D0-TGAUST- 32.00D0* SICON)/177.9D0
      CLG= ( 503.2D0+TGAUST-129.70D0* SICON)/389.1D0
C***  coeficiente de particion del C - no puedo introducirlo en el 
C     archivo *-s.dat debido a que con las cuentas del balance de 
C     masa el valor de CLA va modificandose.
C     
      AKCA= CAL/CLA             
C     
C**** ESTABLISHES VARIABLES OF MODEL 9
C     
      FSOLIDO= 1.0D0-FLIQDO
C     
C**** SOLVES MODEL 9
C     
      TSOCGN= 0.0D0             ! only useful to thermal jacobian matrix
      TSOCGG= 0.0D0
      TSOCGA= 0.0D0
C     
C**** GRAPHITE EUTECTIC SOLIDIFICATION
C     
      GUG    = TEG-TGAUST       ! graphite eutectic undercooling
      ICORREC= 0                ! not used
C
      totalTime=  totalTime+dtimet
      IF(GUG.GT.0.0D0.AND.FSOLG.LT.1.0D0.AND.FSOLID.LT.1.0D0) THEN
C     
C     ***********  ///////////////////////////////  ************
C     ***********    SOLIDIFIED FOR FIRST TIME      ************    
C     ***********  ///////////////////////////////  ************           
C     
         IF(FSOLA.EQ.0.0D0.AND.FSOLG.EQ.0.0D0) THEN ! solidified for first time
            IF(DTEMPT.GT.0.0D0) THEN ! control
               RADT= RADT+DTIMET ! guardo el t en radt para el proximo paso
               GOTO 110         ! no solidifica nada y sale de la subrutina
            ENDIF
C     
            IF(RADT.NE.0.0D0) THEN
               DERTEMP= -(GUG/(RADT+DTIMET)) ! calc average pending
            ELSE
               DTIMET = GUG/(TGAUAT-TGAUST)*DTIMETO ! cal. dt entre Te-Tgaus
               DERTEMP= DTEMPT  ! uso la pendiente que calcula el programa
            ENDIF
            SINDEXG= TGAUST
C     
C**** AUSTENITE NUCLEATION
C     
            DGRA= ANUCA*(-DERTEMP)
            IF(DGRA.EQ.0.0D0) THEN ! control
               RADT= RADT+DTIMET
               GOTO 110
            ENDIF
C***  almaceno el primer tiempo de solidificacion
            tf= tf + dtimet
C     
C**** GRAPHITE NUCLEATION
C     
            IF(INUCMX.EQ.1) THEN ! Su's nucleation
               JO=1
               DNU(1,1)= 0.0D0
               DNU(1,2)= 0.0D0
               DNU(1,3)= ANUCB*GUG**ANUCC*FLIQD*DTIMET 
            ELSEIF(INUCMX.EQ.2) THEN ! Boeri's nucleation law  (DEFAULT)
               JO= 1
               DNU(1,1)= 0.0D0
               DNU(1,2)= 0.0D0
               DNU(1,3)= ANUCB*GUG*DEXP(-ANUCC/GUG)*FLIQD*DTIMET
            ELSEIF(INUCMX.EQ.3) THEN ! Rappaz's nucleation law
               JO=1
               DNU(1,1)=0.0D0  
               DNU(1,2)=0.0D0 
               DNU(1,3)=ANUCB*DEXP(-(GUG-DIFCA)**2.0D0/ANUCC)*GUG
            ELSEIF(INUCMX.EQ.4) THEN ! Stefanescu's nucleation law 
               JO=1
               DNU(1,1)=0.0D0
               DNU(1,2)=0.0D0
               DNU(1,3)=2.0D0*ANUCB*ANUCC*DEXP(-ANUCC/GUG**2.0D0)
     .              /GUG**2.0D0*DTIMET
            ELSE
               CALL RUNENDT('ERROR IN MICROS9: INUCMX DOES NOT EXIST')
            ENDIF
C     initial radius of graphite - inicialmente todos los nodulos
C     tienen el radio inicial dado en el archivo *-s.dat
            RNU(1)  = RNODO
            RNUZ1(1)= RNODO
C     initials fractions of graphite per zone
            UVVZ1= 0.0D0
            UVVZ2= 0.0D0
            UVVZ3= 4.0D0/3.0D0*PI*DNU(1,3)*(RNU(1))**3.0D0
C     
            CPRO=(CA*DENSA-100.0D0*UVVZ3*DENSG)/((1.0D0-UVVZ3)*DENSA)
C     control on CLA
            IF(CLA.LT.CPRO) CLA= CPRO ! control (CLA can't .LE. CPRO)
C     
C**** AUSTENITE GROWING
C     
            IF(IGRAMX.EQ.1) THEN
               VGRA=DIFCL*AMLQ*(CLA-CPRO)**2.0D0 ! Rappaz (DEFAULT)
     .              /(PI**2.0D0*COGTH*(AKCA-1.0D0)*CA)
            ELSEIF(IGRAMX.EQ.2) THEN
               VGRA= AKCA*DIFCL*AMLQ*(CLA-CPRO)**2.0D0 ! Stefanescu
     .              /(2.0D0*PI**2.0D0*COGTH*(AKCA-1)*CA)
            ELSE
               CALL RUNENDT('ERROR IN MICROS9: IGRAMX DOES NOT EXIST')
            ENDIF
C     
            IF(VGRA.EQ.0.0D0) THEN ! control
C     si entra vuelve a solidificar por primera vez
               JO      = 0
               DNU(1,1)= 0.0D0 
               DNU(1,2)= 0.0D0 
               DNU(1,3)= 0.0D0 
               RADT= RADT+DTIMET
               GOTO 110      
            ENDIF
C**** calculo unit cell  radius
            RADT= (3.0D0/(4.0D0*PI*DGRA))**(1.0D0/3.0D0)
C     calculo espesor capa limite
            DELT= 2.0D0*DIFCL/VGRA
C     gradiente de concentracion
            GRARG= 2.0D0*(CLA-CPRO)/DELT
C**** CALCULO DE Rg Y Rn
            DELRG= VGRA*DTIMET
            RADG=  DELRG        ! because is the first time
            IF(RADG.GT.RADT) GOTO 20 ! control - if Rg is greater than Rt
            DELTC= CLA-CPRO
C     
C**** DELRN GROWING
C     
            IF(IDELRN.EQ.1) THEN ! calc with derivadas
               DELRN= ((DIFCL*(1.0D0-UVVZ3)*RADG*RADG*GRARG*DTIMET+
     .              DELTC*(1.0D0-UVVZ3)*RADG**3.0D0/3.0D0)/
     .              ((1.0D0-AKCA)*CLA))**(1.0D0/3.0D0)
            ELSEIF(IDELRN.EQ.2) THEN ! calc with flux
               DELRN= ((3.0D0*DIFCL*RADG*RADG*GRARG*DTIMET+DELTC*
     .              RADG**3.0D0)/((1.0D0-AKCA)*
     .              CLA))**(1.0D0/3.0D0)
            ELSEIF(IDELRN.EQ.3) THEN ! calc sin flujo (DEFAULT)
               DELRN= ((DELTC*RADG**3.0D0)/
     .              ((1.0D0-AKCA)*CLA))**(1.0D0/3.0D0)
            ELSE
               CALL RUNENDT('ERROR IN MICROS9: IDELRN DOES NOT EXIST')
            ENDIF
C     first time DELRN=RN
            RADN= DELRN
C     having Rn - Rg - y Rt -- calculate zone's fractions
            FRAZ1= (RADN/RADT)**3.0D0
            FRAZ2= (RADG**3.0D0-RADN**3.0D0)/RADT**3.0D0
            FRAZ3= (RADT**3.0D0-RADG**3.0D0)/RADT**3.0D0 ! == 1-fraz1-fraz2
C     modificate concentration
            CPRO=(CPRO-CAL*FRAZ1-CLA*FRAZ2)/FRAZ3
            CINF= CA            ! initial concentration of C
            IF(CLA.LT.CPRO) CLA= CPRO ! control
C     initial redistribution  of graphite densities
            DNU(1,1)= DNU(1,3)* FRAZ1
            DNU(1,2)= DNU(1,3)* FRAZ2
            DNU(1,3)= DNU(1,3)* FRAZ3
C     
C     ***********  ///////////////////////////////  ************
C     ***********   NOT SOLIDIFIED FOR FIRST TIME   ************    
C     ***********  ///////////////////////////////  ************           
C     
         ELSE     
C***  aumento el tiempo de solidificacion
            tf= tf + dtimet
C     
            IF(INDEXR.EQ.1) GOTO 50 ! no hay mas zona 3
C     calculate graphite fractions
            FRGZ1= 0.0D0
            FRGZ2= 0.0D0
            FRGZ3= 0.0D0
            DO IJO=1,JO                 
               FRGZ1= FRGZ1+DNU(IJO,1)*RNUZ1(IJO)**3.0D0
               FRGZ2= FRGZ2+DNU(IJO,2)*RNU(IJO)**3.0D0
               FRGZ3= FRGZ3+DNU(IJO,3)*RNU(IJO)**3.0D0
            ENDDO
            FRGZ1= FRGZ1*4.0D0/3.0D0*PI
            FRGZ2= FRGZ2*4.0D0/3.0D0*PI
            FRGZ3= FRGZ3*4.0D0/3.0D0*PI
C     
            FRAZ1= FSOLA+FRGZ1
C     control on zone1's faction  (no lo entiendo, ¿que hace aca?)
            IF(FRAZ1.LT.0.0D0) THEN 
               RADN=  0.0D0
               FRAZ1= 0.0D0
            ELSE
               RADN= ALPHAM(IN+10)
            ENDIF
C     having Rn - Rg - y Rt -- calculate zone's fractions
            FRAZ2= (RADG**3.0D0-RADN**3.0D0)/RADT**3.0D0
            FRAZ3= (RADT**3.0D0-RADG**3.0D0)/RADT**3.0D0
C     calculo vel.de punta de las dendritas
            IF(CLAO.LE.CPRO) THEN ! control
               VGRA= 0.0D0 
            ELSE
               DELT= 2.0D0*DIFCL/VGRA
C     recalculate C_inf from uniform distribution
               IF((RADGO+DELT).GE.RADT) THEN ! capa lim. out of Rt?
                  C1= RADG**3.0D0+4*RADG**2.0D0*DELT+RADG**2.0D0*RADT+
     .                 4.0D0*RADG*RADT*DELT+
     .                 RADG*RADT**2.0D0-3*RADT**3.0D0+
     .                 4.0D0*RADT**2.0D0*DELT
                  C2= RADG**2.0D0+RADG*RADT+RADT**2.0D0
                  C3= RADG**3.0D0+RADT*RADG**2.0D0+RADT**2.0D0*RADG-
     .                 3.0D0*RADT**3.0D0
                  CINF1= (CLA*C1-4.0D0*CPRO*DELT*C2)/C3 
                  CINF = (DELT-(RADT-RADG))*(CLA-CINF1)/DELT+CINF1    
               ELSE
                  C1= (6.0D0*RADGO*RADGO+4.0D0*RADGO*DELT+DELT*DELT)
                  C2= (4.0D0*RADGO**3.0D0+6.0D0*RADGO*RADGO*DELT+
     .                 4.0D0*RADGO*DELT*DELT+DELT**3.0D0-4.0D0*
     .                 RADT**3.0D0)
                  CINF= (CLAO*DELT*C1+4.0D0*CPRO*(RADGO**3.0D0-
     .                 RADT**3.0D0))/C2 
               ENDIF
C     
C**** AUSTENITE GROWING
C     
               IF(IGRAMX.EQ.1) THEN
                  VGRA= DIFCL*AMLQ*(CLAO-CINF)**2.0D0/ ! Rappaz (DEFAULT)
     .                 (PI**2.0D0*COGTH*(AKCA-1.0D0)*CA)
               ELSEIF(IGRAMX.EQ.2) THEN
                  VGRA= AKCA*DIFCL*AMLQ*(CLAO-CINF)**2.0D0/ ! Stefanescu
     .                 (2.0D0*PI**2.0D0*COGTH*(AKCA-1.0D0)*CA)
               ELSE
                  CALL RUNENDT('ERROR IN MICROS9: IGRAMX DOES NOT 
     .                 EXIST')
               ENDIF
C     
               DELT= 2.0D0*DIFCL/VGRA ! espesor capa limite
               GRARG= 2.0D0*(CLAO-CINF)/DELT ! grad de conc. en Rg
            ENDIF               ! (clao <= cpro)
C     
C**** GRAPHITE NUCLEATION
C     
C**** recalscence end nucleation
            IF(INUCAX.EQ.1) THEN ! recalescence nucleation
               IF(JO.GT.0.AND.DTEMPT.GT.0.0D0) INDEXG=1 ! end of nucl.
               IF(INDEXG.EQ.0) THEN ! nucleation
                  JO=JO+1
                  IF(JO.GT.NNUM4T) THEN ! see pointes.f
                     CALL RUNENDT('ERROR IN MICROS9: JO GT NNUM4T')
                  ENDIF
C     
                  IF(INUCMX.EQ.1) THEN ! Su's nucleation law
                     DNU(JO,1)= 0.0D0
                     DNU(JO,2)= ANUCB*GUG*ANUCC*
     .                    (1.0D0-FRAZ1-FRAZ3-FRGZ2)*DTIMET
                     DNU(JO,3)= ANUCB*GUG*ANUCC*
     .                    (1.0D0-FRAZ1-FRAZ2-FRGZ3)*DTIMET
                  ELSEIF(INUCMX.EQ.2) THEN ! Boeri's nucl. law (DEFAULT)
                     DNU(JO,1)= 0.0D0
                     DNU(JO,2)= ANUCB*GUG*DEXP(-ANUCC/GUG)*
     .                    (1.0D0-FRAZ1-FRAZ3-FRGZ2)*DTIMET
                     DNU(JO,3)= ANUCB*GUG*DEXP(-ANUCC/GUG)*
     .                    (1.0D0-FRAZ1-FRAZ2-FRGZ3)*DTIMET
                  ELSEIF(INUCMX.EQ.3) THEN ! Rappaz's nucleation law
                     DNU(JO,1)= 0.0D0 
                     DNU(JO,2)= ANUCB*DEXP(-(GUG-DIFCA)**2/ANUCC)*
     .                    (TGAUAT-TGAUST)*(1.0D0-FRAZ1-FRAZ3-FRGZ2) 
                     DNU(JO,3)= ANUCB*DEXP(-(GUG-DIFCA)**2/ANUCC)*
     .                    (TGAUAT-TGAUST)*(1.0D0-FRAZ1-FRAZ2-FRGZ3)
                  ELSEIF(INUCMX.EQ.4) THEN !  Stefanescu's nucleation law
                     DNU(JO,1)= 0.0D0
                     DNU(JO,2)= 2*ANUCB*ANUCC*DEXP(-ANUCC/GUG**2)/GUG**2
     .                    *(1.0D0-FRAZ1-FRAZ3-FRGZ2)*DTIMET
                     DNU(JO,3)= 2*ANUCB*ANUCC*DEXP(-ANUCC/GUG**2)/GUG**2
     .                    *(1.0D0-FRAZ1-FRAZ2-FRGZ3)*DTIMET
                  ELSE
                     CALL RUNENDT('ERROR IN MICROS9: INUCMX DOES NOT 
     .                    EXIST')
                  ENDIF         ! inucmx.eq...
C     a los nodulos nucleados les asigno el radio inicial
                  RNU(JO)= RNODO 
                  TSOCGNN= 0.0D0 
                  FRGZ2=   0.0D0
                  FRGZ3=   0.0D0
C     calculo fracciones de grafito zona 2 y zona 3
                  DO IJO=1,JO                   
                     FRGZ2= FRGZ2+4.0D0/3.0D0*PI*DNU(IJO,2)*
     .                    RNU(IJO)**3.0D0
                     FRGZ3= FRGZ3+4.0D0/3.0D0*PI*DNU(IJO,3)*
     .                    RNU(IJO)**3.0D0
                  ENDDO
                  SINDEXG= TGAUST ! prueba
               ENDIF            ! indexg.eq... (nuclea el grafito)
C     
C**** T_min end nucleation
C     
            ELSEIF(INUCAX.EQ.2) THEN ! T_min nucleation
               IF(TGAUST.LT.SINDEXG) SINDEXG= TGAUAT
               IF(JO.GT.0.AND.DTEMPT.LT.0.0D0.AND.TGAUST.LE.SINDEXG)THEN
                  JO=JO+1
                  IF(JO.GT.NNUM4T) THEN ! see pointes.f
                     CALL RUNENDT('ERROR IN MICROS9: JO GT NNUM4T')
                  ENDIF
C     
                  IF(INUCMX.EQ.1) THEN ! Su's nucleation law
                     DNU(JO,1)= 0.0D0
                     DNU(JO,2)= ANUCB*GUG*ANUCC*
     .                    (1.0D0-FRAZ1-FRAZ3-FRGZ2)*DTIMET
                     DNU(JO,3)= ANUCB*GUG*ANUCC*
     .                    (1.0D0-FRAZ1-FRAZ2-FRGZ3)*DTIMET
                  ELSEIF(INUCMX.EQ.2) THEN ! Boeri's nucl. law (DEFAULT)
                     DNU(JO,1)= 0.0D0
                     DNU(JO,2)= ANUCB*GUG*DEXP(-ANUCC/GUG)*
     .                    (1.0D0-FRAZ1-FRAZ3-FRGZ2)*DTIMET
                     DNU(JO,3)= ANUCB*GUG*DEXP(-ANUCC/GUG)*
     .                    (1.0D0-FRAZ1-FRAZ2-FRGZ3)*DTIMET
                  ELSEIF(INUCMX.EQ.3) THEN ! Rappaz's nucleation law
                     DNU(JO,1)= 0.0D0 
                     DNU(JO,2)= ANUCB*DEXP(-(GUG-DIFCA)**2/ANUCC)*
     .                    (TGAUAT-TGAUST)*(1.0D0-FRAZ1-FRAZ3-FRGZ2) 
                     DNU(JO,3)= ANUCB*DEXP(-(GUG-DIFCA)**2/ANUCC)*
     .                    (TGAUAT-TGAUST)*(1.0D0-FRAZ1-FRAZ2-FRGZ3)
                  ELSEIF(INUCMX.EQ.4) THEN !  Stefanescu's nucleation law
                     DNU(JO,1)= 0.0D0
                     DNU(JO,2)= 2*ANUCB*ANUCC*DEXP(-ANUCC/GUG**2)/GUG**2
     .                    *(1.0D0-FRAZ1-FRAZ3-FRGZ2)*DTIMET
                     DNU(JO,3)= 2*ANUCB*ANUCC*DEXP(-ANUCC/GUG**2)/GUG**2
     .                    *(1.0D0-FRAZ1-FRAZ2-FRGZ3)*DTIMET
                  ELSE
                     CALL RUNENDT('ERROR IN MICROS9: INUCMX DOES NOT 
     .                    EXIST')
                  ENDIF         ! inucmx.eq...
C     a los nodulos nucleados les asigno el radio inicial
                  RNU(JO)= RNODO 
                  TSOCGNN= 0.0D0 
                  FRGZ2=   0.0D0
                  FRGZ3=   0.0D0
C     calculo fracciones de grafito zona 2 y zona 3
                  DO IJO=1,JO                   
                     FRGZ2= FRGZ2+4.0D0/3.0D0*PI*DNU(IJO,2)*
     .                    RNU(IJO)**3.0D0
                     FRGZ3= FRGZ3+4.0D0/3.0D0*PI*DNU(IJO,3)*
     .                    RNU(IJO)**3.0D0
                  ENDDO
                  SINDEXG= TGAUST ! prueba
               ENDIF            ! sindexg.ne...
            ELSE                ! inucax.eq...
               CALL RUNENDT('ERROR IN MICROS9: INUCAX DOES NOT EXIST')
            ENDIF               ! graphite nucleation ...
C     calculo fraccion de grafito por unidad de zona
            UVVZ1= FRGZ1/FRAZ1
            UVVZ2= FRGZ2/FRAZ2
CCCCC write(94,*) ' UVVZ2_1****   = ', UVVZ2
            UVVZ3= FRGZ3/FRAZ3
C     calculo el deltac
            DELTC= CLA-CLAO     
C     velocidad de la punta de la dendrita
            IF(VGRA.LE.0.0D0) THEN 
               IF(DELTC.LE.0.0D0) THEN
                  CLA= CLAO
                  GOTO 120      ! Rn not growing
               ELSE
                  DELRN= (((1.0D0-AKCA)*CLAO*RADN**3.0D0+
     .                 DELTC*RADG**3.0D0)/ ! ¿RADT O RADGO?
     .                 ((1.0D0-AKCA)*CLAO+DELTC))**(1.0D0/3.0D0)
                  RADN= DELRN   ! Rn growing
                  GOTO 190
               ENDIF            ! deltac.le.0.0D0
            ELSE
            ENDIF               ! fin growth rad. Rn para cdo. vgra.eq.0 ..
C     crece RG
            DELRG= VGRA* DTIMET
            RADG = RADG+ DELRG
C     crecimiento de la austenita eutectica si Rg > Rt
            IF(RADG.GE.RADT) THEN
               RADG= RADT
               DELRN= (((1.0D0-AKCA)*CLAO*RADN**3.0D0+
     .              DELTC*RADT**3.0D0)/ ! ¿RADT O RADGO?
     .              ((1.0D0-AKCA)*CLAO+DELTC))**(1.0D0/3.0D0)
               RADN= DELRN
C     goto 30 -- no hay mas zona 3
               IF(RADNO.LT.RADN.AND.DELTC.GT.0.0D0) GOTO 30
            ELSE
               IF(RADN.EQ.0.0D0) THEN
C     CALCULO DE DELRN
                  DELRN= ((3.0D0*DIFCL*RADGO*RADGO*(1-UAAZ3)*
     .                 GRARG*DTIMET+DELTC*(1.0D0-UVVZ3)*RADG**3.0D0)
     .                 /(((1.0D0-AKCA)*CLA+DELTC)*(1.0D0-UVVZ3))
     .                 )**(1.0D0/3.0D0)
               ELSE
                  DELRN= ((3.0D0*DIFCL*RADGO*RADGO*GRARG*DTIMET+
     .                 (1.0D0-AKCA)*CLAO*RADN**3.0D0+DELTC*RADGO**3.0D0)
     .                 /((1.0D0-AKCA)*CLAO+DELTC))**(1.0D0/3.0D0)
               ENDIF
               RADN= DELRN          
C     goto 170 -- todavia existen las tres zonas
               IF(RADNO.LT.RADN.AND.DELTC.GT.0.0D0) GOTO 170
            ENDIF               ! Rg>Rt
C     
            RADG= RADGO
            RADN= RADNO
C     variacion de Cla y Cpro debido al crecimiento de la austenita
            CLA= (-3.0D0*DIFCL*RADGO*RADGO*GRARG*DTIMET/
     .           (RADGO**3.0D0-RADN**3.0D0))+CLA
C     C            write(94,*) ' cla**   = ', cla
            CPRO= CPRO+3.0D0*DIFCL*RADGO*RADGO*GRARG*DTIMET/
     .           (RADT**3.0D0-RADGO**3.0D0)
            GOTO 120     
 170        CPRO= CPRO+3.0D0*DIFCL*RADGO*RADGO*GRARG*DTIMET/
     .           (RADT**3.0D0-RADGO**3.0D0)
C     como crecio Rn - calculo las nuevas fracciones de zona ...
 190        FRAZ1O= FRAZ1
            FRAZ2O= FRAZ2
            FRAZ3O= FRAZ3
            FRAZ1 = (RADN/RADT)**3.0D0
            FRAZ2 = (RADG**3.0D0-RADN**3.0D0)/RADT**3.0D0
            FRAZ3 = (RADT**3.0D0-RADG**3.0D0)/RADT**3.0D0
            DFRAZ1= (FRAZ1-FRAZ1O)/FRAZ2O
            DFRAZ3= (FRAZ3-FRAZ3O)/FRAZ3O
C     cambio de densidades debido a la variacion de los radios
            DO IJO=1,JO         ! cambio de densidades
               IF(DNU(IJO,1).GT.0.0D0)
     .              RNUZ1(IJO)= ((DNU(IJO,1)*RNUZ1(IJO)**3.0D0+
     .              DNU(IJO,2)*DFRAZ1*RNU(IJO)**3.0D0)/
     .              (DNU(IJO,1)+DNU(IJO,2)*DFRAZ1))**(1.0D0/3.0D0)
               IF(DNU(IJO,1).EQ.0.0D0.AND.DFRAZ1.GT.0.0D0)
     .              RNUZ1(IJO)= RNU(IJO)
               DNU(IJO,1)= DNU(IJO,1)+DNU(IJO,2)*DFRAZ1
               DNU(IJO,2)= DNU(IJO,2)-DNU(IJO,3)*DFRAZ3-DNU(IJO,2)*
     .              DFRAZ1
               DNU(IJO,3)= DNU(IJO,3)+DNU(IJO,3)*DFRAZ3
            ENDDO
C     solidifica por primera vez y no --  aqui se juntan las dos ramas
         ENDIF                  ! end sol. first time or not
C     calculo las nuevas fracciones de gr y volumenes unif. dist.
         FRGZ1= 0.0D0
         FRGZ2= 0.0D0 
         FRGZ3= 0.0D0 
         DO IJO=1,JO   
            FRGZ1= FRGZ1+4.0D0/3.0D0*PI*DNU(IJO,1)
     .           *RNUZ1(IJO)**3.0D0
            FRGZ2= FRGZ2+4.0D0/3.0D0*PI*DNU(IJO,2)
     .           *RNU(IJO)**3.0D0
            FRGZ3= FRGZ3+4.0D0/3.0D0*PI*DNU(IJO,3)
     .           *RNU(IJO)**3.0D0
         ENDDO
C     fracciones de gr por unidad de zona
         UVVZ1= FRGZ1/FRAZ1
         UVVZ2= FRGZ2/FRAZ2
CCCCC write(94,*) ' UVVZ2_2****   = ', UVVZ2
         UVVZ3= FRGZ3/FRAZ3
C     graphite growing in contact with liquid
 120     DRNPRO= 0.0D0
         IF((FRGZ2+FRGZ3).GT.0.0D0) THEN
            DO IJO=1,JO
               IF((DNU(IJO,2)+DNU(IJO,3)).GT.0.0D0) THEN
                  IF(IGRGMX.EQ.1) THEN ! Patri's law
                     DRNUZ2= DIFCL*(CLA-CLG)*DENSA*DTIMET/
     .                    (2.0D0*RNU(IJO)*(100.0D0*DENSG-CLG*DENSA))
                     DRNUZ3= DIFCL*(CINF-CLG)*DENSA*DTIMET/ ! O (CPRO-CLG)
     .                    (RNU(IJO)*(100.0D0*DENSG-CLG*DENSA))
C     la linea que sigue se debe a que al principio cuando se hace cinf=ca 
C     puede que  sea cinf<clg<0 (esto es por cuestiones del diagrama de 
C     equilibrio y de las composiciones iniciales) 
                     IF(DRNUZ3.LT.0.0D0) DRNUZ3= 0.0D0
                  ELSEIF(IGRGMX.EQ.2) THEN ! Zener's law
C     crecimiento de una particula esferica en una matriz sobresaturada 
C     de soluto -- czener -- 
                     DRNUZ2= DIFCL*(CLA-CLG)*DENSA*DTIMET*(1.0D0-FSOLA)/ 
     .                    (2.0D0*RNU(IJO)*DENSG*(100.0D0-CLG))
                     DRNUZ3= DIFCL*(CINF-CLG)*DENSA*DTIMET/ 
     .                    (2.0D0*RNU(IJO)*DENSG*(100.0D0-CLG))
                  ELSE
                     CALL RUNENDT('ERROR IN MICROS9: IGRGMX DOES NOT 
     .                    EXIST')
                  ENDIF
C     nuevo calculo repartiendo el flujo -- radio promedio --
                  DRNPRO=(DNU(IJO,2)*DRNUZ2+DNU(IJO,3)*DRNUZ3)/
     .                 (DNU(IJO,2)+DNU(IJO,3))
                  RNU(IJO)= RNU(IJO)+DRNPRO
               ELSE
               ENDIF            ! (dnu(ijo,2)+dnu(ijo,3)).gt. ...
            ENDDO
         ELSE
         ENDIF                  ! fin crec. del graf. in contact wtih liq.
C     fracciones de grafito por zona
         FRGZ1= 0.0D0
         FRGZ2= 0.0D0
         FRGZ3= 0.0D0
C     
         DO IJO=1,JO                    
            FRGZ1= FRGZ1+4.0D0/3.0D0*PI*DNU(IJO,1)
     .           *RNUZ1(IJO)**3.0D0 ! calculate using average radio
            FRGZ2= FRGZ2+4.0D0/3.0D0*3.14*DNU(IJO,2)
     .           *RNU(IJO)**3.0D0 
            FRGZ3= FRGZ3+4.0D0/3.0D0*PI*DNU(IJO,3)*
     .           RNU(IJO)**3.0D0 ! calculate using average radio
         ENDDO
C     fracciones de gr por unidad de zona
         UVVZ1O= UVVZ1
         UVVZ2O= UVVZ2
         UVVZ3O= UVVZ3
         UVVZ1 = FRGZ1/FRAZ1
         UVVZ2 = FRGZ2/FRAZ2
CCCCC write(94,*) ' UVVZ2_3****   = ', UVVZ2
         UVVZ3 = FRGZ3/FRAZ3
C     new concentration of C in the liquid (z1 and z2) then 
C     graphite growing 
         CPRO= (CPRO*(1.0D0-UVVZ3O)*DENSA+100.0D0*(UVVZ3O-UVVZ3)
     .        *DENSG)/((1.0D0-UVVZ3)* DENSA) 
         CLA = (CLA*(1.0D0-UVVZ2O)*DENSA+100.0D0*(UVVZ2O-UVVZ2)
     .        *DENSG)/((1.0D0-UVVZ2)* DENSA)
CCC   write(91,*) ' cla***   = ', cla
CCC   write(91,*) 'frgz1_33      = ', frgz1
C     
CCC   write(94,*) ' cla**   = ', cla
C     after growth graphite and austenite goto 70
         GOTO 70                ! calculo --> f_sol, f_aust, f_gr, Si
C     from first step of solidified when Rg>=Rt     
 20      RADN  = 0.0D0
         CLAO  = CPRO
         INDEXR= 1              ! Rg alcanzo el valor de Rt
         RADG  = RADT
         FRAZ3 = 0.0D0
         DELTC = CLA-CLAO
         DELRN = ((DELTC*RADT**3.0D0)/ ! falta otra formula
     .        ((1.0D0-AKCA)*CLA))**(1.0D0/3.0D0)
         RADN= DELRN
C     control on Rn
         IF(RADN.GT.RADT) THEN
            FSOLID= 1.0D0
            FSOLG = DNU(1,3)*RNU(1)**3.0D0
            FSOLA = FSOLID-FSOLG
CCCC  write(92,*) 'fsola_1      = ', fsola
            GOTO 110
         ELSE
            FRAZ1   = (RADN/RADT)**3.0D0
            FRAZ2   = (RADG**3.0D0-RADN**3.0D0)/RADT**3.0D0 
            DNU(1,1)= DNU(1,3)*FRAZ1
            DNU(1,2)= DNU(1,3)*FRAZ2    
         ENDIF
C     
         GOTO 150
 30      INDEXR= 1              ! cdo. Rg>Rt por primera vez (first time)
C     calc DELRN
         DELRN= (((1.0D0-AKCA)*CLAO*RADN**3.0D0+
     .        DELTC*RADT**3.0D0)/
     .        ((1.0D0-AKCA)*CLAO+DELTC))**(1.0D0/3.0D0)
         RADN  = DELRN
         FRAZ1O= FRAZ1
         FRAZ2O= FRAZ2
         FRAZ3O= FRAZ3
         FRAZ1 = (RADN/RADT)**3.0D0
         FRAZ2 = (RADG**3.0D0-RADN**3.0D0)/RADT**3.0D0
         FRAZ3 = 0.0D0          ! zone 3 doesn't exist 
         DFRAZ1= (FRAZ1-FRAZ1O)/FRAZ2O
         DFRAZ3= (FRAZ3-FRAZ3O)/FRAZ3O ! esto esta mal !
C     density changes
         DO IJO=1,JO 
            IF(DNU(IJO,1).GT.0.0D0)
     .           RNUZ1(IJO)= ((DNU(IJO,1)*RNUZ1(IJO)**3.0D0+
     .           DNU(IJO,2)*DFRAZ1*RNU(IJO)**3.0D0)/
     .           (DNU(IJO,1)+DNU(IJO,2)*DFRAZ1))**(1.0D0/3.0D0)
C     
            IF(DNU(IJO,1).EQ.0.0D0.AND.DFRAZ1.GT.0.0D0)
     .           RNUZ1(IJO)= RNU(IJO)
            DNU(IJO,1)= DNU(IJO,1)+DNU(IJO,2)*DFRAZ1
            DNU(IJO,2)= DNU(IJO,2)+DNU(IJO,3)-DNU(IJO,2)*DFRAZ1
            DNU(IJO,3)= 0.0D0
         ENDDO
         FRGZ3= 0.0D0           ! zone 3 doesn't exist
         GOTO 150              
C     calc new fractions of graphite
 50      FRGZ1= 0.0D0
         FRGZ2= 0.0D0
         FRGZ3= 0.0D0
         DO IJO=1,JO                   
            FRGZ1= FRGZ1+DNU(IJO,1)*RNUZ1(IJO)**3.0D0
            FRGZ2= FRGZ2+DNU(IJO,2)*RNU(IJO)**3.0D0
         ENDDO
         FRGZ1= FRGZ1*4.0D0/3.0D0*PI
         FRGZ2= FRGZ2*4.0D0/3.0D0*PI
C     final fractions
         FSOLG= FRGZ1+ FRGZ2
C     
         FRAZ1= FSOLID-FRGZ2
         RADN = ALPHAM(IN+10)
         FRAZ2= (RADT**3.0D0-RADN**3.0D0)/RADT**3.0D0 ! == 1.0D0-FRAZ1
         FRAZ3= 0.0D0           ! zone 3 doesn't exist
C     
C**** GRAPHITE NUCLEATION
C     
         IF(INUCAX.EQ.1) THEN   ! recalescence nucleation
            IF(JO.GT.0.AND.DTEMPT.GT.0.0D0) INDEXG=1 ! end of nucl.
            IF(INDEXG.EQ.0) THEN ! nucleation
               JO=JO+1
               IF(JO.GT.NNUM4T) THEN ! see pointes.f
                  CALL RUNENDT('ERROR IN MICROS9: JO GT NNUM4T')
               ENDIF
C     
               IF(INUCMX.EQ.1) THEN ! Su's nucleation law
                  DNU(JO,1)= 0.0D0
                  DNU(JO,2)= ANUCB*GUG*ANUCC*(1.0D0-FSOLID)*DTIMET
                  DNU(JO,3)= 0.0D0
               ELSEIF(INUCMX.EQ.2) THEN ! Boeri's nucleation law (DEFAULT)
                  DNU(JO,1)= 0.0D0
                  DNU(JO,2)= ANUCB*GUG*DEXP(-ANUCC/GUG)*
     .                 (1.0D0-FSOLID)*DTIMET
                  DNU(JO,3)= 0.0D0
               ELSEIF(INUCMX.EQ.3) THEN ! Rappaz's nucleation law
                  DNU(JO,1)= 0.0D0 
                  DNU(JO,2)= ANUCB*DEXP(-(GUG-DIFCA)**2/ANUCC)*
     .                 (TGAUAT-TGAUST)*(1.0D0-FSOLID)**2.0D0
                  DNU(JO,3)= 0.0D0
               ELSEIF(INUCMX.EQ.4) THEN !  Stefanescu's nucleation law
                  DNU(JO,1)= 0.0D0
                  DNU(JO,2)= 2*ANUCB*ANUCC*DEXP(-ANUCC/GUG**2)/
     .                 GUG**2.0D0*(1.0D0-FSOLID)*DTIMET
                  DNU(JO,3)= 0.0D0
               ELSE
                  CALL RUNENDT('ERROR IN MICROS9: INUCMX DOES NOT 
     .                 EXSIT')
               ENDIF            ! inucmx.eq...
C     a los nodulos nucleados (por cualquiera de los dos criterios)
C     les asigno el radio inicial
               RNU(JO)= RNODO
C     calculo fracciones de grafito zona 2 (zone3=0)
               FRGZ2= 0.0D0
               FRGZ3= 0.0D0
               DO IJO=1,JO                   
                  FRGZ2= FRGZ2+4.0D0/3.0D0*PI*DNU(IJO,2)
     .                 *RNU(IJO)**3.0D0
               ENDDO
               SINDEXG= TGAUST  ! prueba
            ENDIF               ! indexg.eq...
         ELSEIF(INUCAX.EQ.2) THEN ! T_min nucleation (default)
            IF(TGAUST.LT.SINDEXG) SINDEXG=TGAUAT
            IF(JO.GT.0.AND.DTEMPT.LT.0.0D0.AND.TGAUST.LE.SINDEXG) THEN
               JO=JO+1
               IF(JO.GT.NNUM4T) THEN ! see pointes.f
                  CALL RUNENDT('ERROR IN MICROS9: JO GT NNUM4T')
               ENDIF
C     
               IF(INUCMX.EQ.1) THEN ! Su's nucleation law
                  DNU(JO,1)= 0.0D0
                  DNU(JO,2)= ANUCB*GUG*ANUCC*(1.0D0-FSOLID)*DTIMET
                  DNU(JO,3)= 0.0D0
               ELSEIF(INUCMX.EQ.2) THEN ! Boeri's nucleation law (DEFAULT)
                  DNU(JO,1)= 0.0D0
                  DNU(JO,2)= ANUCB*GUG*DEXP(-ANUCC/GUG)*
     .                 (1.0D0-FSOLID)*DTIMET
                  DNU(JO,3)= 0.0D0
               ELSEIF(INUCMX.EQ.3) THEN ! Rappaz's nucleation law
                  DNU(JO,1)= 0.0D0 
                  DNU(JO,2)= ANUCB*DEXP(-(GUG-DIFCA)**2.0D0/ANUCC)*
     .                 (TGAUAT-TGAUST)*(1.0D0-FSOLID)**2.0D0
                  DNU(JO,3)= 0.0D0
               ELSEIF(INUCMX.EQ.4) THEN !  Stefanescu's nucleation law
                  DNU(JO,1)= 0.0D0
                  DNU(JO,2)= 2*ANUCB*ANUCC*DEXP(-ANUCC/GUG**2.0D0)/
     .                 GUG**2.0D0*(1.0D0-FSOLID)*DTIMET
                  DNU(JO,3)= 0.0D0
               ELSE
                  CALL RUNENDT('ERROR IN MICROS9: INUCMX DOES NOT 
     .                 EXIST')
               ENDIF            ! inucmx.eq...
C     a los nodulos nucleados (por cualquiera de los dos criterios)
C     les asigno el radio inicial
               RNU(JO)= RNODO
C     calculo fracciones de grafito zona 2 (zone3=0)
               FRGZ2= 0.0D0
               FRGZ3= 0.0D0
               DO IJO=1,JO                   
                  FRGZ2= FRGZ2+4.0D0/3.0D0*PI*DNU(IJO,2)
     .                 *RNU(IJO)**3.0D0
               ENDDO
               SINDEXG= TGAUST  ! prueba
            ENDIF               ! sindexg.ne...
         ELSE                   ! inucax.eq...
            CALL RUNENDT('ERROR IN MICROS9: INUCAX DOES NOT EXIST')
         ENDIF                  ! graphite nucleation ...
C     calc. el crecimiento promedio de los nodulos nucleados en zona 2 y 3
         DELTC= CLA-CLAO
         IF(DELTC.LE.0.0D0) then
            CLA= CLAO
            GOTO 180
         ELSE
         ENDIF
C     calculo del DELRN promedio -- ver micros9.f --  deberia elegir.
         DELRN= (((1.0D0-AKCA)*CLAO*RADN**3.0D0+DELTC*RADT**3.0D0)/
     .        ((1.0D0-AKCA)*CLAO+DELTC))**(1.0D0/3.0D0)
         RADN= DELRN
C     si el radio Rn alcanzo el radio Rt nos vamos a 60 donde assign RN=RT
         IF(RADN.GE.RADT) GOTO 60 ! calc. la fraccion solida
         FRAZ1O= FRAZ1
         FRAZ2O= FRAZ2
         FRAZ1 = (RADN/RADT)**3.0D0
         FRAZ2 = (RADT**3.0D0-RADN**3.0D0)/RADT**3.0D0 ! == 1.0D0-FRAZ1
         DFRAZ1= (FRAZ1-FRAZ1O)/FRAZ2O 
C     crecimiento de Rn -- varian las densidades de zona
         DO IJO=1,JO            ! growth 
            IF(DNU(IJO,1).GT.0.0D0)
     .           RNUZ1(IJO)= ((DNU(IJO,1)*RNUZ1(IJO)**3.0D0+
     .           DNU(IJO,2)*DFRAZ1*RNU(IJO)**3.0D0)/
     .           (DNU(IJO,1)+DNU(IJO,2)*DFRAZ1))**(1.0D0/3.0D0)
C     
            IF(DNU(IJO,1).EQ.0.0D0.AND.DFRAZ1.GT.0.0D0)
     .           RNUZ1(IJO)= RNU(IJO)
            DNU(IJO,1)= DNU(IJO,1)+DNU(IJO,2)*DFRAZ1
            DNU(IJO,2)= DNU(IJO,2)-DNU(IJO,2)*DFRAZ1
         ENDDO
C     fracciones gr zona 1 y zona 2
 150     FRGZ1= 0.0D0
         FRGZ2= 0.0D0 
         DO IJO=1,JO  
            FRGZ1= FRGZ1+DNU(IJO,1)*RNUZ1(IJO)**3.0D0
            FRGZ2= FRGZ2+DNU(IJO,2)*RNU(IJO)**3.0D0
         ENDDO
         FRGZ1= FRGZ1*4.0D0/3.0D0*PI
         FRGZ2= FRGZ2*4.0D0/3.0D0*PI
C     fracciones de gr por unidad de zona 1 y 2
         UVVZ1= FRGZ1/FRAZ1
         UVVZ2= FRGZ2/FRAZ2
CCCC  write(94,*) ' UVVZ2_4****   = ', UVVZ2
CCCC  write(94,*) ' fragz2_4****  = ', frgz2
CCCC  write(94,*) ' fraz2_4****   = ', fraz2
C     
C**** graphite growing in contact with liquid -- only in zone 2
C     
 180     IF(IGRGMX.EQ.1) THEN   ! Patri's law (default)
            DO IJO=1,JO
               DRNUZ2=DIFCL*(CLA-CLG)*DENSA*DTIMET*(1.0D0-FSOLA)/
     .              (2.0D0*RNU(IJO)*(100.0D0*DENSG-CLG*DENSA))
C     nuevo calculo solo zona 2
               RNU(IJO)= RNU(IJO)+ DRNUZ2
            ENDDO
         ELSEIF(IGRGMX.EQ.2) THEN ! Zener's law
            DO IJO=1,JO
               DRNUZ2= DIFCL*(CLA-CLG)*DENSA*DTIMET*(1.0D0-FSOLID)/ 
     .              (2.0D0*RNU(IJO)*DENSG*(100.0D0-CLG))
C     nuevo calculo solo zona 2
               RNU(IJO)= RNU(IJO)+ DRNUZ2
            ENDDO
         ELSE
            CALL RUNENDT('ERROR IN MICROS9: IGRGMX DOES NOT EXIST')
         ENDIF
C     calculamos las fracciones de grafito de zona 1 y 2 --
         FRGZ1= 0.0D0
         FRGZ2= 0.0D0 
         DO IJO=1,JO                    
            FRGZ1= FRGZ1+DNU(IJO,1)*RNUZ1(IJO)**3.0D0
            FRGZ2= FRGZ2+DNU(IJO,2)*RNU(IJO)**3.0D0
         ENDDO
         FRGZ1= FRGZ1*4.0D0/3.0D0*PI
         FRGZ2= FRGZ2*4.0D0/3.0D0*PI
C     fracciones gr only in zone 2
         UVVZ2O= UVVZ2
         UVVZ2 = FRGZ2/FRAZ2
C     
CCC   write(95,*) ' cla****_a    = ', cla
CCC   write(95,*) ' UVVZ2O****_a    = ', UVVZ2O
CCC   write(95,*) ' UVVZ2****_a    = ', UVVZ2         
C     calculamos la nueva concentracion CLA
CCCC  write(93,*) ' cla_a****    = ', cla
CCCC  write(93,*) ' UVVZ2O_a**** = ', UVVZ2O
CCCC  write(93,*) ' radt_a****   = ', radt
CCCC  write(93,*) ' UVVZ2_a****  = ', UVVZ2         
C     calculamos la nueva concentracion CLA
C     
         CLA= (CLA*(1.0D0-UVVZ2O)*DENSA+100.0D0*(UVVZ2O-UVVZ2)
     .        *DENSG)/((1.0D0-UVVZ2)* DENSA)
CCCC  write(93,*) ' cla_d****    = ', cla  
CCCC  write(94,*) ' cla****   = ', cla
CCCC  write(94,*) ' uvvz2****   = ', uvvz2
CCCC  write(94,*) ' frgz2****   = ', frgz2     
CCC   write(95,*) ' DENSA     = ', densa
CCC   write(95,*) ' DENSG     = ', densg
CCC   write(95,*) ' fsolid = ', fsolid
C     calculate f_sol - f_aust - f_gr - Si
         GOTO 70                ! calculo --> f_sol, f_aust, f_gr, Si
C     
C**** goto 60 - we have only zone 1
C     
 60      RADN  = RADT           ! assign radn to radius of unit cell - radt
         FRAZ1O= FRAZ1
         FRAZ1 = 1.0D0
         DFRAZ1=(FRAZ1-FRAZ1O)/FRAZ2
C     growth only in zone 1
         DO IJO=1,JO            ! growth 
            IF(DNU(IJO,1).GT.0.0D0)
     .           RNUZ1(IJO)= ((DNU(IJO,1)*RNUZ1(IJO)**3.0D0+
     .           DNU(IJO,2)*DFRAZ1*RNU(IJO)**3.0D0)/
     .           (DNU(IJO,1)+DNU(IJO,2)*DFRAZ1))**(1.0D0/3.0D0)
            IF(DNU(IJO,1).EQ.0.0D0.AND.DFRAZ1.GT.0.0D0)
     .           RNUZ1(IJO)= RNU(IJO)
            DNU(IJO,1)= DNU(IJO,1)+DNU(IJO,2)
            DNU(IJO,2)= 0.0D0
         ENDDO
C     
         FRAZ2= 0.0D0
         FRGZ1= 0.0D0
         FRGZ2= 0.0D0 
         FRGZ3= 0.0D0
         DO IJO=1,JO                    
            FRGZ1= FRGZ1+DNU(IJO,1)*RNUZ1(IJO)**3.0D0
         ENDDO
         FRGZ1= FRGZ1*4.0D0/3.0D0*PI
C     
         FSOLG = FRGZ1
         FSOLID= 1.0D0
         FSOLA = FSOLID-FSOLG
CCCC  write(92,*) 'fsola_2      = ', fsola
         FLIQD = 0.0D0
C     
         GOTO 160
C     
C**** end only zone 1
C     
C     only useful for thermal jacobian
 70      TSOCGN= 0.0D0
         TSOCGG= 0.0D0
         TSOCGA= 0.0D0
C     calculate graphite_fraction & austenite_fraction % solid_fraction
 90      FSOLG = FRGZ1+ FRGZ2+ FRGZ3
         FSOLA = FRAZ1-FRGZ1
CCCC  write(92,*) 'fsola_3 = ', fsola
CCCC  write(92,*) 'frgz1_3 = ', frgz1
CCCC  write(92,*) 'frgz2_3 = ', frgz2
CCCC  write(92,*) 'frgz3_3 = ', frgz3
CCCC  write(92,*) 'fraz1_3 = ', fraz1
         FSOLID= FRAZ1+ FRGZ2+ FRGZ3
C     control
         IF(FSOLID.GE.0.9990D0) THEN
            FSOLID= 1.0D0
            FSOLA = FSOLID-FSOLG
CCCC  write(92,*) 'fsola_4      = ', fsola
         ENDIF
      ENDIF                     ! gug.gt.0.0.and....
C     
C**** LIQUID FRACTION
C     
 110  FLIQD= 1.0D0-FSOLID       ! liquid fraction
C     controls
      IF(FSOLID.EQ.1.0D0) FLIQD= 0.0D0
      IF(FLIQD.LT.0.0D0)  FLIQD= 0.0D0 ! control
      IF(FLIQD.GT.1.0D0)  FLIQD= 1.0D0 ! control
C     
C**** DEFINES THE MICROSEGREGATION MODELS 
C**   en esta subrutina calculo las nuevas concentraciones 
C**   en base a la particion de los elementos de aleacion
C     
C     modelos de segregacion implementados
C     iMicrSegr= 1 -- LEVER RULE --
C     iMicrSegr= 2 -- SCHEIL RULE (1942) --
C     iMicrSegr= 3 -- BRODY-FLEMINGS RULE (1996) --
C     iMicrSegr= 4 -- OHNAKA RULE (1986) --
C     iMicrSegr= 5 -- WITHOUT SEGREGATION --
C     
C     hago segregar hasta un 90 % de la fraccion solida
C     
C***  CALCULATE SECONDARY DENDRITE ARM SPACING (SDAS)
C     C         if(fsolid.ge.1.0D-4) then
C     call sdasGrowthRate(isdasgr,difcl,cogth,cla,ca,amlq,akca,tf,
 160  call sdasGrowthRate(isdasgr,difcl,cogth,cla,ca,amlq,akca,
     .     totalTime,sdas, sdasini)
C     endif
      if(fsolid.le.9.99D-1) then
CCCC  write(92,*) 'sicon_antes = ', sicon         
CCCC  write(92,*) 'pdas_antes  = ', pdas
C     SI CONTENT (SI's MICROSEGREGATION)  silicium
         call micrSegrL(SICONX, ISILSMX, AKSIL, SI, FSOLID, tf, 
     .        DIFCA, PDAS, SDAS)
C     
CCCCC SICON= SI*(1.0D0-FSOLID)**(AKSIL-1.0D0)
C     
CCCC  write(92,*) 'si_despues   = ', si         
CCCC  write(92,*) 'sicon_despues= ', sicon         
CCCC  write(92,*) 'fsolid       = ', fsolid
CCCC  write(92,*) 'fsola        = ', fsola
C**** elementos que segregand e forma directa
C     
C     CU CONTENT (Cu's MICROSEGREGATION) copper
         call micrSegrL(CUCONX, ICULSMX, AKCUL, CU, FSOLID, tf, 
     .        DIFCA, PDAS, SDAS)
C     NI CONTENT (NI's MICROSEGREGATION)
         call micrSegrL(QICONX, IQILSMX, AKQIL, QI, FSOLID, tf, 
     .        DIFCA, PDAS, SDAS)
C     
C***  para los elementos que segregan de forma inversa, adopto este 
C     criterio porque sino cuando fsolid=1.0D0 y al ser el coeficiente 
C     de particion menor que uno, el exponente da negativo y 
C     (1.0D0-fsolid)=0.0D0 con lo cual me da negativo
C     
C     MN CONTENT (Mn's MICROSEGREGATION) manganese
         call micrSegrL(QNCONX, IQNLSMX, AKQNL, QN, FSOLID, tf,
     .        DIFCA, PDAS, SDAS)
C     NB CONTENT (NB's MICROSEGREGATION) niobium
         call micrSegrL(QBCONX, IQBLSMX, AKQBL, QB, FSOLID, tf, 
     .        DIFCA, PDAS, SDAS)
C     SN CONTENT (SN's MICROSEGREGATION) estano
         call micrSegrL(SNCONX, ISNLSMX, AKSNL, SN, FSOLID, tf, 
     .        DIFCA, PDAS, SDAS)
C     CH CONTENT (CH's MICROSEGREGATION) chromiun
         call micrSegrL(CRCONX, ICRLSMX, AKCRL, CR, FSOLID, tf, 
     .        DIFCA, PDAS, SDAS)
C     MO CONTENT (MO's MICROSEGREGATION)
         call micrSegrL(QOCONX, IQOLSMX, AKQOL, QO, FSOLID, tf, 
     .        DIFCA, PDAS, SDAS)
C**   assign values to local variables
         QOCON= QOCONX
         CRCON= CRCONX
         SNCON= SNCONX
         QBCON= QBCONX
         QNCON= QNCONX
         QICON= QICONX
         CUCON= CUCONX
         SICON= SICONX
      endif
C     
C**** DEFINES MICROSTRUCTURAL-DEPENDENT "MACROSCOPICAL" PROPERTIES
C     for the moment conductivity doesn't depend of micro.
      IF(IKMICX.EQ.1) THEN     
         IF(IMICOUP.EQ.1) THEN  ! Boeri's
            IF(FLIQD.LE.0.0D0)                     BASKK(1)=BASKS
            IF(FLIQD.GT.0.0D0.AND.FLIQD.LT.1.0D0)  BASKK(1)=BASKM
            IF(FLIQD.GE.1.0D0)                     BASKK(1)=BASKL
         ELSEIF(IMICOUP.EQ.2) THEN ! old fl (DEFAULT)     
            IF(FLIQDO.LE.0.0D0)                    BASKK(1)= BASKS
            IF(FLIQDO.GT.0.0D0.AND.FLIQD.LT.1.0D0) BASKK(1)= BASKM
            IF(FLIQDO.GE.1.0D0)                    BASKK(1)= BASKL
         ELSEIF(IMICOUP.EQ.3) THEN ! smooth
            IF(FLIQD.LE.0.0D0)                     BASKK(1)= BASKS
            IF(FLIQD.GT.0.0D0.AND.FLIQD.LT.1.0D0)  BASKK(1)= BASKS+
     .           (BASKL-BASKS)*FLIQD
            IF(FLIQD.GE.1.0D0)                     BASKK(1)= BASKL
         ELSE
            CALL RUNENDT('ERROR IN MICROS9: IMICOUP DOES NOT EXIST')
         ENDIF
      ENDIF
C     
C**** ESTABLISHES LATENT HEAT * PHASE-CHANGE FUNCTION (at time t & t+dt)
C     
      TSOE1(IPLAT)= FLIQD       ! f_pc at time t+dt
      TSOE2(IPLAT)= FLIQDO      ! f_pc at time t
      TSOE2(IPLAT)= TSOE2(IPLAT)*HENER ! L*f_pc at time t
      TSOE1(IPLAT)= TSOE1(IPLAT)*HENER ! L*f_pc at time t+dt
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
      GOTO (1,2,3,4,5) IFPCDT
C     
 1    ICACOT= 0
      VELTOT= 1.D-10
      VELaA1T= DABS(DTEMPT)
      IF(VELA1T.GT.VELTOT) THEN
         IF(IITERT.GT.0)
     .        TSOC1(IPLAT)=(TSOE1(IPLAT)-TSOE2(IPLAT))/(DTEMPT*DTIMET)
      ENDIF
      GOTO 10
C     
 2    ICACOT=1
      VELTOT=1.D-10
      VELA1T=DABS(DTEMPT)
      IF(VELA1T.GT.VELTOT) THEN
         IF(IITERT.GT.0)
     .        TSOC1(IPLAT)=(TSOE1(IPLAT)-TSOE2(IPLAT))/(DTEMPT*DTIMET)
      ENDIF
      GOTO 10
C     
 3    ICACOT=0
      TSOC1(IPLAT)=TSOCGN+TSOCGG+TSOCGA
      TSOC1(IPLAT)=TSOC1(IPLAT)*HENER
      GOTO 10
C     
 4    ICACOT=1
      TSOC1(IPLAT)=TSOCGN+TSOCGG+TSOCGA
      TSOC1(IPLAT)=TSOC1(IPLAT)*HENER
      GOTO 10
C     
 5    ICACOT=1
      TSOC1(IPLAT)=0.0D0
      GOTO 10
C     
 10   CONTINUE
C     
C**** TRANSFER MICROSCOPICAL VARIABLES TO "ALPHAM"
C     
      ALPHAM(IN+1)= FLIQD
C     
      ALPHAM(IN+2)= FSOLA
      ALPHAM(IN+3)= FSOLG
      ALPHAM(IN+4)= SICON       ! with Si segregation or DGRA
      ALPHAM(IN+5)= RADT
      ALPHAM(IN+6)= RADG
      ALPHAM(IN+7)= CLA
C     
      IF(VGRA.EQ.0.0D0) THEN
         ALPHAM(IN+8)= VGRAO    ! pq sino en el prox. paso delta --> inf.
      ELSE
         ALPHAM(IN+8)= VGRA
      ENDIF
C**   
      ALPHAM(IN+9) = CPRO
      ALPHAM(IN+10)= RADN
C***  new variables
      ALPHAM(IN+11)= RADC       ! RADC -- cellular radius (see Rivera's Thesis)
      ALPHAM(IN+12)= PDAS       ! ISDASGR - - i model of sdas 
      ALPHAM(IN+13)= totalTime  ! total time from pouring
      ALPHAM(IN+14)= SDAS       ! SDAS - secondary dendrite arm spacing
C     SDAS - secondary dendrite arm spacing -- initial for every step
      ALPHAM(IN+15)= SDASINI
      ALPHAM(IN+16)= PHCON      ! phosporus content (current) - PHCON
      ALPHAM(IN+17)= CUCON      ! cooper content (current) - CUCON
      ALPHAM(IN+18)= QNCON      ! manganese content (current) - QNCON
      ALPHAM(IN+19)= QGCON      ! magnesium (magnesio) content (current)- HGCON
      ALPHAM(IN+20)= QBCON      ! niobium content (current) - NBCON
      ALPHAM(IN+21)= SNCON      ! tin (estano) content (current)- SNCON
      ALPHAM(IN+22)= CRCON      ! chromium  (cromo) content (current)- CRCON
      ALPHAM(IN+23)= QOCON      ! molybdenum (molibdeno) content (current)- QOCON
      ALPHAM(IN+24)= QICON      ! nickel (nikel) content (current)- QICON
      ALPHAM(IN+25)= NUY1       ! not used yet 1
      ALPHAM(IN+26)= NUY2       ! not used yet 2
C     
      IF(JO.GT.0) THEN
         DO INU=1,JO
            ALPHAM(IN+27+INU-1)             = DNU(INU,1)
            ALPHAM(IN+28+ NNUM4T-1+INU-1)   = DNU(INU,2)
            ALPHAM(IN+29+2*(NNUM4T-1)+INU-1)= DNU(INU,3)
            ALPHAM(IN+30+3*(NNUM4T-1)+INU-1)= RNU(INU)
            ALPHAM(IN+31+4*(NNUM4T-1)+INU-1)= RNUZ1(INU)
         ENDDO
      ENDIF
C     
      ALPHAM(IN+32+5*(NNUM4T-1))= FLOAT(JO) ! graphite nucl. index - jo
      ALPHAM(IN+33+5*(NNUM4T-1))= FLOAT(INDEXG) ! recalescence nucleation 
      ALPHAM(IN+34+5*(NNUM4T-1))= SINDEXG !  temp. minima nucleacion 
      ALPHAM(IN+35+5*(NNUM4T-1))= FLOAT(INDEXR) ! grap. Rg growth index INDEXR=1
      ALPHAM(IN+36+5*(NNUM4T-1))= TF ! local final time of solidification
C     
      DTIMET= DTIMETO
C     
C**** INCREMENTS "ALPHAM" INDEX
C     
      INUPC= INUPC+26+5*NNUM4T+5 ! less than NBASES; see pointes.f
C     
      RETURN
      END
