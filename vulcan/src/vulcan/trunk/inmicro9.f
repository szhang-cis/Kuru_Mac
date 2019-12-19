      SUBROUTINE INMICRO9(ALPHAM,INUPC,IPLAT)
C***********************************************************************
C     
C**** THIS ROUTINE INITIALITES MICROSTRUCTURAL PARAMETERS FOR MODEL 9
C     
C***********************************************************************
C     
C     Index of variables
C     
C     ALPHAM= array of microstructural (microscopical) variables
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
C     
C**** TRANSFERS "VPLAT" TO MICROSCOPICAL VARIABLES
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
      CR    = VPLAT(IPLAT,46)
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
C     
      INUCMX= INT(VPLAT(IPLAT,61))
      ANUCA = VPLAT(IPLAT,62)
      ANUCB = VPLAT(IPLAT,63)
      ANUCC = VPLAT(IPLAT,64)
C     
      INUCAX= INT(VPLAT(IPLAT,65))
C     
      IGRGMX= INT(VPLAT(IPLAT,66))
      DIFCL = VPLAT(IPLAT,67)
      RNODA = VPLAT(IPLAT,68)
      RNODO = VPLAT(IPLAT,69)
      AUSGR = VPLAT(IPLAT,70)
      DENSA = VPLAT(IPLAT,71)
      DENSG = VPLAT(IPLAT,72)
C     
      IGRAMX= INT(VPLAT(IPLAT,73))
      DIFCA = VPLAT(IPLAT,74)
      AMLQ  = VPLAT(IPLAT,75)
      AKCA  = VPLAT(IPLAT,76)
      COGTH = VPLAT(IPLAT,77)
C     
      IDELRN = INT(VPLAT(IPLAT,78))
      ISDASGR= INT(VPLAT(IPLAT,79))
      IMICOUP= INT(VPLAT(IPLAT,80))
C     
      IKMICX = INT(VPLAT(IPLAT,81))
      IKAUX= 0
      IF(IKMICX.EQ.1) THEN
         IF(IMICOUP.EQ.1.OR.IMICOUP.EQ.2.OR.IMICOUP.EQ.3) THEN
            IKAUX= 3
            VPLAT(IPLAT,82)= BASKS
            VPLAT(IPLAT,83)= BASKM
            VPLAT(IPLAT,84)= BASKL
         ENDIF
      ENDIF
C     
      IFPCDT= INT(VPLAT(IPLAT,82+IKAUX))
      IAFLOJ= INT(VPLAT(IPLAT,83+IKAUX))
C     
      IV= 83+IKAUX              ! number of VPLAT defined in input data
C     
C**** INITIALITES MICROSTRUCTURAL VARIABLES ("ALPHAM")
C     
      IN           = INUPC
      ALPHAM(IN+ 1)= 1.0D0      ! assumed liq. state as initial condition
      ALPHAM(IN+ 2)= 0.0D0      ! assumed zero initial austenite fraction
      ALPHAM(IN+ 3)= 0.0D0      ! assumed zero initial graphite fraction
      ALPHAM(IN+ 4)= SI         ! SICON -- Si content
      ALPHAM(IN+ 5)= 0.0D0      ! total grain radius - RADT
      ALPHAM(IN+ 6)= 0.0D0      ! grain radius Rg - RADG
      ALPHAM(IN+ 7)= 0.0D0      ! last step C1/a concentration - Cl/a
      ALPHAM(IN+ 8)= 0.0D0      ! dendrite tip growth rate - VGRA
      ALPHAM(IN+ 9)= 0.0D0      ! inter-grain liq. C concentration - Cpro
      ALPHAM(IN+10)= 0.0D0      ! Rn - RADN
C     
      ALPHAM(IN+11)= 0.0D0      ! RADC -- cellular radius (see Rivera's Thesis)
      ALPHAM(IN+12)= 0.0D0      ! PDAS - primary dendrite arm spacing
      ALPHAM(IN+13)= 0.0D0      ! total time from pouring
      ALPHAM(IN+14)= 0.0D0      ! SDAS - secondary dendrite arm spacing
      ALPHAM(IN+15)= 0.0D0      ! SDASINI - secondary dendrite arm spacing
      ALPHAM(IN+16)= PH         ! phosporus content (current) - PHCON
      ALPHAM(IN+17)= CU         ! cooper content (current) - CUCON
      ALPHAM(IN+18)= QN         ! manganese content (current) - QNCON
      ALPHAM(IN+19)= QG         ! magnesium (magnesio) content (current)- HGCON
      ALPHAM(IN+20)= QB         ! niobium content (current) - NBCON
      ALPHAM(IN+21)= SN         ! tin (estano) content (current) - SNCON
      ALPHAM(IN+22)= CR         ! chromium  (cromo) content (current) - CHCON
      ALPHAM(IN+23)= QO         ! molybdenum (molibdeno) content (current) - QOCON
      ALPHAM(IN+24)= QI         ! nickel (nikel) content (current) - QICON
      ALPHAM(IN+25)= 0.0D0      ! not used yet 1
      ALPHAM(IN+26)= 0.0D0      ! not used yet 2
C     
      DO INU=27,28+3*(NNUM4T-1) ! NNUM4T different graphite grain
         ALPHAM(IN+INU)= 0.0D0  ! densities (z1/2/3)
      ENDDO
      DO INU=30+3*(NNUM4T-1),
     .     31+5*(NNUM4T-1)      ! NNUM4T different graphite nod. radii 
         ALPHAM(IN+INU)= 0.0D0
      ENDDO
C     
      ALPHAM(IN+32+5*(NNUM4T-1))= 0.0D0 ! graphite nucl. index - jo
      ALPHAM(IN+33+5*(NNUM4T-1))= 0.0D0 ! grap. nucl. index for rec.- INDEXG
      ALPHAM(IN+34+5*(NNUM4T-1))= 0.0D0 ! grap. nucl. index for T_min - SINDEXG
      ALPHAM(IN+35+5*(NNUM4T-1))= 0.0D0 ! grap. Rg growth index INDEXR=1
      ALPHAM(IN+36+5*(NNUM4T-1))= 0.0D0 ! local final time of solidification
C     
C**** INCREMENTS "ALPHAM" INDEX
C     
      INUPC= INUPC+26+5*NNUM4T+5 ! less than NBASES; see pointes.f
C     
      RETURN
      END
      
