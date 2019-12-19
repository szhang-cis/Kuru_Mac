      SUBROUTINE INMICRO12(ALPHAM,INUPC,IPLAT)
C***********************************************************************
C     
C**** THIS ROUTINE INITIALITES MICROSTRUCTURAL PARAMETERS FOR MODEL 12
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
C     
      TEFER = VPLAT(IPLAT, 6)   ! stable temperature (ferrite)
      TEPER = VPLAT(IPLAT, 7)   ! metastable temperature (perlite)
      HENERF= VPLAT(IPLAT, 8)   ! latent heat of ferrite (eutectoid)
      HENERP= VPLAT(IPLAT, 9)   ! latent heat of pearlite (eutectoid)
      GNF   = VPLAT(IPLAT,10)   ! initial numbers of ferrite grains
      RINIFE= VPLAT(IPLAT,11)   ! initial radius of ferrite grains
C     
C***  ferrite growth
C     
      IGRMFX= INT(VPLAT(IPLAT,12))
      IF(IGRMFX.EQ.1.OR.IGRMFX.EQ.2.OR.IGRMFX.EQ.3) THEN
         DENSF = VPLAT(IPLAT,13) ! ferrite density
         DENSA = VPLAT(IPLAT,14) ! austenite density
         DIFCA = VPLAT(IPLAT,15) ! diffusion coef. C in austenite
         HENERF= VPLAT(IPLAT,16) ! latent heat of ferrite (eutectoid)
      ENDIF
C     
C***  pearlite nucleation
C     
      INUCMPX= INT(VPLAT(IPLAT,17))
C     lacaze, own or indues model
      IF(INUCMPX.EQ.1.OR.INUCMPX.EQ.2.OR.INUCMPX.EQ.3) THEN 
         GNUCCP = VPLAT(IPLAT,18) ! nucleation coefficient of pearlite
         GNUCEP = VPLAT(IPLAT,19) ! nucleation exponent of pearlite
         IGNUCOP= INT(VPLAT(IPLAT,20)) ! no assigned
      ENDIF
      INUCAX= INT(VPLAT(IPLAT,21)) ! nucleation arrest criterion index
C     
C***  pearlite growth
C     
      IGROMPX= INT(VPLAT(IPLAT,22)) ! growth model index
C     lacaze or indues model
      IF(IGROMPX.EQ.1.OR.IGROMPX.EQ.2.OR.IGROMPX.EQ.3) THEN
         GROCP = VPLAT(IPLAT,23) ! growth coefficient of pearlite
         GROEP = VPLAT(IPLAT,24) ! growth exponent of pearlite
         HENERP= VPLAT(IPLAT,25) ! latent heat of pearlite (eutectoid)
         IGROMP= INT(VPLAT(IPLAT,26)) ! pearlite growth model 
      ENDIF
C     
      IKMICX= INT(VPLAT(IPLAT,27)) ! index for micro.-dep. conductivity
      IKAUX = 0
      IF(IKMICX.EQ.1.OR.IKMICX.EQ.2.OR.IKMICX.EQ.3) THEN
         IKAUX= 3
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
C**** INITIALITES MICROSTRUCTURAL VARIABLES ("ALPHAM")
C     
      IN= INUPC
      ALPHAM(IN+1)= 0.0D0       ! ferrite fraction
      ALPHAM(IN+2)= 0.0D0       ! pearlite fraction
      ALPHAM(IN+3)= 0.0D0       ! carbon content in austenite in contact with graphite
      ALPHAM(IN+4)= 0.0D0       ! carbon content in austenite (mass equilibrium)
C     
C***  at the first time, we assume that boundary layer take the typicall
C     value published in 
C     Stefanescus's Solidification Book (chapter 2 pp.21)      
C     DELTAG
      DO INU=5,5+(NNUM4T-1)     ! NNUM4T different ferrite grain radius
         ALPHAM(IN+INU)= 5.0D-10 ! (micro m) == 5e-9 m (not used now)
      ENDDO        
C***  DELTAF
      DO INU=6+(NNUM4T-1),      ! NNUM4T different ferrite grain radius      
     .     6+2*(NNUM4T-1)       ! (micro m) == 5e-9 m 
         ALPHAM(IN+INU)= 5.0D-10
      ENDDO        
C***  RRAFGGN
      DO INU=7+2*(NNUM4T-1),    ! NNUM4T different ferrite grain radius
     .     7+3*(NNUM4T-1)
         ALPHAM(IN+INU)= 0.0D0  
      ENDDO        
C***  FNU
      DO INU=8+3*(NNUM4T-1),      
     .     8+4*(NNUM4T-1)       ! NNUM4T different pearlite grain densities
         ALPHAM(IN+INU)= 1.0D-7 ! it is set in micros12.f
      ENDDO        
C***  DPERC
      DO JNU=9+4*(NNUM4T-1),
     .     9+5*(NNUM4T-1)       ! NNUM4T different pearlite grain densities
         ALPHAM(IN+JNU)= 0.0D0
      ENDDO
C***  PNU
      DO JNU=10+5*(NNUM4T-1),
     .     10+6*(NNUM4T-1)      ! NNUM4T different pearlite grain radius
         ALPHAM(IN+JNU)= 0.0D0
      ENDDO
C     carbon content corresponding to S' point in Fe-C-Si diagram
      ALPHAM(IN+11+6*(NNUM4T-1))= 0.0D0 
C     initial fraction of austenite
      ALPHAM(IN+12+6*(NNUM4T-1))= 0.0D0
C     initial relation of austenite fraction transformed
      ALPHAM(IN+13+6*(NNUM4T-1))= 0.0D0 
C     initial pearlites colonies densities (for instantaneus nucleation)
      ALPHAM(IN+14+6*(NNUM4T-1))= 0.0D0 
C     pearlite nucl. indexC (pearlite nucleation index)
      ALPHAM(IN+15+6*(NNUM4T-1))= 0.0D0 
C     index to set TSP and CAFO
      ALPHAM(IN+16+6*(NNUM4T-1))= 0.0D0 
C     temperature corresponding to S prime point in FE-C-Si phase diagram
      ALPHAM(IN+17+6*(NNUM4T-1))= 0.0D0 
C     temperature corresponding to S prime point in FE-C-Si phase diagram
      ALPHAM(IN+18+6*(NNUM4T-1))= 0.0D0 
C     graphite volumetric fraction
C     ALPHAM(IN+19+6*(NNUM4T-1))= 0.0D0
C     index to set RFA for eutectoid transformation
      ALPHAM(IN+19+6*(NNUM4T-1))= 0.0D0
C     index to set CFG
      ALPHAM(IN+20+6*(NNUM4T-1))= 0.0D0
C     
      IF(INUCAX.EQ.1)
     .     ALPHAM(IN+21+6*(NNUM4T-1))= 0.0D0 ! index: nucleation up to Trecal
      IF(INUCAX.EQ.2)
     .     ALPHAM(IN+21+6*(NNUM4T-1))= 2.0D5 ! index: nucleation up to Tmin
C     
C**** INCREMENTS "ALPHAM" INDEX
C     
      INUPC= INUPC+4+6*NNUM4T+11 ! less than NBASES; see pointes.f
C     
      RETURN
      END
      
