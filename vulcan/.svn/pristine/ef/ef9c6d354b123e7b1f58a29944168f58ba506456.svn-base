      SUBROUTINE IDEPROS12(PROPST,IPLAT,IA1)
C***********************************************************************
C     
C**** THIS ROUTINE ORDERS THE MICROSTRUCTURAL PROPERTIES OF MODEL
C     NUMBER 12 (IPCMO=12) OF RATE PHASE-CHANGE FORMULATIONS
C     
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C     
C**** THERMAL VARIABLES
C     
      INCLUDE 'auxl_omt.f'
      INCLUDE 'prob_omt.f'
C     
      DIMENSION PROPST(*)
C     
      IA2=IA1+2                 ! 2=ipcfo,ipcmo
C     
      TEFER = PROPST(IA2+1)     ! stable temperature (ferrite)
      TEPER = PROPST(IA2+2)     ! metastable temperature (perlite)
      HENERF= PROPST(IA2+3)     ! latent heat of ferrite (eutectoid)
      HENERP= PROPST(IA2+4)     ! latent heat of pearlite (eutectoid)
      GNF   = PROPST(IA2+5)     ! kinetics params of ferrite
      RINIFE= PROPST(IA2+6)     ! diffusion coef. C in ferrite
C     
C***  ferrite growth
C     
      IGRMFX= INT(PROPST(IA2+7)) ! ferrite's growth model
      IF(IGRMFX.EQ.1.OR.IGRMFX.EQ.2.OR.IGRMFX.EQ.3) THEN
         DENSF = PROPST(IA2+ 8) ! ferrite density
         DENSA = PROPST(IA2+ 9) ! austenite density
         DIFCA = PROPST(IA2+10) ! diffusion coef. C in austenite
         HENERF= PROPST(IA2+11) ! latent heat of ferrite (eutectoid)
      ENDIF
C     pearlite nucleation
      INUCMPX= INT(PROPST(IA2+12))
      IF(INUCMPX.EQ.1.OR.INUCMPX.EQ.2.OR.INUCMPX.EQ.3) THEN 
         GNUCCP = PROPST(IA2+13) ! nucleation coefficient of pearlite
         GNUCEP = PROPST(IA2+14) ! nucleation exponent of pearlite
         IGNUCOP= INT(PROPST(IA2+15)) ! continuos or instantaneous nucleation
      ENDIF
C     pearlite arrest criterion
      INUCAX= INT(PROPST(IA2+16))
C     pearlite growth
      IGROMPX=INT(PROPST(IA2+17)) ! pearlite's growth model index
      IF(IGROMPX.EQ.1.OR.IGROMPX.EQ.2.OR.IGROMPX.EQ.3) THEN
         GROCP = PROPST(IA2+18) ! growth coefficient of pearlite
         GROEP = PROPST(IA2+19) ! growth exponent of pearlite
         HENERP= PROPST(IA2+20) ! latent heat of pearlite (eutectoid)
         IGROMP= INT(PROPST(IA2+21)) ! pearlite's growth model 
      ENDIF
C     
      IKMICX= INT(PROPST(IA2+22))
      IKAUX = 0
      IF(IKMICX.EQ.1.OR.IKMICX.EQ.2.OR.IKMICX.EQ.3) THEN
         IKAUX= 3
         BASKS= PROPST(IA2+23)
         BASKM= PROPST(IA2+24)
         BASKL= PROPST(IA2+25)
      ENDIF
C     
      IFPCDT= INT(PROPST(IA2+23+IKAUX))
      IAFLOJ= INT(PROPST(IA2+24+IKAUX))
C     
      VPLAT(IPLAT, 6)= TEFER    ! stable temperature (ferrite)
      VPLAT(IPLAT, 7)= TEPER    ! metastable temperature (perlite)
      VPLAT(IPLAT, 8)= HENERF   ! latent heat of ferrite (eutectoid)
      VPLAT(IPLAT, 9)= HENERP   ! latent heat of pearlite (eutectoid)
C     
C***  graphite growth
C     
      VPLAT(IPLAT,10)= GNF      ! kinetics params of ferrite
      VPLAT(IPLAT,11)= RINIFE   ! diffusion coef. C in ferrite
C     
C***  ferrite growth
C     
      VPLAT(IPLAT,12)= FLOAT(IGRMFX) ! ferrite's growth model
      IF(IGRMFX.EQ.1.OR.IGRMFX.EQ.2.OR.IGRMFX.EQ.3) THEN
         VPLAT(IPLAT,13)= DENSF ! ferrite density
         VPLAT(IPLAT,14)= DENSA ! austenite density
         VPLAT(IPLAT,15)= DIFCA ! diffusion coef. C in austenite
         VPLAT(IPLAT,16)= HENERF ! latent heat of ferrite (eutectoid)
      ENDIF
C     pearlite nucleation
      VPLAT(IPLAT,17)= FLOAT(INUCMPX) ! pearlite nucleation model
      IF(INUCMPX.EQ.1.OR.INUCMPX.EQ.2.OR.INUCMPX.EQ.3) THEN
         VPLAT(IPLAT,18)= GNUCCP ! nucleation coefficient of pearlite
         VPLAT(IPLAT,19)= GNUCEP ! nucleation exponent of pearlite
         VPLAT(IPLAT,20)= FLOAT(IGNUCOP) ! not assigned yet
      ENDIF
C     arrest criterion
      VPLAT(IPLAT,21)= FLOAT(INUCAX) ! arrest criterion
C     pearlite growth
      VPLAT(IPLAT,22)= FLOAT(IGROMPX) !  growth model
      IF(IGROMPX.EQ.1.OR.IGROMPX.EQ.2.OR.IGROMPX.EQ.3) THEN
         VPLAT(IPLAT,23)= GROCP ! growth coefficient of pearlite
         VPLAT(IPLAT,24)= GROEP ! growth exponent of pearlite
         VPLAT(IPLAT,25)= HENERP ! latent heat of pearlite (eutectoid)
         VPLAT(IPLAT,26)= FLOAT(IGROMP) ! not assigned yet
      ENDIF
C     
      VPLAT(IPLAT,27)= FLOAT(IKMICX) ! index for micro.-dep. conductivity
      IKAUX=0
      IF(IKMICX.EQ.1.OR.IKMICX.EQ.2.OR.IKMICX.EQ.3) THEN
         IKAUX=3
         VPLAT(IPLAT,28)= BASKS ! solid conductivity
         VPLAT(IPLAT,29)= BASKM ! mushy conductivity
         VPLAT(IPLAT,30)= BASKL ! liquid conductivity
      ENDIF
C     
      VPLAT(IPLAT,28+IKAUX)= FLOAT(IFPCDT) ! index for temp. derivative
      VPLAT(IPLAT,29+IKAUX)= FLOAT(IAFLOJ) ! index for fraction correction
C     
      IMODE= 24+IKAUX           ! imode= total number of prop. of model 12
      IA1  = IA2+IMODE
C     
      RETURN
      END
      
