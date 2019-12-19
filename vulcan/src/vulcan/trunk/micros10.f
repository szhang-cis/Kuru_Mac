      SUBROUTINE MICROS10(TGAUAT,TGAUST,TGAUIT,DTEMPT,
     .                     BASMM, BASCC, BASKK,
     .                     TSOE2, TSOE1, TSOC1,
     .                     IPLAT,
     .                    ALPHAM, INUPC)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE PHASE-CHANGE FUNCTION ACCORDING TO THE
C     MICROSTRUCTURAL MODEL NUMBER 10 (IPCMO=10)
C
C     COPPER (CONTINUOUS NUCLEATION AND GROWTH) MICROSTRUCTURAL MODEL
C
C     Notes:
C     THIS ROUTINE UTILIZE KELVIN DEGREES FOR TEMPERATURE
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
C     ALPHAM(IN+1)=liquid fraction
C     ALPHAM(IN+2)=dendritic solid fraction 
C     ALPHAM(IN+3)=eutectic solid fraction
C     ALPHAM(IN+4)=nuclei density
C     ALPHAM(IN+5)=radii of nuclei
C     ALPHAM(IN+6)=nucleation index
C     ALPHAM(IN+7)=dendritic internal solid fraction
C     ALPHAM(IN+8)=eutectic internal solid fraction
C     ALPHAM(IN+9)=radial globular fraction (FSOPC)
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
C**** CHECKS LIQUID FRACTION
C
      IN=INUPC
      FLIQD=ALPHAM(IN+1)       ! liquid fraction
C
      IF(FLIQD.EQ.0.0D0) THEN  ! gets out if liquid fraction is zero
       INUPC=INUPC+12
       RETURN
      ENDIF
C
      TWOPIT=6.283185307179586D0
      PI=TWOPIT/2.0D0
C
C**** TRANSFERS "VPLAT" TO MICROSCOPICAL VARIABLES
C
      TEINF=VPLAT(IPLAT,1)     ! not used
      TESUP=VPLAT(IPLAT,2)     ! not used
      HENER=VPLAT(IPLAT,3)
      DELEE=TESUP-TEINF        ! not used
C
      CTOK=VPLAT(IPLAT,6)
C
      OX =VPLAT(IPLAT, 7)
      PH =VPLAT(IPLAT, 8)
      AS =VPLAT(IPLAT, 9)
      SB =VPLAT(IPLAT,10)
      FE =VPLAT(IPLAT,11)
      SE =VPLAT(IPLAT,12)
      TE =VPLAT(IPLAT,13)
      SU =VPLAT(IPLAT,14)
      ZN =VPLAT(IPLAT,15)
      PB =VPLAT(IPLAT,16)
C
C     TLI=VPLAT(IPLAT,17)          ! Only to test
      TLD=VPLAT(IPLAT,17)             
      EU =VPLAT(IPLAT,18)
      OP =VPLAT(IPLAT,19)
C
      INUCMX=INT(VPLAT(IPLAT,20))  ! nucleation model index
      IF(INUCMX.EQ.3) THEN         ! Oldfield model
       ANUCQ=VPLAT(IPLAT,21)       ! A parameter
       ANUCR=VPLAT(IPLAT,22)       ! n coefficient
       ICX=2
      ENDIF
C
      IGROMX=INT(VPLAT(IPLAT,21+ICX))       ! growth model index
      IF(IGROMX.EQ.2) THEN         ! Pero - Sanz growth model
       GROCF=VPLAT(IPLAT,22+ICX)   ! Z coefficient
       IGX=1
      ENDIF
C
      IFPCDT=INT(VPLAT(IPLAT,22+ICX+IGX))   ! index for temp. derivative
      IAFLOJ=INT(VPLAT(IPLAT,23+ICX+IGX))   ! not used
C
      IV=23+ICX+IGX              ! number of VPLAT defined in input data
C
C**** TRANSFERS (REST OF) "ALPHAM" TO MICROSCOPICAL VARIABLES
C
      FSOLD=ALPHAM(IN+2)           ! dendritic solid fraction
      FSOEU=ALPHAM(IN+3)           ! eutectic solid fraction
      DNA  =ALPHAM(IN+4)           ! nuclei density
      RNA  =ALPHAM(IN+5)           ! radii of nuclei
C
      INDEC=INT(ALPHAM(IN+6))      ! nucleation index                   
      FINTD=ALPHAM(IN+7)           ! dendritic internal solid fraction
      FINTDE=ALPHAM(IN+8)          ! eutectic internal solid fraction
      FSOPC=ALPHAM(IN+9)           ! radial globular fraction (FSOPC)
      OPE=ALPHAM(IN+10)            ! Effective part. coeff
      TSOCCP=ALPHAM(IN+11)         ! solo prueba, Primary phase
      TSOCCE=ALPHAM(IN+12)         ! solo prueba, Eutectic phase
C
C**** LOAD LAST CONVERGED VALUES
C
      FLIQDO=FLIQD
      FSOLID=FSOLIDO
      FSOPCO=FSOPC
      FSOEUO=FSOEU
      FSOLDO=FSOLD
      FINTDO=FINTD
      FINTDEO=FINTDE
      FSOPCO=FSOPC
C
C**** COMPUTES EQUILIBRIUM LIQUIDUS TEMPERATURE & OTHER PARAMETERS
C
C     Nomenclature:
C     TLD: equilibrium liquidus temperature [K]
C     TSD: equilibrium solidus temperature [K]
C     TSS: equilibrium solvus temperature [K]
C     SOLIDOX: equilibrium oxygen concentration for solidus temperature
C     SOLVUX: equilibrium oxygen concentration for solvus temperature
C     TOPUN: maximum undercooling according to the oxygen concentration
C
      TGAUST=TGAUST+CTOK
      TGAUAT=TGAUAT+CTOK
      TGAUIT=TGAUIT+CTOK
C
C     TLD=1356.0D0-(41.782D0*OX)
      IF(TGAUST.LE.1356.0D0.AND.TGAUST.GT.1338.0D0) THEN
       SOLIDOX=(1356.0D0-TGAUST)/2400.0D0
C      TSD=1356.0D0-(2400.0D0*OX)     ! valid only for %O < 0.0075%
      END IF
C
      IF(TGAUST.LT.1338.0D0) THEN
       SOLVUX=(TGAUST-298.0D0)/138666.0D0
C      TSS=298.0D0+(138666.67D0*OX)   ! valid only for %O < 0.0075%
      END IF
C
C**** ESTABLISHES VARIABLES OF MODEL 12
C
      FSOLID=1.0D0-FLIQD        ! solid fraction
      FSOLIDO=1.0D0-FLIQDO
C
C**** SOLVES MODEL 12
C
C*****SCHEIL SOLIDIFICACTION EXTREME CASE, ONLY TO TEST*****************
C
C     PCU=TLD-TGAUST
C     DIFOL=6.51836E-27*TGAUST**5.74954D0  ! diffusivity of O in liq. Cu
C     IF(RNA.EQ.0.0D0) THEN
C      CAP=0.0D0
C      OPE=OP
C     ELSE
C      CAP=2.0D0*DIFOL/(GROCF*PCU/TGAUST)
C      OPE=OP/(OP+(1.0D0-OP)*EXP(-RNA*CAP/DIFOL)) ! Eff. partition coef.
C     END IF
C     CLIQ=OX*((1-FSOLID)**(OPE-1.0D0))
C     CSOL=OPE*OX*((1-FSOLID)**(OPE-1.0D0))
C     TLD=TLI-CLIQ*41.782      ! Initial liquidus temp., TLI; case equil
C     TLD=TLI-(CLIQ-OX)*41.782 ! Initial liquidus temp., TLI; case norma
C     EU=EUI-(CLIQ-OX)*41.782  ! Initial liquidus temp., TL
C***********************************************************************
C
      TSOCCP=0.0D0             ! only useful to thermal jacobian matrix
      TSOCCE=0.0D0
C
C**** PRIMARY COPPER SOLIDIFICATION
C
      ICORREC=0                           ! not used
      IF(RNA.LE.0.0D0) THEN
       CONLIQ=0.0D0
       DIFOL=0.0D0
       ZETAX=0.0D0
       SOBQUIM=0.0D0
c      FLIQD=1.0D0
c      FSOEU=0.0D0
c      FINTD=0.0D0
c      FINTDE=0.0D0
      ELSE 
       CONLIQ=(TGAUST-1356.0D0)/-41.782D0  ! oxygen conc. at liquidus T
       DIFOL=6.51836E-27*TGAUST**5.74954D0 ! diffusivity of O in liq. Cu
      END IF
C     FINTX=FINTD+FINTDE
C
C**** NUCLEATION LAWS
C
C     Notes:
C     TSOCN=-derivative of nucleation law with respect to temperature
C
      PCU=TLD-TGAUST                       ! primary copper undercooling    
      IF(OX.GT.0.0075D0) THEN              ! anodic copper (ox<0.0075)
       IF(PCU.GT.0.0D0.AND.FSOLID.LT.1.0D0) THEN
        IF(DNA.GT.0.0D0.AND.DTEMPT.GT.0.0D0) INDEC=1      ! end of nucl.
        IF(INDEC.EQ.0) THEN
         IF(INUCMX.EQ.3) THEN              ! Oldfield model
          DNAO=DNA
          DNA=ANUCQ*(PCU**ANUCR)*DTIMET    ! nucleation, instantaneous
          TSOCN=ANUCQ*ANUCR*(PCU**(ANUCR-1.0D0))*DTIMET
C
          IF(DNA.LT.DNAO) THEN             ! control
           DNA=DNAO
           TSOCN=0.0D0
          ENDIF
C
          FSOPC=2.0D0/3.0D0*TWOPIT*DNA*RNA**3.0D0
C         IF (FSOPC.LT.FSOPCO) FSOPC=FSOPCO          ! No remelting
C
          IF(FSOPC.GT.1.0D0) THEN        ! assumption valid for small dt
           FSOPC=1.0D0       
           DNA=FSOPC/(2.0D0/3.0D0*TWOPIT*RNA**3.0D0)
           TSOCN=0.0D0                   ! assumption
           FSOPC=2.0D0/3.0D0*TWOPIT*DNA*RNA**3.0D0
          ENDIF
         ENDIF      ! inucmx=3
        END IF      ! indec=0
C
C**** GROWTH LAWS
C
C     Notes:
C     TSOCR=-derivative of growth law with respect to temperature
C
        IF(IGROMX.EQ.2) THEN             ! Del Pero - Sanz model
         RNAO=RNA
         RNA=RNA+(GROCF*PCU/TGAUST)*DTIMET
         TSOCR=(GROCF*((TGAUST+PCU)/(TGAUST**2)))*DTIMET
C
         IF(RNA.LT.RNAO) THEN            ! control
          RNA=RNAO
          TSOCR=0.0D0
         ENDIF
C
         FSOPC=2.0D0/3.0D0*TWOPIT*DNA*RNA**3.0D0
C        IF (FSOPC.LT.FSOPCO) FSOPC=FSOPCO          ! No remelting
         IF(FSOPC.GT.1.0D0) THEN         ! assumption valid for small dt
          FSOPC=1.0D0
          RNA=(FSOPC/(2.0D0/3.0D0*TWOPIT*DNA))**(1.0D0/3.0D0)
          TSOCR=0.0D0                    ! assumption
          FSOPC=2.0D0/3.0D0*TWOPIT*DNA*RNA**3.0D0
         ENDIF
C
         SOBQUIM=(CONLIQ-OX)/(CONLIQ*(1.0D0-OP))  ! Chemical supersatur.
         ZETAX=2.0D0*DIFOL/(GROCF*PCU*RNA/TGAUST) ! Fac. of solute diff.
         GZETAX=(1.0D0+(1.5D0*ZETAX)+
     .          (ZETAX**2)+(ZETAX**3/4.0D0)) ! less than 1.0
         IF(GZETAX.GT.1.0D0) GZETAX=1.0D0    ! Control, int.sol.fraction
         FINTD=SOBQUIM*GZETAX                ! Internal solid fraction
         TSOCF=((TGAUST-1356.0D0)**2-(41.782D0*OX))*GZETAX/
     .         ((OP-1.0D0)*(TGAUST-1356.0D0)**2.0D0)
C
         IF(FINTD.LT.FINTDO) FINTD=FINTDO    ! control (no remelting)
         FINTX=FINTD+FINTDE                  ! control
         IF(FINTX.GT.1.0D0) THEN         ! assumption valid for small dt
          FINTD=1.0D0-FINTDE
          TSOCF=0.0D0                    ! assumption
         ENDIF
C
         FSOLD=FSOPC*FINTD               ! Dendritic solid fraction
         TSOCCP=2.0D0*TWOPIT/3.0D0*(TSOCN*RNA**3.0D0*FINTD+
     .          DNA*3.0D0*RNA**2.0D0*TSOCR*FINTD+
     .          DNA*RNA**3.0D0*TSOCF)
         FSOLID=FSOLD+FSOEU
         IF(FSOLID.GT.1.0D0) THEN
          FSOLID=1.0D0
          FSOLD=1.0D0-FSOEU
          TSOCCP=0.0D0
         END IF       
        END IF    ! igromx=2
       END IF     ! pcu>0 & fsolid=0
C
C**** EUTECTIC SOLIDIFICATION
C
C     Notes:
C     TSOCE=-derivative of eutectic fraction with respect to temperature
C
       PCE=EU-TGAUST                     ! eutectic undercooling
       IF(PCE.GE.0.0D0.AND.FSOLID.LT.1.0D0) THEN
        FINTDE=PB*OX*(EU-TGAUST)
        IF(FINTDE.LT.FINTDEO) FINTDE=FINTDEO   ! control (no remelting)
        TSOCFE=-PB*OX
        FINTX=FINTD+FINTDE                     ! control
        IF(FINTX.GT.1.0D0) THEN          ! assumption valid for small dt
         FINTDE=1.0D0-FINTD
         TSOCFE=0.0D0
        ENDIF
        FSOEU=FSOPC*FINTDE                     ! Eutectic solid fraction
        TSOCCE=2.0D0*TWOPIT/3.0D0*(TSOCN*RNA**3.0D0*FINTDE+
     .         DNA*3.0D0*RNA**2.0D0*TSOCR*FINTDE+
     .         DNA*RNA**3.0D0*TSOCFE)
        FSOLID=FSOLD+FSOEU
        IF(FSOLID.GT.1.0D0) THEN
         FSOLID=1.0D0
         FSOEU=1.0D0-FSOLD
        END IF
       END IF
C
C**** GLOBAL OUTPUTS: LIQUID FRACTION
C
       FLIQD=1.000D0-FSOLID                 ! liquid fraction
       IF(FLIQD.LT.0.0D0) FLIQD=0.0D0       ! control
       IF(FLIQD.GT.1.0D0) FLIQD=1.000D0     ! control
C
      ELSE          ! electrolitic copper (ox<0.0075)
C
       IF(PCU.GT.0.0D0.AND.FSOPC.LT.1.0D0.AND.FSOLID.LT.1.0D0) THEN
        IF(DNA.GT.0.0D0.AND.DTEMPT.GT.0.0D0) INDEC=1      ! end of nucl.
        IF(INDEC.EQ.0) THEN
         IF(INUCMX.EQ.3) THEN            ! Oldfield model
          DNAO=DNA
          DNA=ANUCQ*(PCU**ANUCR)*DTIMET  ! nucleation, instantaneous
          TSOCN=ANUCQ*ANUCR*(PCU**(ANUCR-1.0D0))*DTIMET
C
          IF(DNA.LT.DNAO) THEN           ! control
           DNA=DNAO
           TSOCN=0.0D0
          ENDIF
C
          FSOPC=2.0D0/3.0D0*TWOPIT*DNA*RNA**3.0D0
          IF(FSOPC.GT.1.0D0) THEN        ! assumption valid for small dt
           FSOPC=1.0D0       
           DNA=FSOPC/(2.0D0/3.0D0*TWOPIT*RNA**3.0D0)
           TSOCN=0.0D0                   ! assumption
           FSOPC=2.0D0/3.0D0*TWOPIT*DNA*RNA**3.0D0
          ENDIF
         ENDIF      ! inucmx=3
        ENDIF       ! indec=0
C
C**** GROWTH LAW
C
        IF(IGROMX.EQ.2) THEN             ! Del Pero - Sanz model
         RNAO=RNA  
         RNA=RNA+(GROCF*PCU/TGAUST)*DTIMET
         TSOCR=(GROCF*((TGAUST+PCU)/(TGAUST**2.0D0)))*DTIMET
C
         IF(RNA.LT.RNAO) THEN            ! control
          RNA=RNAO
          TSOCR=0.0D0
         ENDIF
C
         FSOPC=2.0D0/3.0D0*TWOPIT*DNA*RNA**3.0D0
         IF(FSOPC.GT.1.0D0) THEN         ! assumption valid for small dt
          FSOPC=1.0D0
          RNA=(FSOPC/(2.0D0/3.0D0*TWOPIT*DNA))**(1.0D0/3.0D0)
          TSOCR=0.0D0                    ! assumption
          FSOPC=2.0D0/3.0D0*TWOPIT*DNA*RNA**3.0D0
         ENDIF
C
         SOBQUIM=(CONLIQ-OX)/(CONLIQ*(1.0D0-OP))  ! Chemical supersat.
         ZETAX=2.0D0*DIFOL/(GROCF*PCU*RNA/TGAUST) ! Fac. of solute diff.
         GZETAX=(1.0D0+(1.5D0*ZETAX)+
     .          (ZETAX**2)+(ZETAX**3/4.0D0))      ! less than 1.0
         IF(GZETAX.GT.1.0D0) GZETAX=1.0D0    ! Control, int.sol.fraction
         FINTD=SOBQUIM*GZETAX                ! Internal solid fraction
         TSOCF=((TGAUST-1356.0D0)**2-(41.782D0*OX))*GZETAX/
     .         ((OP-1.0D0)*(TGAUST-1356.0D0)**2.0D0)
C
         IF(FINTD.LT.FINTDO) THEN 
          FINTD=FINTDO                   ! control (no remelting)
          TSOCF=0.0D0                    ! assumption
         ENDIF
C
         FSOLD=FSOPC*FINTD               ! Dendritic solid fraction
         TSOCCP=2.0D0*TWOPIT/3.0D0*(TSOCN*RNA**3.0D0*FINTD+
     .          DNA*3.0D0*RNA**2.0D0*TSOCR*FINTD+
     .          DNA*RNA**3.0D0*TSOCF)
         FSOLID=FSOLD          
         IF(FSOLID.GT.1.0D0) THEN
          FSOLID=1.0D0
          TSOCCP=0.0D0
         END IF
        END IF      ! igromx=2
C
C**** GLOBAL OUTPUTS
C
        FLIQD=1.0D0-FSOLID               ! liquid fraction
        IF(FLIQD.LT.0.0D0) FLIQD=0.0D0   ! control
        IF(FLIQD.GT.1.0D0) FLIQD=1.0D0   ! control
       END IF       ! pcu>0,...
      END IF        ! ox>0.0075
C
C**** ESTABLISHES LATENT HEAT * PHASE-CHANGE FUNCTION (at time t & t+dt)
C
      TSOE1(IPLAT)=FLIQD                 ! f_pc at time t+dt
      TSOE2(IPLAT)=FLIQDO                ! f_pc at time t
      TSOE2(IPLAT)=TSOE2(IPLAT)*HENER    ! L*f_pc at time t
      TSOE1(IPLAT)=TSOE1(IPLAT)*HENER    ! L*f_pc at time t+dt
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
      ENDIF
      GO TO 10
C
    2 ICACOT=1
      VELTOT=1.D-10
      VELA1T=DABS(DTEMPT)
      IF(VELA1T.GT.VELTOT) THEN
       IF(IITERT.GT.0)
     .  TSOC1(IPLAT)=(TSOE1(IPLAT)-TSOE2(IPLAT))/(DTEMPT*DTIMET)
      ENDIF
      GO TO 10
C
    3 ICACOT=0
      TSOC1(IPLAT)=TSOCCP+TSOCCE
      TSOC1(IPLAT)=TSOC1(IPLAT)*HENER
      GO TO 10
C
    4 ICACOT=1
      TSOC1(IPLAT)=TSOCPC+TSOCCP+TSOCCE
      TSOC1(IPLAT)=TSOC1(IPLAT)*HENER
      GO TO 10
C
    5 ICACOT=1
      TSOC1(IPLAT)=0.0D0
      GO TO 10
C
   10 CONTINUE
C
C**** TRANSFER MICROSCOPICAL VARIABLES TO "ALPHAM"
C
      ALPHAM(IN+1)=FLIQD
      ALPHAM(IN+2)=FSOLD
      ALPHAM(IN+3)=FSOEU
      ALPHAM(IN+4)=DNA
      ALPHAM(IN+5)=RNA
C
      ALPHAM(IN+6)=FLOAT(INDEC)
      ALPHAM(IN+7)=FINTD 
      ALPHAM(IN+8)=FINTDE
      ALPHAM(IN+9)=FSOPC  
      ALPHAM(IN+10)=OPE   
      ALPHAM(IN+11)=TSOCCP
      ALPHAM(IN+12)=TSOCCE
C
      TGAUST=TGAUST-CTOK
      TGAUAT=TGAUAT-CTOK
      TGAUIT=TGAUIT-CTOK
C
C**** INCREMENTS "ALPHAM" INDEX
C
      INUPC=INUPC+12                   ! less than NBASES; see pointes.f
C
      RETURN
      END
