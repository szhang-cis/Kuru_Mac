      SUBROUTINE MICROS5(TGAUAT,TGAUST,TGAUIT,DTEMPT,
     .                    BASMM, BASCC, BASKK,
     .                    TSOE2, TSOE1, TSOC1,
     .                    IPLAT,
     .                   ALPHAM,INUPC)  
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE PHASE-CHANGE FUNCTION ACCORDING TO THE
C     MICROSTRUCTURAL MODEL NUMBER 5 (IPCMO=5)
C
C     S.G. CAST IRON MICROSTRUCTURAL MODEL: SU'S MODEL (UNI-NODULAR) FOR
C     EUTECTIC & HYPEREUTECTIC COMPOSITIONS
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
C     ALPHAM(IN+2)=austenite fraction (solid)
C     ALPHAM(IN+3)=graphite fraction (solid)
C     ALPHAM(IN+4)=silicon content (current)
C     NNUM4T=number of different grain densities and radii
C     ALPHAM(IN+5:IN+5+(NNUM4T-1))=graphite grain density
C     ALPHAM(IN+6+(NNUM4T-1):IN+6+2*(NNUM4T-1))=graphite grain radius
C     ALPHAM(IN+7+2*(NNUM4T-1):IN+7+3*(NNUM4T-1))=austenite grain radius
C
C     Auxiliar microstructural variables (not printed; see pointes.f):
C     ALPHAM(IN+8+3*(NNUM4T-1))=graphite nucleation index
C     ALPHAM(IN+9+3*(NNUM4T-1))=grap. nucleation arrest index
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
C**** CHECKS LIQUID FRACTION
C
      IN=INUPC
      FLIQD=ALPHAM(IN+1)
C
      IF(FLIQD.EQ.0.0D0) THEN      ! gets out if liquid fraction is zero
       INUPC=INUPC+4+3*NNUM4T+2
       RETURN
      ENDIF
C
C**** TRANSFERS "VPLAT" TO MICROSCOPICAL VARIABLES
C
      TEINF=VPLAT(IPLAT,1)          ! not used
      TESUP=VPLAT(IPLAT,2)          ! not used
      HENER=VPLAT(IPLAT,3)
      DELEE=TESUP-TEINF             ! not used
C
      CA =VPLAT(IPLAT, 6)
      SIO=VPLAT(IPLAT, 7)
      PH =VPLAT(IPLAT, 8)
      CU =VPLAT(IPLAT, 9)
      QM =VPLAT(IPLAT,10)
      QG =VPLAT(IPLAT,11)
C
      AKSI=VPLAT(IPLAT,12)
      IHELAC=INT(VPLAT(IPLAT,13))
C
      INUCMX=INT(VPLAT(IPLAT,14))  ! nucleation model index
      IF(INUCMX.EQ.1.OR.           ! Su, Boeri & Lacaze nuc. models
     .   INUCMX.EQ.2.OR.INUCMX.EQ.3) THEN
       ANUCB=VPLAT(IPLAT,15)       ! nucleation coeff. A
       ANUCC=VPLAT(IPLAT,16)       ! nucleation coeff. n
      ENDIF
      INUCAX=INT(VPLAT(IPLAT,17))  ! nucleation arrest criterion index
      INUCON=INT(VPLAT(IPLAT,18))  ! control minimum N index
      ANUCON=VPLAT(IPLAT,19)       ! minimum N
C
      IGROMX=INT(VPLAT(IPLAT,20))  ! growth model index
      IF(IGROMX.EQ.1.OR.           ! Su (2) & Lacaze growth models
     .   IGROMX.EQ.2.OR.IGROMX.EQ.3) THEN
       DIFCA=VPLAT(IPLAT,21)       ! diffusion coef. of C in austenite
       RSHEL=VPLAT(IPLAT,22)       ! initial r of austenite shell
       RNODO=VPLAT(IPLAT,23)       ! initial radius of nodules
       AUSGR=VPLAT(IPLAT,24)       ! auste. shell radius/graphite radius
      ENDIF
C
      IKMICX=INT(VPLAT(IPLAT,25))  ! index for micro.-dep. conductivity
      IKAUX=0
      IF(IKMICX.EQ.1.OR.IKMICX.EQ.2.OR.IKMICX.EQ.3) THEN
       IKAUX=3
       BASKS=VPLAT(IPLAT,26)       ! solid conductivity
       BASKM=VPLAT(IPLAT,27)       ! mushy conductivity
       BASKL=VPLAT(IPLAT,28)       ! liquid conductivity
      ENDIF
C
      IFPCDT=INT(VPLAT(IPLAT,26+IKAUX))  ! index for temp. derivative
      IAFLOJ=INT(VPLAT(IPLAT,27+IKAUX))  ! not used
C
      IV=27+IKAUX                ! number of VPLAT defined in input data
C
C**** TRANSFERS (REST OF) "ALPHAM" TO MICROSCOPICAL VARIABLES
C
      FSOLA=ALPHAM(IN+2)
      FSOLG=ALPHAM(IN+3)
C
      SICON=ALPHAM(IN+4)
C
      JO=INT(ALPHAM(IN+8+3*(NNUM4T-1)))
      IF(JO.GT.0) THEN
       DO INU=1,JO
        DNU(INU,1)=ALPHAM(IN+5+INU-1)
        RNU(INU)  =ALPHAM(IN+6+NNUM4T-1+INU-1)
        ANU(INU)  =ALPHAM(IN+7+2*(NNUM4T-1)+INU-1)
       ENDDO
      ENDIF
C
      IF(INUCAX.EQ.1)
     . INDEXG=INT(ALPHAM(IN+9+3*(NNUM4T-1)))  ! lt NBASES, see pointes.f
      IF(INUCAX.EQ.2)
     . SINDEXG=ALPHAM(IN+9+3*(NNUM4T-1))      ! lt NBASES, see pointes.f
C
C**** LOAD LAST CONVERGED VALUES
C
      FLIQDO=FLIQD
      FSOLAO=FSOLA
C
C**** COMPUTES EQUILIBRIUM EUTECTIC TEMPERATURE & OTHER PARAMETERS
C
C     Nomenclature:
C     TEG: equilibrium graphite eutectic temperature
C     TGL: equilibrium graphite liquidus temperature
C     CLA: liquid %C at the liquid-austenite interface
C     CAL: austenite %C at the liquid-austenite interface
C     CLG: liquid %C at the liquid-graphite interface
C     CAG: austenite %C at the austenite-graphite interface
C
      TGL=389.1D0*(CA+0.33D0*SICON)-503.2D0
C
      IF(IHELAC.EQ.1) THEN                               ! Heine
       TL= 1569.00D0- 97.3D0*(CA+0.25D0*SICON)           ! hypoeutectic
c      TGL=389.10D0*(CA+0.31D0*SICON+0.33D0*PH-1.3D0)    ! hypereutectic
       TGL=389.10D0*(CA+SICON/3.0D0)-503.20D0            ! other option
c      TGL=389.10D0*(CA+0.31D0*SICON)-505.80D0           ! other option
       TEG=1154.6D0+6.5D0*SICON
C
       CE=4.26D0-0.317D0*SICON
c      CEQ=CA+SICON/3.0D0                    ! equivalent carbon content
       CEQ=CA+SIO/3.0D0                      ! equivalent carbon content
       CEQREF=4.26D0                         ! ref. equiv. C content
C
       GUG=TEG-TGAUST                   ! graphite eutectic undercooling
       IF(CEQ.GT.CEQREF) GUG=TGL-TGAUST ! graphite liquidus undercooling
C
       REAG=(26.26D0+0.087D0*SICON)/(2.16D0-0.101D0*SICON)
C
       CLA=(1569.0D0-TGAUST- 24.32D0*SICON)/ 97.3D0
       CAL=(1528.4D0-TGAUST- 32.00D0*SICON)/177.9D0
       CLG=( 503.2D0+TGAUST-129.70D0*SICON)/389.1D0
       CAG=(-1154.6D0+TGAUST-  6.50D0*SICON)*(1.5D0-0.216D0*SICON)/
     .                         (354.6D0+6.5D0*SICON)+2.1D0-0.216D0*SICON
      ENDIF
      IF(IHELAC.EQ.2) THEN                               ! Lacaze
       TL= 1576.00D0- 97.3D0*(CA+0.23D0*SICON)           ! hypoeutectic
c      TGL=389.10D0*(CA+0.31D0*SICON+0.33D0*PH-1.3D0)    ! hypereutectic
       TGL=389.10D0*(CA+0.29D0*SICON)-534.50D0           ! other option
c      TGL=389.10D0*(CA+0.31D0*SICON)-505.80D0           ! other option
       TEG=1153.7D0+4.865D0*SICON
C
       CE=4.34D0-0.28D0*SICON
c      CEQ=CA+0.28D0*SICON                   ! equivalent carbon content
       CEQ=CA+0.28D0*SIO                     ! equivalent carbon content
       CEQREF=4.33D0                         ! ref. equiv. C content
C
       GUG=TEG-TGAUST                   ! graphite eutectic undercooling
       IF(CEQ.GT.CEQREF) GUG=TGL-TGAUST ! graphite liquidus undercooling
C
       REAG=(26.26D0+0.087D0*SICON)/(2.16D0-0.101D0*SICON)
C
       CLA=(1576.0D0-TGAUST- 22.38D0*SICON)/ 97.3D0
       CAL=(1528.4D0-TGAUST- 32.02D0*SICON)/177.9D0
       CLG=( 534.5D0+TGAUST-112.84D0*SICON)/389.1D0
       CAG=(-1153.7D0+TGAUST-  4.8650D0*SICON)*(1.51D0-0.21D0*SICON)/
     .                       (353.7D0+4.865D0*SICON)+2.11D0-0.21D0*SICON
      ENDIF
C
C**** ESTABLISHES VARIABLES OF MODEL 5
C
      FSOLID=1.0D0-FLIQD        ! solid fraction
      FSOLIDO=1.0D0-FLIQDO
C
C**** SOLVES MODEL 5
C
      TSOCGN=0.0D0              ! only useful to thermal jacobian matrix
      TSOCGG=0.0D0              ! not computed
      TSOCGA=0.0D0
C
C**** GRAPHITE PRO & EUTECTIC SOLIDIFICATION
C
      ICORREC=0                         ! not used
C
      IF(GUG.GT.0.0D0.AND.FSOLG.LT.1.0D0.AND.FSOLID.LT.1.0D0) THEN
       TSOCGNN=0.0D0
       IF(INUCAX.EQ.1) THEN             ! nucleation up to Trecal
        IF(JO.GT.0.AND.DTEMPT.GT.0.0D0) INDEXG=1          ! end of nucl.
        IF(INDEXG.EQ.0) THEN            ! nucleation
         JO=JO+1
         IF(JO.GT.NNUM4T)               ! see pointes.f
     .    CALL RUNENDT('ERROR IN MICROS5: JO GT NNUM4T')
         IF(INUCMX.EQ.1) THEN           ! Su model
          DNU(JO,1)=ANUCB*GUG**ANUCC*DTIMET
          TSOCGNN=ANUCB*ANUCC*GUG**(ANUCC-1.0D0)*
     .            DTIMET                ! temperature derivative of -DNU
         ENDIF
         IF(INUCMX.EQ.2) THEN           ! Boeri model
          DNU(JO,1)=ANUCB*GUG*DEXP(-ANUCC/GUG)*FLIQD*DTIMET
          TSOCGNN=ANUCB*DEXP(-ANUCC/GUG)*FLIQD*DTIMET*
     .            (1.0D0+ANUCC/GUG)     ! temperature derivative of -DNU
         ENDIF
         IF(INUCMX.EQ.3) THEN           ! Lacaze model
          DNU(JO,1)=ANUCB*GUG**ANUCC*FLIQD*DTIMET
          TSOCGNN=ANUCB*ANUCC*GUG**(ANUCC-1.0D0)*FLIQD*
     .            DTIMET                ! temperature derivative of -DNU
         ENDIF
         RNU(JO)=RNODO
         ANU(JO)=RSHEL
         IF(INUCON.EQ.1) THEN           ! control minimum N
          IF(DNU(JO,1).LT.ANUCON) THEN
           DNU(JO,1)=0.0D0
           RNU(JO)=0.0D0
           ANU(JO)=0.0D0
           JO=JO-1
          ENDIF
         ENDIF
        ENDIF
       ENDIF                            ! inucax.eq.1
       IF(INUCAX.EQ.2) THEN             ! nucleation up to Tmin
        IF(TGAUST.LT.SINDEXG) SINDEXG=TGAUAT
        IF(DTEMPT.LT.0.0D0.AND.TGAUST.LT.SINDEXG) THEN   ! nucleation
         IF(JO.EQ.0) SINDEXG=TGAUST
         JO=JO+1
         IF(JO.GT.NNUM4T)               ! see pointes.f
     .    CALL RUNENDT('ERROR IN MICROS5: JO GT NNUM4T')
         IF(INUCMX.EQ.1) THEN           ! Su model
          DNU(JO,1)=ANUCB*GUG**ANUCC*DTIMET
          TSOCGNN=ANUCB*ANUCC*GUG**(ANUCC-1.0D0)*
     .            DTIMET                ! temperature derivative of -DNU
         ENDIF
         IF(INUCMX.EQ.2) THEN           ! Boeri model
          DNU(JO,1)=ANUCB*GUG*DEXP(-ANUCC/GUG)*FLIQD*DTIMET
          TSOCGNN=ANUCB*DEXP(-ANUCC/GUG)*FLIQD*DTIMET*
     .            (1.0D0+ANUCC/GUG)     ! temperature derivative of -DNU
         ENDIF
         IF(INUCMX.EQ.3) THEN           ! Lacaze model
          DNU(JO,1)=ANUCB*GUG**ANUCC*FLIQD*DTIMET
          TSOCGNN=ANUCB*ANUCC*GUG**(ANUCC-1.0D0)*FLIQD*
     .            DTIMET                ! temperature derivative of -DNU
         ENDIF
         RNU(JO)=RNODO
         ANU(JO)=RSHEL
         IF(INUCON.EQ.1) THEN           ! control minimum N
          IF(DNU(JO,1).LT.ANUCON) THEN
           DNU(JO,1)=0.0D0
           RNU(JO)=0.0D0                
           ANU(JO)=0.0D0
           JO=JO-1
          ENDIF
         ENDIF
        ENDIF
       ENDIF                            ! inucax.eq.2
C
       FSOLGX=0.0D0
       FSOLIX=0.0D0
       TSOCGNX=0.0D0
       TSOCGGX=0.0D0
       TSOCGIX=0.0D0
       AKA=(CAL-CAG)*DIFCA/(CLA-CAL)    ! for austenite shell growth
       AKG=(CAL-CAG)*DIFCA*3.646D0/
     .     (100.0D0-CAG)                ! for graphite nodule growth
       IF(JO.GT.0) THEN
        DO IJO=1,JO                     ! growth
         RRAG=AUSGR                                   ! constant
         IF(RRAG.EQ.0.0D0) RRAG=ANU(IJO)/RNU(IJO)     ! variable
         IF(IGROMX.EQ.1) THEN           ! Su model
          RNU(IJO)=RNU(IJO)+AKG*DTIMET/(RNU(IJO)*(1.0D0-1.0D0/RRAG))
          FIMPING=1.0D0
          IF(FLIQD.LT.0.5D0) FIMPING=(FLIQD/0.5D0)**(2.0D0/3.0D0)
          ANU(IJO)=ANU(IJO)+AKA*DTIMET/(ANU(IJO)*(RRAG-1.0D0))*FIMPING
          TSOCGGG=0.0D0                 ! temperature derivative of -RNU
          TSOCGIG=0.0D0                 ! temperature derivative of -ANU
         ENDIF
         IF(IGROMX.EQ.2) THEN           ! Su model (without impignment)
          RNU(IJO)=RNU(IJO)+AKG*DTIMET/(RNU(IJO)*(1.0D0-1.0D0/RRAG))
          FIMPING=1.0D0
          ANU(IJO)=ANU(IJO)+AKA*DTIMET/(ANU(IJO)*(RRAG-1.0D0))*FIMPING
          TSOCGGG=0.0D0                 ! temperature derivative of -RNU
          TSOCGIG=0.0D0                 ! temperature derivative of -ANU
         ENDIF
         IF(IGROMX.EQ.3) THEN           ! Lacaze model
          FIMPING=FLIQD                 ! simplification
          IF(FLIQD.LT.1.0D-2) FIMPING=1.0D-2
          RNU(IJO)=RNU(IJO)+AKG*DTIMET/(RNU(IJO)*(1.0D0-1.0D0/RRAG))*
     .             FIMPING
          AKL=1.0D0+(CLA-CAL)*2.646D0/
     .              (100.0D0-CAG)       ! for graphite nodule growth
          ANU(IJO)=ANU(IJO)+AKA*DTIMET/(ANU(IJO)*(RRAG-1.0D0))*FIMPING*
     .                      AKL
          TSOCGGG=0.0D0                 ! temperature derivative of -RNU
          TSOCGIG=0.0D0                 ! temperature derivative of -ANU
         ENDIF
         FSOLGX=FSOLGX+2.0D0/3.0D0*TWOPIT*DNU(IJO,1)*RNU(IJO)**3.0D0
         FSOLIX=FSOLIX+2.0D0/3.0D0*TWOPIT*DNU(IJO,1)*ANU(IJO)**3.0D0
C
         IF(IJO.EQ.JO)
     .    TSOCGNX=TSOCGNN*2.0D0/3.0D0*TWOPIT*ANU(IJO)**3.0D0
         TSOCGGX=TSOCGGX+TSOCGGG*2.0D0*TWOPIT*DNU(IJO,1)*RNU(IJO)**2.0D0
         TSOCGIX=TSOCGIX+TSOCGIG*2.0D0*TWOPIT*DNU(IJO,1)*ANU(IJO)**2.0D0
        ENDDO                           ! ijo=1,jo
       ENDIF                            ! jo gt 0
       IF(FSOLGX.GT.0.0D0) THEN
        FSOLG=FSOLGX
        FSOLA=FSOLIX-FSOLG
        TSOCGN=TSOCGNX
        TSOCGA=TSOCGIX
        IF(FSOLA.LT.FSOLAO) THEN        ! control
         FSOLA=FSOLAO
         FSOLG=FSOLIX-FSOLA
         TSOCGA=TSOCGGX                 ! assumption
        ENDIF
       ENDIF
C
       FSOLIX=FSOLA+FSOLG               ! control
       IF(FSOLIX.GT.1.0D0) THEN         ! assumption valid for small dt
        REAGX=1.0D0/FSOLIX              ! excess factor
        FSOLG=FSOLG*REAGX               ! (no control on DNU, RNU & ANU)
        FSOLA=FSOLA*REAGX
        TSOCGN=0.0D0                    ! assumption (it works!)
        TSOCGA=0.0D0
        IF(FSOLA.LT.FSOLAO) THEN        ! control
         FSOLA=FSOLAO
         FSOLG=1.0D0-FSOLA
        ENDIF
       ENDIF
      ENDIF                             ! gug.gt.0.0.and....
C
      FSOLID=FSOLA+FSOLG                ! always
C
C**** LIQUID FRACTION
C
      FLIQD=1.0D0-FSOLID                ! liquid fraction
      IF(FLIQD.LT.0.0D0) FLIQD=0.0D0    ! control
      IF(FLIQD.GT.1.0D0) FLIQD=1.0D0    ! control
C
C**** DEFINES THE SI CONTENT (MICROSEGREGATION)
C
      SICON=SIO*FLIQD**(AKSI-1.0D0)
C
C**** DEFINES MICROSTRUCTURAL-DEPENDENT "MACROSCOPICAL" PROPERTIES
C
      IF(IKMICX.EQ.1) THEN
       IF(FLIQD.LE.0.0D0)                     BASKK(1)=BASKS   ! Boeri's
       IF(FLIQD.GT.0.0D0.AND.FLIQD.LT.1.0D0)  BASKK(1)=BASKM
       IF(FLIQD.GE.1.0D0)                     BASKK(1)=BASKL
      ENDIF
      IF(IKMICX.EQ.2) THEN
       IF(FLIQDO.LE.0.0D0)                    BASKK(1)=BASKS   ! old fl
       IF(FLIQDO.GT.0.0D0.AND.FLIQD.LT.1.0D0) BASKK(1)=BASKM
       IF(FLIQDO.GE.1.0D0)                    BASKK(1)=BASKL
      ENDIF
      IF(IKMICX.EQ.3) THEN
       IF(FLIQD.LE.0.0D0)                     BASKK(1)=BASKS   ! smooth
       IF(FLIQD.GT.0.0D0.AND.FLIQD.LT.1.0D0)  BASKK(1)=BASKS+
     .                                               (BASKL-BASKS)*FLIQD
       IF(FLIQD.GE.1.0D0)                     BASKK(1)=BASKL
      ENDIF
C
C**** ESTABLISHES LATENT HEAT * PHASE-CHANGE FUNCTION (at time t & t+dt)
C
      TSOE1(IPLAT)=FLIQD                  ! f_pc at time t+dt
      TSOE2(IPLAT)=FLIQDO                 ! f_pc at time t
      TSOE2(IPLAT)=TSOE2(IPLAT)*HENER     ! L*f_pc at time t
      TSOE1(IPLAT)=TSOE1(IPLAT)*HENER     ! L*f_pc at time t+dt
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
   10 CONTINUE
C
C**** TRANSFER MICROSCOPICAL VARIABLES TO "ALPHAM"
C
      ALPHAM(IN+1)=FLIQD
C
      ALPHAM(IN+2)=FSOLA
      ALPHAM(IN+3)=FSOLG
C
      ALPHAM(IN+4)=SICON
C
      IF(JO.GT.0) THEN
       DO INU=1,JO
        ALPHAM(IN+5+INU-1)=DNU(INU,1)
        ALPHAM(IN+6+NNUM4T-1+INU-1)=RNU(INU)
        ALPHAM(IN+7+2*(NNUM4T-1)+INU-1)=ANU(INU)
       ENDDO
      ENDIF
C
      ALPHAM(IN+8+3*(NNUM4T-1))=FLOAT(JO)
      IF(INUCAX.EQ.1)
     . ALPHAM(IN+9+3*(NNUM4T-1))=FLOAT(INDEXG)
      IF(INUCAX.EQ.2)
     . ALPHAM(IN+9+3*(NNUM4T-1))=SINDEXG
C
C**** INCREMENTS "ALPHAM" INDEX
C
      INUPC=INUPC+4+3*NNUM4T+2         ! less than NBASES; see pointes.f
C
      RETURN
      END
