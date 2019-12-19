      SUBROUTINE MICROS4(TGAUAT,TGAUST,TGAUIT,DTEMPT,
     .                    BASMM, BASCC, BASKK,
     .                    TSOE2, TSOE1, TSOC1,
     .                    IPLAT,
     .                   ALPHAM,INUPC)  
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE PHASE-CHANGE FUNCTION ACCORDING TO THE
C     MICROSTRUCTURAL MODEL NUMBER 4 (IPCMO=4)
C
C     S.G. CAST IRON MICROSTRUCTURAL MODEL: BOERI'S MODEL
C                                                        (MULTI-NODULAR)
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
C      IMNMT=1 for manganese microsegregation (0 otherwise)
C     ALPHAM(IN+4+IMNMT)=manganese content (current)
c      I5A=5+IMNMT
C      NNUM4T=number of different grain densities and radii
C      I5=I5A+IMNMT+NNUM4T-1
C     ALPHAM(IN+I5A:IN+I5)=graphite grain density
C      I6A=I5+1
C      INNUM4T=1 for nodularity computation (0 otherwise)
C      I6B=I6A+(1+INNUM4T)*NNUM4T-1
C     ALPHAM(IN+I6A:IN+I6B)=graphite grain radius
C                       (two different radii for nodularity computation)
C      IWEUT=1 for white eutectic computation (0 otherwise)
C      I7=I6B+IWEUT
C     ALPHAM(IN+I7)=white eutectic fraction (solid)
C      I8=I7+IWEUT
C     ALPHAM(IN+I8)=white eutectic grain density
C      I9=I8+IWEUT
C     ALPHAM(IN+I9)=white eutectic radius
C
C     Auxiliar microstructural variables (not printed; see pointes.f):
C      I10=I6B+3*IWEUT+1
C     ALPHAM(IN+I10)=number of families
C      I11=I10+1
C     ALPHAM(IN+I11)=graphite nucleation index
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
      I2=2+INNUM4T
      FLIQD=ALPHAM(IN+1)
C
      IF(FLIQD.EQ.0.0D0) THEN          ! gets out if liquid fraction=0
       INUPC=INUPC+4+IMNMT+I2*NNUM4T+  ! less than NBASES; see pointes.f
     .       3*IWEUT+2
       RETURN
      ENDIF
C
      TWOPIT=6.283185307179586D0
C
C**** TRANSFERS "VPLAT" TO MICROSCOPICAL VARIABLES
C
      TEINF=VPLAT(IPLAT,1)         ! not used
      TESUP=VPLAT(IPLAT,2)         ! not used
      HENER=VPLAT(IPLAT,3)
      DELEE=TESUP-TEINF            ! not used
C
      CA =VPLAT(IPLAT, 6)
      SIO=VPLAT(IPLAT, 7)
      PH =VPLAT(IPLAT, 8)
      CU =VPLAT(IPLAT, 9)
      QMO=VPLAT(IPLAT,10)
      QG =VPLAT(IPLAT,11)
      CR =VPLAT(IPLAT,12)
C
      AKSI=VPLAT(IPLAT,13)
      AKQM=VPLAT(IPLAT,14)
      IHELAC=INT(VPLAT(IPLAT,15))
      IX0=0
      IF(IHELAC.EQ.4) THEN                    ! Thermocalc
       TEG=VPLAT(IPLAT,16)
       TEC=VPLAT(IPLAT,17)
       IX0=2
      ENDIF
C
      INUCMX=INT(VPLAT(IPLAT,16+IX0))         ! nucleation model index
      IF(INUCMX.EQ.1.OR.INUCMX.EQ.2) THEN     ! Su & Boeri nuc. models
       ANUCB=VPLAT(IPLAT,17+IX0)   ! nucleation coeff. A
       ANUCC=VPLAT(IPLAT,18+IX0)   ! nucleation coeff. n
      ENDIF
      INUCAX=                      ! nucleation arrest criterion index
     .       INT(VPLAT(IPLAT,19+IX0))
      INUCON=INT(VPLAT(IPLAT,20+IX0))         ! control minimum N index
      ANUCON=VPLAT(IPLAT,21+IX0)   ! minimum N
C
      IGROMX=INT(VPLAT(IPLAT,22+IX0))         ! growth model index
      IF(IGROMX.EQ.1) THEN                    ! Boeri growth model
       DIFCL=VPLAT(IPLAT,23+IX0)   ! diffusion coef. of C in liquid
       DIFCA=VPLAT(IPLAT,24+IX0)   ! diffusion coef. of C in austenite
       RNODA=VPLAT(IPLAT,25+IX0)   ! r of nodule enveloped by austenite
       RNODO=VPLAT(IPLAT,26+IX0)   ! initial radius of nodules
       AUSGR=VPLAT(IPLAT,27+IX0)   ! auste. shell radius/graphite radius
       IX1=0
      ENDIF
      IF(IGROMX.EQ.2) THEN                    ! Sandra growth model
       DIFCL=VPLAT(IPLAT,23+IX0)   ! diffusion coef. of C in liquid
       DIFCA=VPLAT(IPLAT,24+IX0)   ! diffusion coef. of C in austenite
       RNODO=VPLAT(IPLAT,25+IX0)   ! initial radius of nodules
       AUSGR=VPLAT(IPLAT,26+IX0)   ! auste. shell radius/graphite radius
       FLLOWN=VPLAT(IPLAT,27+IX0)  ! lower liquid fraction bound
       FLUPPN=VPLAT(IPLAT,28+IX0)  ! upper liquid fraction bound
       IX1=1
      ENDIF
C
      AWHIT=VPLAT(IPLAT,28+IX0+IX1)
      BWHIT=VPLAT(IPLAT,29+IX0+IX1)
C
      IKMICX=                      ! index for micro.-dep. conduct.
     .       INT(VPLAT(IPLAT,30+IX0+IX1))
      IX2=0
      IF(IKMICX.EQ.1.OR.IKMICX.EQ.2.OR.IKMICX.EQ.3) THEN
       IX2=3
       BASKS=VPLAT(IPLAT,31+IX0+IX1)            ! solid conductivity
       BASKM=VPLAT(IPLAT,32+IX0+IX1)            ! mushy conductivity
       BASKL=VPLAT(IPLAT,33+IX0+IX1)            ! liquid conductivity
      ENDIF
C
      IFPCDT=INT(VPLAT(IPLAT,31+IX0+IX1+IX2))   ! index for T derivative
      IAFLOJ=INT(VPLAT(IPLAT,32+IX0+IX1+IX2))   ! not used
C
      IV=32+IX0+IX1+IX2          ! number of VPLAT defined in input data
C
C**** TRANSFERS (REST OF) "ALPHAM" TO MICROSCOPICAL VARIABLES
C
      FSOLA=ALPHAM(IN+2)
      FSOLG=ALPHAM(IN+3)
C
      SICON=ALPHAM(IN+4)
C
      QMCON=QMO                  ! initial manganese content
      I5A=5
      IF(IMNMT.EQ.1) THEN        ! manganese microsegregation
       QMCON=ALPHAM(IN+5)
       I5A=6
      ENDIF
C
      I5=I5A+NNUM4T-1
      I6A=I5+1
      I6B=I6A+(1+INNUM4T)*NNUM4T-1
      I10=I6B+3*IWEUT+1
      I11=I10+1
C
      JO=INT(ALPHAM(IN+I10))
      IF(JO.GT.0) THEN
       DO INU=1,JO
        DNU(INU,1)=ALPHAM(IN+I5A+INU-1)
        RNU(INU)  =ALPHAM(IN+I6A+INU-1)
        IF(IGROMX.EQ.2) RNUZ1(INU)=ALPHAM(IN+I6A+NNUM4T+INU-1)
       ENDDO
      ENDIF
C
      FSOLC=0.0D0
      IF(IWEUT.EQ.1) THEN
       I7=I6B+1
       FSOLC=ALPHAM(IN+I7  )
       XXNCN=ALPHAM(IN+I7+1)
       RRRCM=ALPHAM(IN+I7+2)
      ENDIF
C
      IF(INUCAX.EQ.1)
     . INDEXG=INT(ALPHAM(IN+I11))              ! < NBASES, see pointes.f
      IF(INUCAX.EQ.2)
     . SINDEXG=ALPHAM(IN+I11)                  ! < NBASES, see pointes.f
C
C**** LOAD LAST CONVERGED VALUES
C
      FLIQDO=FLIQD
      FSOLAO=FSOLA
C
      IF(IWEUT.EQ.1) THEN
       XXNCNO=XXNCN
       RRRCMO=RRRCM
      ENDIF
C
C**** COMPUTES EQUILIBRIUM EUTECTIC TEMPERATURE & OTHER PARAMETERS
C
C     Nomenclature:
C     TGL: equilibrium graphite liquidus temperature
C     TEG: equilibrium graphite eutectic temperature
C     TEC: "equilibrium" cementite eutectic temperature
C     REAG: austenite volume/graphite volume (eutectic)
C     CLA: liquid %C at the liquid-austenite interface
C     CAL: austenite %C at the liquid-austenite interface
C     CLG: liquid %C at the liquid-graphite interface
C     CAG: austenite %C at the austenite-graphite interface
C
      IF(IHELAC.EQ.1) THEN                               ! Heine
       TL= 1569.00D0- 97.3D0*(CA+0.25D0*SICON)           ! hypoeutectic
c      TGL=389.10D0*(CA+0.31D0*SICON+0.33D0*PH-1.3D0)    ! hypereutectic
       TGL=389.10D0*(CA+SICON/3.0D0)-503.20D0            ! other option
c      TGL=389.10D0*(CA+0.31D0*SICON)-505.80D0           ! other option
       TEG=1154.6D0+6.5D0*SICON
       TEC=1148.0D0-15.0D0*SICON+3.0D0*QMCON-37.0D0*PH   ! Stefanescu
C
       CE=4.26D0-0.317D0*SICON
c      CEQ=CA+SICON/3.0D0                    ! equivalent carbon content
       CEQ=CA+SI0/3.0D0                      ! equivalent carbon content
       CEQREF=4.26D0                         ! ref. equiv. C content
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
       TEC=1148.0D0-15.0D0*SICON+3.0D0*QMCON-37.0D0*PH   ! Stefanescu
C
       CE=4.34D0-0.28D0*SICON
c      CEQ=CA+0.28D0*SICON                   ! equivalent carbon content
       CEQ=CA+0.28D0*SIO                     ! equivalent carbon content
       CEQREF=4.33D0                         ! ref. equiv. C content
C
       REAG=(26.26D0+0.087D0*SICON)/(2.16D0-0.101D0*SICON)
C
       CLA=(1576.0D0-TGAUST- 22.38D0*SICON)/ 97.3D0
       CAL=(1528.4D0-TGAUST- 32.02D0*SICON)/177.9D0
       CLG=( 534.5D0+TGAUST-112.84D0*SICON)/389.1D0
       CAG=(-1153.7D0+TGAUST-  4.8650D0*SICON)*(1.51D0-0.21D0*SICON)/
     .                       (353.7D0+4.865D0*SICON)+2.11D0-0.21D0*SICON
      ENDIF
      IF(IHELAC.EQ.3) THEN                               ! Escobar
       TL= 1576.3D0-97.3D0*CA-23.0D0*SICON-57.8D0*PH-5.66D0*QMCON-
     .                                                         2.71D0*CR
       TGL=389.15D0*CA+113.2D0*SICON+89.6D0*PH-2.4D0*QMCON-13.14D0*CR-
     .                                                           534.5D0
c      TCL=                             ! falta (cementita proeutectica)
       TEG=1154.0D0+4.24D0*SICON-28.61D0*PH-4.98D0*QMCON-4.75D0*CR
       TEC=1148.0D0-15.0D0*SICON-37.0D0*PH+3.0D0*QMCON+10.0D0*CR
C
       CE=4.33D0-0.28D0*SICON-0.30D0*PH-0.007D0*QMCON+0.021D0*CR
       CTE=2.08D0-0.11D0*SICON-0.35D0*PH+0.006D0*QMCON+0.07D0*CR
       CEQ=CA+0.28D0*SIO+                    ! equivalent carbon content
     .     0.30D0*PH-0.007D0*QMO+0.033D0*CR  ! (with initial Si and Mn)
       CEQREF=4.33D0                         ! ref. equiv. C content
C
       REAG=(100.0D0-CE)/3.646D0/(CE-CTE)
C
       CLA=(1576.3D0-TGAUST- 23.0D0*SICON-57.8D0*PH-5.66D0*QMCON-
     .                                                2.71D0*CR)/ 97.3D0
       CAL=(1528.4D0-TGAUST- 32.02D0*SICON)/177.9D0
       CLG=( 534.5D0+TGAUST-113.20D0*SICON-89.6D0*PH+2.4D0*QMCON+
     .                                               13.14D0*CR)/389.1D0
       CAG=(-1154.0D0+TGAUST-4.24D0*SICON+28.61D0*PH+4.98D0*QMCON+
     .      4.75D0*CR)*(1.48D0-0.11D0*SICON-0.35D0*PH+0.006D0*QMCON-
     .                                                       0.07D0*CR)/
     .         (354.0D0+4.24D0*SICON-28.61D0*PH-4.98D0*QMCON-4.75D0*CR)+
     .             2.08D0-0.11D0*SICON-0.35D0*PH+0.006D0*QMCON-0.07D0*CR
      ENDIF
      IF(IHELAC.EQ.4) THEN                               ! Thermocalc
       TL= 1576.3D0-97.3D0*CA-23.0D0*SICON-57.8D0*PH-5.66D0*QMCON-
     .                                                         2.71D0*CR
       TGL=389.15D0*CA+113.2D0*SICON+89.6D0*PH-2.4D0*QMCON-13.14D0*CR-
     .                                                           534.5D0
c      TCL=                             ! falta (cementita proeutectica)
       CE=4.33D0-0.28D0*SICON-0.30D0*PH-0.007D0*QMCON+0.021D0*CR
       CTE=2.08D0-0.11D0*SICON-0.35D0*PH+0.006D0*QMCON+0.07D0*CR
       CEQ=CA+0.28D0*SIO+                    ! equivalent carbon content
     .     0.30D0*PH-0.007D0*QMO+0.033D0*CR  ! (with initial Si and Mn)
       CEQREF=4.33D0                         ! ref. equiv. C content
C
       REAG=(100.0D0-CE)/3.646D0/(CE-CTE)
C
       CLA=(1576.3D0-TGAUST- 23.0D0*SICON-57.8D0*PH-5.66D0*QMCON-
     .                                                2.71D0*CR)/ 97.3D0
       CAL=(1528.4D0-TGAUST- 32.02D0*SICON)/177.9D0
       CLG=( 534.5D0+TGAUST-113.20D0*SICON-89.6D0*PH+2.4D0*QMCON+
     .                                               13.14D0*CR)/389.1D0
       CAG=(-1154.0D0+TGAUST-4.24D0*SICON+28.61D0*PH+4.98D0*QMCON+
     .      4.75D0*CR)*(1.48D0-0.11D0*SICON-0.35D0*PH+0.006D0*QMCON-
     .                                                       0.07D0*CR)/
     .         (354.0D0+4.24D0*SICON-28.61D0*PH-4.98D0*QMCON-4.75D0*CR)+
     .             2.08D0-0.11D0*SICON-0.35D0*PH+0.006D0*QMCON-0.07D0*CR
      ENDIF
C
C**** ESTABLISHES VARIABLES OF MODEL 4
C
      FSOLID=1.0D0-FLIQD        ! solid fraction
      FSOLIDO=1.0D0-FLIQDO
C
C**** SOLVES MODEL 4
C
      TSOCGN=0.0D0              ! only useful to thermal jacobian matrix
      TSOCGG=0.0D0
      TSOCGA=0.0D0
      TSOCCN=0.0D0
      TSOCCG=0.0D0
C
C**** GRAPHITE EUTECTIC SOLIDIFICATION
C
      GUG=TEG-TGAUST                    ! graphite eutectic undercooling
      ICORREC=0                         ! not used
C
      IF(GUG.GT.0.0D0.AND.FSOLG.LT.1.0D0.AND.FSOLID.LT.1.0D0) THEN
       TSOCGNN=0.0D0
       IF(INUCAX.EQ.1) THEN             ! nucleation up to Trecal
        IF(JO.GT.0.AND.DTEMPT.GT.0.0D0) INDEXG=1          ! end of nucl.
        IF(INDEXG.EQ.0) THEN            ! nucleation
         JO=JO+1
         IF(JO.GT.NNUM4T)               ! see pointes.f
     .    CALL RUNENDT('ERROR IN MICROS4: JO GT NNUM4T')
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
         RNU(JO)=RNODO
         IF(IGROMX.EQ.2) RNUZ1(JO)=RNODO
         IF(INUCON.EQ.1) THEN           ! control minimum N
          IF(DNU(JO,1).LT.ANUCON) THEN
           DNU(JO,1)=0.0D0
           RNU(JO)=0.0D0
           RNUZ1(JO)=0.0D0
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
     .    CALL RUNENDT('ERROR IN MICROS4: JO GT NNUM4T')
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
         RNU(JO)=RNODO
         IF(IGROMX.EQ.2) RNUZ1(JO)=RNODO
         IF(INUCON.EQ.1) THEN           ! control minimum N
          IF(DNU(JO,1).LT.ANUCON) THEN
           DNU(JO,1)=0.0D0
           RNU(JO)=0.0D0
           RNUZ1(JO)=0.0D0
           JO=JO-1
          ENDIF
         ENDIF
        ENDIF
       ENDIF                            ! inucax.eq.2
C
       FSOLGX=0.0D0
       TSOCGGX=0.0D0
       IF(IGROMX.EQ.1) THEN             ! Boeri growth
        IF(JO.GT.0) THEN
         DO IJO=1,JO
          IF(RNU(IJO).LT.RNODA) THEN
           AKK=(CLA-CLG)*DIFCL*3.646D0/ ! growth in contact with liquid
     .                                 (100.0D0-CLG)
           RNU(IJO)=RNU(IJO)+(-2.0D0*RNU(IJO)+((2.0D0*RNU(IJO))**2.0D0+
     .                       4.0D0*AKK*DTIMET)**0.5D0)/2.0D0
C
           TSOCGGG=0.0D0                ! temperature derivative of -RNU
c          CDERCRL=                     ! * (see below; to be checked !)
c    .           -DIFCL*3.646D0*DTIMET*(-0.00642376285/(98.70675919-
c    .            0.002570033410*TGAUST+0.333333333333*SICON)+
c    .            0.001285016705*(14.83214459-0.0128475257*TGAUST+
c    .            0.0833847208*SICON)/(98.70675919-0.002570033410*
c    .            TGAUST+0.333333333333*SICON)**2)
c          TSOCGGG=CDERCRL              ! temperature derivative of -RNU
          ELSE
           AKK=(CAL-CAG)*DIFCA*3.646D0/ ! growth in contact with auste.
     .                                 (100.0D0-CAG)
           RNU(IJO)=RNU(IJO)+0.9D0*AKK*FLIQD**0.66D0*DTIMET/
     .                       (RNU(IJO)*(1.0D0-1.0D0/AUSGR))
C
           TSOCGGG=0.0D0                ! temperature derivative of -RNU
c          CDERCRA=                     ! *
c    .           -1.911*DIFCA*3.646D0*DTIMET*FLIQD**(2/3)*
c    .           ((-0.005621135469-(1.5-0.216*SICON)/(354.6+6.5*SICON))/
c    .           (97.9-(TGAUST-1154.6-6.5*SICON)*(1.5-0.216*SICON)/
c    .           (354.6+6.5*SICON)+0.21*SICON)+(6.491343451-
c    .           0.005621135*TGAUST+0.03012665*SICON-(TGAUST-1154.6-
c    .           6.5*SICON)*(1.5-0.216*SICON)/(354.6+6.5*SICON))*
c    .           (1.5-0.216*SICON)/((97.9-((TGAUST-1154.6-6.5*SICON)*
c    .           (1.5-0.216*SICON)/(354.6+6.5*SICON)+0.21*SICON)**2*
c    .           (354.6+6.5*SICON)))
c          TSOCGGG=CDERCRA              ! temperature derivative of -RNU
          ENDIF
          FSOLGX=FSOLGX+2.0D0/3.0D0*TWOPIT*DNU(IJO,1)*RNU(IJO)**3.0D0
C
          IF(IJO.EQ.JO)
     .     TSOCGNX=TSOCGNN*2.0D0/3.0D0*TWOPIT*RNU(IJO)**3.0D0
          TSOCGGX=TSOCGGX+
     .                   TSOCGGG*2.0D0*TWOPIT*DNU(IJO,1)*RNU(IJO)**2.0D0
         ENDDO                          ! ijo=1,jo
        ENDIF                           ! jo gt 0
       ENDIF                            ! igromx.eq.1
       IF(IGROMX.EQ.2) THEN             ! Sandra growth
        IF(JO.GT.0) THEN
         DO IJO=1,JO
          IF(FLIQD.GT.FLUPPN) THEN
           AKK=(CLA-CLG)*DIFCL*3.646D0/ ! growth in contact with liquid
     .                                 (100.0D0-CLG)
           DRNUL=(-2.0D0*RNU(IJO)+((2.0D0*RNU(IJO))**2.0D0+
     .                       4.0D0*AKK*DTIMET)**0.5D0)/2.0D0
           RNU(IJO)=RNU(IJO)+DRNUL
           RNUZ1(IJO)=RNUZ1(IJO)+DRNUL
C
           TSOCGGG=0.0D0                ! temperature derivative of -RNU
          ELSEIF(FLIQD.LE.FLUPPN.AND.FLIQD.GT.FLLOWN) THEN
           AKK=(CLA-CLG)*DIFCL*3.646D0/ ! growth in contact with liquid
     .                                 (100.0D0-CLG)
           DRNUL=(-2.0D0*RNU(IJO)+((2.0D0*RNU(IJO))**2.0D0+
     .                       4.0D0*AKK*DTIMET)**0.5D0)/2.0D0
           RNU(IJO)=RNU(IJO)+DRNUL
           AKK=(CAL-CAG)*DIFCA*3.646D0/ ! growth in contact with auste.
     .                                 (100.0D0-CAG)
           DRNUA=+0.9D0*AKK*FLIQD**0.66D0*DTIMET/
     .                       (RNUZ1(IJO)*(1.0D0-1.0D0/AUSGR))
           RNUZ1(IJO)=RNUZ1(IJO)+DRNUA
C
           TSOCGGG=0.0D0                ! temperature derivative of -RNU
          ELSE
           AKK=(CAL-CAG)*DIFCA*3.646D0/ ! growth in contact with auste.
     .                                 (100.0D0-CAG)
           DRNUA=+0.9D0*AKK*FLIQD**0.66D0*DTIMET/
     .                       (RNUZ1(IJO)*(1.0D0-1.0D0/AUSGR))
           RNU(IJO)=RNU(IJO)+DRNUA
           RNUZ1(IJO)=RNUZ1(IJO)+DRNUA
C
           TSOCGGG=0.0D0                ! temperature derivative of -RNU
          ENDIF
          FSOLGX=FSOLGX+1.0D0/3.0D0*TWOPIT*DNU(IJO,1)*
     .                               (RNU(IJO)**3.0D0+RNUZ1(IJO)**3.0D0)
C
          IF(IJO.EQ.JO)
     .     TSOCGNX=TSOCGNN*1.0D0/3.0D0*TWOPIT*
     .                               (RNU(IJO)**3.0D0+RNUZ1(IJO)**3.0D0)
          TSOCGGX=TSOCGGX+
     .                   TSOCGGG*1.0D0*TWOPIT*DNU(IJO,1)*
     .                               (RNU(IJO)**2.0D0+RNUZ1(IJO)**3.0D0)
         ENDDO                          ! ijo=1,jo
        ENDIF                           ! jo gt 0
       ENDIF                            ! igromx.eq.2
C
       IF(FSOLGX.GT.0.0D0) THEN
        FSOLG=FSOLGX
        FSOLA=REAG*FSOLG
        IF(FSOLA.LT.FSOLAO) THEN        ! control
         FSOLA=FSOLAO
         FSOLG=FSOLA/REAG
        ENDIF
        TSOCGN=TSOCGNX
        TSOCGG=TSOCGGX
        TSOCGA=REAG*(TSOCGN+TSOCGG)
       ENDIF
C
       FSOLIX=FSOLA+FSOLG+FSOLC         ! control
       IF(FSOLIX.GT.1.0D0) THEN         ! assumption valid for small dt
        FSOLG=1.0D0/(1.0D0+REAG)        ! (no control on DNU & RNU)
        FSOLA=REAG*FSOLG
        TSOCGN=0.0D0                    ! assumption (it works!)
        TSOCGA=0.0D0
        IF(FSOLA.LT.FSOLAO) THEN        ! control
         FSOLA=FSOLAO
         FSOLG=1.0D0-FSOLA
        ENDIF
       ENDIF
       FSOLID=FSOLA+FSOLG+FSOLC         ! always
      ENDIF                             ! gug.gt.0.0.and....
C
C**** WHITE EUTECTIC SOLIDIFICATION
C
C     Assumption: the white eutectic nucleation arrest model is the same
C                 as that of the graphite eutectic
C
      IF(IWEUT.EQ.1) THEN
       GUC=TEC-TGAUST                   ! white eutectic undercooling
       ICORREC=0
C
       IF(GUC.GT.0.0D0.AND.FSOLC.LT.1.0D0.AND.FSOLID.LT.1.0D0) THEN
        INDEXC=1
        IF(INUCAX.EQ.1) THEN
         IF(INDEXG.EQ.0) INDEXC=0
        ENDIF
        IF(INUCAX.EQ.2) THEN
         IF(DTEMPT.LT.0.0D0.AND.TGAUST.LT.SINDEXG) INDEXC=0
        ENDIF
        IF(INDEXC.EQ.0) THEN
         XXNCN=AWHIT*GUC*GUC            ! nucleation
         TSOCCN=AWHIT*2.0D0*GUC
        ENDIF
C
        IF(XXNCN.LT.XXNCNO) THEN        ! control
         XXNCN=XXNCNO
         TSOCCN=0.0D0
        ENDIF
C
        FSOLC=2.0D0/3.0D0*TWOPIT*XXNCN*RRRCM**3.0D0            ! control
        FSOLIX=FSOLA+FSOLG+FSOLC
        IF(FSOLIX.GT.1.0D0) THEN        ! assumption valid for small dt
         FSOLC=1.0D0-FSOLA-FSOLG
         F1=XXNCN
         XXNCN=3.0D0*FSOLC/(2.0D0*TWOPIT*RRRCM**3.0D0)
         TSOCCN=TSOCCN*XXNCN/F1
         IF(XXNCN.LT.XXNCNO) THEN
          XXNCN=XXNCNO
          TSOCCN=0.0D0
         ENDIF
        ENDIF
C
        IF(XXNCN.GT.0.0D0) THEN
         RRRCM=RRRCM+                   ! growth
     .               BWHIT*GUC*GUC*DTIMET
         TSOCCG=BWHIT*2.0D0*GUC*DTIMET
C
         IF(RRRCM.LT.RRRCMO) THEN       ! control
          RRRCM=RRRCMO
          TSOCCG=0.0D0
         ENDIF
C
         FSOLC=2.0D0/3.0D0*TWOPIT*XXNCN*RRRCM**3.0D0           ! control
         FSOLIX=FSOLA+FSOLG+FSOLC
         IF(FSOLIX.GT.1.0D0) THEN       ! assumption valid for small dt
          FSOLC=1.0D0-FSOLA-FSOLG
          F1=RRRCM
          RRRCM=(3.0D0*FSOLC/(2.0D0*TWOPIT*XXNCN))**(1.0D0/3.0D0)
          TSOCCG=TSOCCG*RRRCM/F1
          IF(RRRCM.LT.RRRCMO) THEN
           RRRCM=RRRCMO
           TSOCCG=0.0D0
          ENDIF
         ENDIF
C
         TSOCCN=TSOCCN*2.0D0/3.0D0*TWOPIT*RRRCM**3.0D0
         TSOCCG=TSOCCG*2.0D0*TWOPIT*XXNCN*RRRCM**2.0D0
        ENDIF                           ! xxncn.gt.0.0
        FSOLID=FSOLA+FSOLG+FSOLC        ! always
       ENDIF                            ! guc.gt.0.0.and....
      ENDIF                             ! iweut=1
C
C**** LIQUID FRACTION
C
      FLIQD=1.0D0-FSOLID                   ! liquid fraction
      IF(FLIQD.LT.0.0D0) FLIQD=0.0D0       ! control
      IF(FLIQD.GT.1.0D0) FLIQD=1.0D0       ! control
C
C**** DEFINES THE SI & MN CONTENTS (MICROSEGREGATION)
C
      SICON=SIO*FLIQD**(AKSI-1.0D0)
      IF(IMNMT.EQ.1) QMCON=QMO*FLIQD**(1.0D0-AKQM)
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
      TSOC1(IPLAT)=TSOCGN+TSOCGG+TSOCGA+TSOCCN+TSOCCG
      TSOC1(IPLAT)=TSOC1(IPLAT)*HENER
      GO TO 10
C
    4 ICACOT=1
      TSOC1(IPLAT)=TSOCGN+TSOCGG+TSOCGA+TSOCCN+TSOCCG
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
      IF(IMNMT.EQ.1) ALPHAM(IN+5)=QMCON
C
      IF(JO.GT.0) THEN
       DO INU=1,JO
        ALPHAM(IN+I5A+INU-1)=DNU(INU,1)
        ALPHAM(IN+I6A+INU-1)=RNU(INU)
        IF(IGROMX.EQ.2) ALPHAM(IN+I6A+NNUM4T+INU-1)=RNUZ1(INU)
       ENDDO
      ENDIF
C
      IF(IWEUT.EQ.1) THEN
       ALPHAM(IN+I7  )=FSOLC
       ALPHAM(IN+I7+1)=XXNCN
       ALPHAM(IN+I7+2)=RRRCM
      ENDIF
C
      ALPHAM(IN+I10)=FLOAT(JO)
      IF(INUCAX.EQ.1)
     . ALPHAM(IN+I11)=FLOAT(INDEXG)
      IF(INUCAX.EQ.2)
     . ALPHAM(IN+I11)=SINDEXG
C
C**** INCREMENTS "ALPHAM" INDEX
C
      INUPC=INUPC+4+IMNMT+I2*NNUM4T+   ! less than NBASES; see pointes.f
     .      3*IWEUT+2
C
      RETURN
      END
