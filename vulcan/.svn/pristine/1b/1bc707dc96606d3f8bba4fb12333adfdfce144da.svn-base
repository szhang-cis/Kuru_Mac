      SUBROUTINE MICROS3(TGAUAT,TGAUST,TGAUIT,DTEMPT,
     .                    BASMM, BASCC, BASKK,
     .                    TSOE2, TSOE1, TSOC1,
     .                    IPLAT,
     .                   ALPHAM,INUPC)  
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE PHASE-CHANGE FUNCTION, DENSITY,
C     CAPACITY & CONDUCTIVITY ACCORDING WITH THE MICROSTRUCTURAL
C     MODEL NUMBER 3 (IPCMO=3)
C
C     ALUMINIUM ALLOY MICROSTRUCTURAL MODEL:
C
C     Nucleation and growth models for dendritic & eutectic
C     solidification
C
C***********************************************************************
C
C     Index of variables
C
C     Input:
C     TGAUAT= Temperature at time t
C     TGAUST= Temperature at time t+dt
C     DTEMPT= Temperature rate
C
C     Input/output:
C     BASMM = Density
C     BASCC = Capacity coefficient
C     BASKK = Conductivity
C
C     Output:
C     TSOE2 = each L*phase-change function at time t
C     TSOE1 = each L*phase-change function at time t+dt
C
C     ALPHAM= array of microstructural (microscopical) variables
C
C     ALPHAM(IN+1)=liquid fraction
C
C     ALPHAM(IN+2)=dendritic solid fraction (dendritic solidification)
C     ALPHAM(IN+3)=dendritic grain solid fraction (dendritic solid.)
C     ALPHAM(IN+4)=dendritic internal solid fraction (dendritic solid.)
C     ALPHAM(IN+5)=dendritic grain density (dendritic solidification)
C     ALPHAM(IN+6)=dendritic grain radius (dendritic solidification)
C
C     ALPHAM(IN+7)= eutectic solid fraction (eutectic solidification)
C     ALPHAM(IN+8)= eutectic grain solid fraction (eutectic solid.)
C     ALPHAM(IN+9)= interdendritic eutectic solid fraction (eut. solid.)
C     ALPHAM(IN+10)=eutectic internal solid fraction (eutectic solid.)
C     ALPHAM(IN+11)=eutectic grain density
C     ALPHAM(IN+12)=eutectic grain radius
C
C     Auxiliar microstructural variables (not printed; see pointes.f):
C     ALPHAM(IN+13)=index for dendritic nucleation
C     ALPHAM(IN+14)=index for eutectic nucleation
C
C
C     Dendritic solidification
C
C     Nucleation models: 1=instantaneous (Gaussian distribution)
C                        2=idem model 1 with end of nucleation at
C                          maximum undercooling with zero temperature
C                          rate
C
C     Growth models:     1=radius rate in terms of undercooling &
C                          internal fraction equals to one
C                        2=radius rate in terms of undercooling &
C                          internal rate given by linear law
C                        3=radius rate in terms of undercooling &
C                          internal fraction given by Scheil's equation
C                        4=radius rate in terms of undercooling &
C                          internal rate given by lever rule
C                        5=radius rate in terms of undercooling &
C                          ???? (to be implemented!, Thevoz model)
C
C
C     Eutectic solidification
C
C     Nucleation models: 1=instantaneous (Gaussian distribution)
C                        2=idem model 1 with end of nucleation at
C                          maximum undercooling with zero temperature
C                          rate
C
C     Growth models:     1=zero (eutectic) grain solid fraction &
C                          internal fraction equals to zero
C                        2=zero (eutectic) grain solid fraction &
C                          internal fraction given by linear law
C                        3=zero (eutectic) grain solid fraction &
C                          internal fraction given by Scheil's equation
C                        4=zero (eutectic) grain solid fraction &
C                          internal fraction given by lever rule
C                        5=radius rate in terms of undercooling &
C                          internal fraction equals to zero
C                        6=radius rate in terms of undercooling &
C                          grain and interdendritic fractions derived
C                          in a proportional form respectively
C
C
C     Notes:
C
C     IAFLOJ=fraction correction form (=0, relaxed correction;
C                                      =1, full correction)
C
C     Improvements to be done:
C     A better input of the non-equilibrium partition coefficient
C     (January/2001)
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
      TEINF=VPLAT(IPLAT,1)
      TESUP=VPLAT(IPLAT,2)
      HENER=VPLAT(IPLAT,3)
      DELEE=TESUP-TEINF
C
      SICON=VPLAT(IPLAT,6)
      SBCON=VPLAT(IPLAT,7)
      SICONE=VPLAT(IPLAT,8)
      NLICO=INT(VPLAT(IPLAT,9))
      PARCO=VPLAT(IPLAT,10)
C
      IMNUC=INT(VPLAT(IPLAT,11))
      IMGRO=INT(VPLAT(IPLAT,12))
      IMNUCE=INT(VPLAT(IPLAT,13))
      IMGROE=INT(VPLAT(IPLAT,14))
C
      IF(IMNUC.EQ.1.OR.IMNUC.EQ.2) THEN
       TEMAV=VPLAT(IPLAT,14+1)
       TEMDE=VPLAT(IPLAT,14+2)
       GRDEM=VPLAT(IPLAT,14+3)
       IVN=3
      ENDIF
C
      IF(IMGRO.EQ.1.OR.IMGRO.EQ.2.OR.IMGRO.EQ.3.OR.IMGRO.EQ.4.OR.
     .   IMGRO.EQ.5) THEN
       NRAUN=INT(VPLAT(IPLAT,14+IVN+1))
       DIFLI=VPLAT(IPLAT,14+IVN+2)   ! non-equilibrium partition coef.
       HENER=VPLAT(IPLAT,14+IVN+3)   ! redefines latent heat (dendritic)
       IVG=3
      ENDIF
C
      IF(IMNUCE.EQ.1.OR.IMNUCE.EQ.2) THEN
       TEMAVE=VPLAT(IPLAT,14+IVN+IVG+1)
       TEMDEE=VPLAT(IPLAT,14+IVN+IVG+2)
       GRDEME=VPLAT(IPLAT,14+IVN+IVG+3)
       IVNE=3
      ENDIF
C
      IF(IMGROE.EQ.1.OR.IMGROE.EQ.2.OR.IMGROE.EQ.3.OR.
     .   IMGROE.EQ.4.OR.IMGROE.EQ.5.OR.IMGROE.EQ.6) THEN
       NRAUNE=INT(VPLAT(IPLAT,14+IVN+IVG+IVNE+1))
       HENERE=VPLAT(IPLAT,14+IVN+IVG+IVNE+2) ! redef. lat. heat (eutec.)
       IVGE=2
      ENDIF
C
      IFPCDT= INT(VPLAT(IPLAT,14+IVN+IVG+IVNE+IVGE+1))
      IAFLOJ= INT(VPLAT(IPLAT,14+IVN+IVG+IVNE+IVGE+2))
      NSUGUEM=INT(VPLAT(IPLAT,14+IVN+IVG+IVNE+IVGE+3))
      EXTGUE=     VPLAT(IPLAT,14+IVN+IVG+IVNE+IVGE+4)
      IVFPC=4
C
      IV=                        ! number of VPLAT defined in input data
     .   14+IVN+IVG+IVNE+IVGE+IVFPC
C
C**** TRANSFER "ALPHAM" TO MICROSCOPICAL VARIABLES
C
      IN=INUPC
      FLIQD=ALPHAM(IN+1)
C
      FSOLD=ALPHAM(IN+2)
      FSOLG=ALPHAM(IN+3)
      FINTD=ALPHAM(IN+4)
      DENNU=ALPHAM(IN+5)
      RADGR=ALPHAM(IN+6)
C
      FSOLDE= ALPHAM(IN+7)
      FSOLGE= ALPHAM(IN+8)
      FSOLGDE=ALPHAM(IN+9)
      FINTDE= ALPHAM(IN+10)
      DENNUE= ALPHAM(IN+11)
      RADGRE= ALPHAM(IN+12)
C
      INDEXD= INT(ALPHAM(IN+13))
      INDEXE= INT(ALPHAM(IN+14))           ! less than 19; see pointes.f
C
      IF(IMNUC.EQ.1)  INDEXD= 0
      IF(IMNUCE.EQ.1) INDEXE= 0
C
C**** LOAD LAST CONVERGED VALUES
C
      FLIQDO=FLIQD
C
      FSOLDO=FSOLD
      FSOLGO=FSOLG
      FINTDO=FINTD
      DENNUO=DENNU
      RADGRO=RADGR
C
      FSOLDEO= FSOLDE
      FSOLGEO= FSOLGE
      FSOLGDEO=FSOLGDE
      FINTDEO= FINTDE
      DENNUEO= DENNUE
      RADGREO= RADGRE
C
      INDEXDO=INDEXD
      INDEXEO=INDEXE
C
C**** ESTABLISH VARIABLES OF MODEL 3
C
      FSOLID=1.0D0-FLIQD        ! solid fraction
      FSOLIDO=1.0D0-FLIQDO
C
      FINTF=FINTD+FINTDE        ! total internal fraction
C
      IF(SICON.LT.(SICONE*PARCO)) THEN
       INDEX1=1                 ! dendritic solidification
       INDEX2=0                 ! no eutectic solidification
      ENDIF
      IF(SICON.GE.(SICONE*PARCO).AND.SICON.LT.SICONE) THEN
       INDEX1=1                 ! dendritic solidification
       INDEX2=1                 ! eutectic solidification
      ENDIF
      IF(SICON.EQ.SICONE) THEN
       INDEX1=0                 ! no dendritic solidification
       INDEX2=1                 ! eutectic solidification
      ENDIF
      IF(SICON.GT.SICONE) THEN
       CALL RUNENDT('ERROR: HYPEREUTECTIC NOT IMPLEMENTED - MODEL 3')
      ENDIF
C
      IF(INDEX1.EQ.1) THEN              ! dendritic solidification
       IF(IMGRO.EQ.1.OR.IMGRO.EQ.2.OR.IMGRO.EQ.3.OR.IMGRO.EQ.4.OR.
     .    IMGRO.EQ.5) THEN
        CALL LIQCOMT(SICON,PROPST,TESUP,
     .               NLICO,    1)       ! redefines liquidus temperature
c       ISBCON=1                        ! better as input
c       IF(ISBCON.EQ.1) TESUP=TESUP-300.0D0*SBCON
C
        SICOX=0.0D0
        CALL LIQCOMT(SICOX,PROPST,TESUPF,
     .               NLICO,    1)       ! obtains melting T pure metal
       ENDIF
       IF(IMGRO.EQ.2.OR.IMGRO.EQ.3.OR.IMGRO.EQ.4) THEN
c       CALL LIQCOMT(SICON,PROPST,TEINF,    !!!!!!! to be improved !!!!
c    .               NLICO,    8)       ! redefines solidus temperature
        DELEE=TESUP-TEINF
       ENDIF
      ENDIF
      IF(INDEX2.EQ.1) THEN              ! eutectic solidification
       IF(IMGROE.EQ.1.OR.IMGROE.EQ.2.OR.IMGROE.EQ.3.OR.
     .    IMGROE.EQ.4.OR.IMGROE.EQ.5.OR.IMGROE.EQ.6) THEN
        CALL LIQCOMT(SICONE,PROPST,TESUPE,
     .               NLICO,    1)       ! obtains eutectic temperature
c       ISBCON=1                        ! better as input
c       IF(ISBCON.EQ.1) TESUPE=TESUPE-300.0D0*SBCON
       ENDIF
       IF(IMGROE.EQ.2.OR.IMGROE.EQ.3.OR.IMGROE.EQ.4) THEN
c       CALL LIQCOMT(SICON,PROPST,TEINF,    !!!!!!! to be improved !!!!
c    .               NLICO,    8)       ! redefines solidus temperature
        DELEE=TESUP-TEINF
       ENDIF
      ENDIF
C
C**** SOLVES MODEL 3
C
      DELNU=0.0D0
      DELRA=0.0D0
C
      DELNUE=0.0D0
      DELRAE=0.0D0
C
      TSOCN=0.0D0
      TSOCR=0.0D0
      TSOCF=0.0D0
C
      TSOCNE=0.0D0
      TSOCRE=0.0D0
      TSOCFE=0.0D0
C
C**** DENDRITIC SOLIDIFICATION
C
      IF(INDEX1.EQ.1) THEN
       GU=TESUP-TGAUST            ! dendritic undercooling
       ICORREC=0
C
C**** NUCLEATION MODEL
C
       IF(DENNU.GT.0.0D0.AND.DTEMPT.GT.0.0D0) INDEXD=1    ! end of nucl.
       IF(IMNUC.EQ.1.OR.IMNUC.EQ.2) THEN           ! instantaneous nucl.
        IF(GU.GT.0.0D0.AND.DTEMPT.LT.0.0D0.AND.FSOLG.LT.1.0D0
     .                .AND.FSOLD.LT.1.0D0.AND.FSOLID.LT.1.0D0
     .                .AND.INDEXD.EQ.0) THEN
C
         NSUGUE=1                                    ! subincrementation
         IF(EXTGUE.GT.0.0D0) THEN
          IF(TEMAV.GT.0.0D0) THEN
           NSUGUE=INT(GU/(EXTGUE*TEMAV))
          ELSE
           NSUGUE=INT(GU/EXTGUE)
          ENDIF
          IF(NSUGUE.LT.1)       NSUGUE=1
          IF(NSUGUE.GT.NSUGUEM) NSUGUE=NSUGUEM
         ENDIF
         IF(NSUGUE.EQ.1) THEN
          AAAAA=((GU-TEMAV)/(2.0D0*TEMDE))**2.0D0
          FNUCL=1.0D0/(DSQRT(TWOPIT))*DEXP(-AAAAA)
          DELNU=-GRDEM*FNUCL*DTEMPT
          DENNU=DENNU+DELNU*DTIMET
          IF(DENNU.GT.GRDEM) DENNU=GRDEM
         ELSE
          DTEMPX=DTEMPT/NSUGUE
          TGAUSX=TGAUAT
          DO I=1,NSUGUE
           TGAUSX=TGAUSX+DTEMPX*DTIMET
           GUX=TESUP-TGAUSX
           IF(GUX.GT.0.0D0) THEN
            AAAAA=((GUX-TEMAV)/(2.0D0*TEMDE))**2.0D0
            FNUCL=1.0D0/(DSQRT(TWOPIT))*DEXP(-AAAAA)
            DELNU=-GRDEM*FNUCL*DTEMPX
            DENNU=DENNU+DELNU*DTIMET
            IF(DENNU.GT.GRDEM) DENNU=GRDEM
           ENDIF
          ENDDO
         ENDIF
C
         FSOLG=2.0D0/3.0D0*TWOPIT*DENNU*RADGR**3.0D0           ! control
         FSOLIX=FSOLG+FSOLGE
         IF(FSOLIX.GT.1.0D0) THEN        ! assumption valid for small dt
          FSOLG=1.0D0-FSOLGE
          DENNU=3.0D0*FSOLG/(2.0D0*TWOPIT*RADGR**3.0D0)
          IF(DENNU.LT.DENNUO) DENNU=DENNUO
          FNUCL=-(DENNU-DENNUO)/(DTIMET*GRDEM*DTEMPT)
         ENDIF
C
         TSOCN=GRDEM*FNUCL
         ICORREC=1
        ENDIF                            ! gu.gt.0.and....
       ENDIF                             ! imnuc.eq.1.or....
C
       IF(ICORREC.EQ.1) THEN             ! correction
        IF(IAFLOJ.EQ.1) THEN             ! strict conditions
         FSOLD=FSOLG*FINTD
         FSOLGDE=FSOLG*FINTDE
         FSOLDE=FSOLGDE+FSOLGE
         FSOLID=FSOLD+FSOLDE
        ENDIF                            ! iafloj.eq.1
       ENDIF                             ! icorrec.eq.1
C
C**** GROWTH MODEL
C
C     a) radius growth
C     b) internal fraction growth
C
       ICORREC=0                                                    ! a)
       IF(IMGRO.EQ.1.OR.IMGRO.EQ.2.OR.IMGRO.EQ.3.OR.IMGRO.EQ.4.OR.
     .    IMGRO.EQ.5) THEN
        IF(GU.GT.0.0D0.AND.FSOLG.LT.1.0D0.AND.FSOLD.LT.1.0D0
     .                .AND.FSOLID.LT.1.0D0.AND.DENNU.GT.0.0D0) THEN
         CALL LIQCOMT(DELRA,PROPST,GU,   ! computes grain velocity (v)
     .                NRAUN,    4)
         RADGR=RADGR+DELRA*DTIMET
C
         FSOLG=2.0D0/3.0D0*TWOPIT*DENNU*RADGR**3.0D0           ! control
         FSOLIX=FSOLG+FSOLGE
         IF(FSOLIX.GT.1.0D0) THEN        ! assumption valid for small dt
          FSOLG=1.0D0-FSOLGE
          RADGR=(3.0D0*FSOLG/(2.0D0*TWOPIT*DENNU))**(1.0D0/3.0D0)
          IF(RADGR.LT.RADGRO) RADGR=RADGRO
C                                        ! no control on TSOCR
         ENDIF
C
         CALL LIQCOMT(DEDELRA,PROPST,GU, ! computes dv/dT
     .                NRAUN,    5)
         TSOCR=DEDELRA*DTIMET
         ICORREC=1
        ENDIF                            ! gu.gt.0.0.and....
       ENDIF                             ! imgro.eq.1.or....
C
       IF(ICORREC.EQ.1) THEN             ! correction
        IF(IAFLOJ.EQ.1) THEN             ! strict conditions
         FSOLD=FSOLG*FINTD
         FSOLGDE=FSOLG*FINTDE
         FSOLDE=FSOLGDE+FSOLGE
         FSOLID=FSOLD+FSOLDE
        ENDIF                            ! iafloj.eq.1
       ENDIF                             ! icorrec.eq.1
C
       ICORREC=0                                                    ! b)
       IF(IMGRO.EQ.1) THEN               ! internal fracion=1.0
        IF(GU.GT.0.0D0.AND.FINTD.LT.1.0D0.AND.FINTF.LT.1.0D0) THEN
         IF(FSOLG.GT.0.0D0) THEN         ! dendritic internal fraction
          FINTD=1.0D0
          FINTX=FINTD+FINTDE                                   ! control
          IF(FINTX.GT.1.0D0) THEN        ! assumption valid for small dt
           FINTD=1.0D0-FINTDE
C                                        ! no control on TSOCF
          ENDIF
          TSOCF=0.0D0
          ICORREC=1
         ENDIF                           ! fsolg.lt.1.0
        ENDIF                            ! gu.gt.0.and....
       ENDIF                             ! imgro.eq.1
       IF(IMGRO.EQ.2) THEN               ! linear model
        IF(GU.GT.0.0D0.AND.FINTD.LT.1.0D0.AND.FINTF.LT.1.0D0) THEN
         IF(FSOLG.GT.0.0D0) THEN         ! dendritic internal fraction
          IF(INDEX2.EQ.0) THEN           ! only dendritic solidification
           IF(TGAUST.LE.TEINF) FINTD=1.0D0
           IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .      FINTD=1.0D0-((TGAUST-TEINF)/(TESUP-TEINF))
           IF(TGAUST.GT.TESUP) FINTD=0.0D0
C
           IF(TGAUST.LE.TEINF) TSOCF=0.0D0
           IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .      TSOCF=1.0D0/(TESUP-TEINF)
           IF(TGAUST.GT.TESUP) TSOCF=0.0D0
c          TSOCF=DABS((FINTD-FINTDO)/(DTEMPT*DTIMET))           ! simpl.
          ELSE                           ! also eutectic solidification
           IF(TGAUST.LE.TESUPE)          ! instead of 1.0
     .      FINTD=1.0D0-((TESUPE-TEINF)/(TESUP-TEINF))
           IF(TGAUST.GT.TESUPE.AND.TGAUST.LE.TESUP)
     .      FINTD=1.0D0-((TGAUST-TEINF)/(TESUP-TEINF))
           IF(TGAUST.GT.TESUP) FINTD=0.0D0
C
           IF(TGAUST.LE.TESUPE) TSOCF=0.0D0
           IF(TGAUST.GT.TESUPE.AND.TGAUST.LE.TESUP)
     .      TSOCF=1.0D0/(TESUP-TEINF)
           IF(TGAUST.GT.TESUP) TSOCF=0.0D0
c          TSOCF=DABS((FINTD-FINTDO)/(DTEMPT*DTIMET))           ! simpl.
          ENDIF
          IF(FINTD.LT.FINTDO) FINTD=FINTDO ! control (no remelting)
          FINTX=FINTD+FINTDE                                   ! control
          IF(FINTX.GT.1.0D0) THEN        ! assumption valid for small dt
           FINTD=1.0D0-FINTDE
C                                        ! no control on TSOCF
          ENDIF
          ICORREC=1
         ENDIF                           ! fsolg.lt.1.0
        ENDIF                            ! gu.gt.0.0.and.fintd.lt.1.0...
       ENDIF                             ! imgro.eq.2
       IF(IMGRO.EQ.3) THEN               ! Scheil's eq.
        IF(GU.GT.0.0D0.AND.FINTD.LT.1.0D0.AND.FINTF.LT.1.0D0) THEN
         IF(FSOLG.GT.0.0D0) THEN         ! dendritic internal fraction
c         PARTX=1.0D0/(PARCO-1.0D0)      ! Exponent (see micphas.f) eq.
          PARTX=1.0D0/(DIFLI-1.0D0)      ! Exponent (see micphas.f) neq.
C
          IF(INDEX2.EQ.0) THEN           ! only dendritic solidification
           IF(TGAUST.LE.TEINF) FINTD=1.0D0
           IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .      FINTD=1.0D0-((TGAUST-TESUPF)/(TESUP-TESUPF))**PARTX
           IF(TGAUST.GT.TESUP) FINTD=0.0D0
           IF(FINTD.LT.FINTDO) FINTD=FINTDO ! control (no remelting)
C
           IF(TGAUST.LE.TEINF) TSOCF=0.0D0
           IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .      TSOCF=PARTX*(((TGAUST-TESUPF)/
     .                    (TESUP-TESUPF))**(PARTX-1.0D0))/(TESUP-TESUPF)
           IF(TGAUST.GT.TESUP) TSOCF=0.0D0
          ELSE                           ! also eutectic solidification
           IF(TGAUST.LE.TESUPE)          ! instead of 1.0
     .      FINTD=1.0D0-((TESUPE-TESUPF)/(TESUP-TESUPF))**PARTX
           IF(TGAUST.GT.TESUPE.AND.TGAUST.LE.TESUP)
     .      FINTD=1.0D0-((TGAUST-TESUPF)/(TESUP-TESUPF))**PARTX
           IF(TGAUST.GT.TESUP) FINTD=0.0D0
           IF(FINTD.LT.FINTDO) FINTD=FINTDO ! control (no remelting)
C
           IF(TGAUST.LE.TESUPE) TSOCF=0.0D0
           IF(TGAUST.GT.TESUPE.AND.TGAUST.LE.TESUP)
     .      TSOCF=PARTX*(((TGAUST-TESUPF)/
     .                    (TESUP-TESUPF))**(PARTX-1.0D0))/(TESUP-TESUPF)
           IF(TGAUST.GT.TESUP) TSOCF=0.0D0
          ENDIF                          ! index2.eq.0
          IF(FINTD.LT.FINTDO) FINTD=FINTDO ! control (no remelting)
          FINTX=FINTD+FINTDE                                   ! control
          IF(FINTX.GT.1.0D0) THEN        ! assumption valid for small dt
           FINTD=1.0D0-FINTDE
C                                        ! no control on TSOCF
          ENDIF
          ICORREC=1
         ENDIF                           ! fsolg.lt.1.0
        ENDIF                            ! gu.gt.0.0.and.fintd.lt.1.0...
       ENDIF                             ! imgro.eq.3
       IF(IMGRO.EQ.4) THEN               ! lever rule
        IF(GU.GT.0.0D0.AND.FINTD.LT.1.0D0.AND.FINTF.LT.1.0D0) THEN
         IF(FSOLG.GT.0.0D0) THEN         ! dendritic internal fraction
          TEMEL=TESUPF
          PARTX=((TESUP-TEMEL)/(TGAUST-TEMEL))  ! factor (see micphas.f)
C
          IF(INDEX2.EQ.0) THEN           ! only dendritic solidification
           IF(TGAUST.LE.TEINF) FINTD=1.0D0
           IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .      FINTD=1.0D0-((TGAUST-TEINF)/(TESUP-TEINF))*PARTX
           IF(TGAUST.GT.TESUP) FINTD=0.0D0
C
           IF(TGAUST.LE.TEINF) TSOCF=0.0D0
           IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .      TSOCF=((TEMEL-TESUP)*(TEMEL-TEINF))/
     .                     ((TGAUST-TEMEL)*(TGAUST-TEMEL)*(TESUP-TEINF))
           IF(TGAUST.GT.TESUP) TSOCF=0.0D0
c          TSOCF=DABS((FINTD-FINTDO)/(DTEMPT*DTIMET))           ! simpl.
          ELSE                           ! also eutectic solidification
           PARTE=((TESUP-TEMEL)/(TESUPE-TEMEL)) ! factor (see micphas.f)
           IF(TGAUST.LE.TESUPE)          ! instead of 1.0
     .      FINTD=1.0D0-((TESUPE-TEINF)/(TESUP-TEINF))*PARTE
           IF(TGAUST.GT.TESUPE.AND.TGAUST.LE.TESUP)
     .      FINTD=1.0D0-((TGAUST-TEINF)/(TESUP-TEINF))*PARTX
           IF(TGAUST.GT.TESUP) FINTD=0.0D0
C
           IF(TGAUST.LE.TESUPE) TSOCF=0.0D0
           IF(TGAUST.GT.TESUPE.AND.TGAUST.LE.TESUP)
     .      TSOCF=((TEMEL-TESUP)*(TEMEL-TEINF))/
     .                     ((TGAUST-TEMEL)*(TGAUST-TEMEL)*(TESUP-TEINF))
           IF(TGAUST.GT.TESUP) TSOCF=0.0D0
c          TSOCF=DABS((FINTD-FINTDO)/(DTEMPT*DTIMET))           ! simpl.
          ENDIF                         ! index2.eq.0
          IF(FINTD.LT.FINTDO) FINTD=FINTDO      ! control (no remelting)
          FINTX=FINTD+FINTDE                                   ! control
          IF(FINTX.GT.1.0D0) THEN       ! assumption valid for small dt
           FINTD=1.0D0-FINTDE
C                                       ! no control on TSOCF
          ENDIF
          ICORREC=1
         ENDIF                          ! fsolg.lt.1.0
        ENDIF                           ! gu.gt.0.0.and.fintd.lt.1.0....
       ENDIF                            ! imgro.eq.4
       IF(IMGRO.EQ.5) THEN
        call runendt('error in micros3: imgro=5 not implemented')
         IF(GU.GT.0.0D0.AND.FINTF.LT.1.0D0) THEN
          IF(FSOLG.LT.1.0D0) THEN        ! solute diffusion model
           CALL LIQCOMT(COMPL,PROPST,TGAUST,
     .                  NLICO,    2)     ! computes liquid composition
           OMEGAC=(COMPL-SICON)/(COMPL*(1.0D0-PARCO))
           ZETAG=(2.0D0*DIFLI)/(DELRA*RADGR)
           GFUNCT=1.0D0+1.5D0*ZETAG+ZETAG**2.0D0+0.25D0*ZETAG**3.0D0
           FINTF=OMEGAC*GFUNCT
           IF(FINTF.LT.FINTO) FINTF=FINTO ! control (no remelting)
           IF(FINTF.GT.1.0D0) FINTF=1.0D0 !assumption valid for small dt
C
           CALL LIQCOMT(SICON,PROPST,PENDCL,
     .                 NLICO,    3)      ! computes slope of T-C curve
           PENDCL=-PENDCL
           TSOCF=GFUNCT/                 ! approximation 1
     .     (PENDCL*COMPL*(1.0D0-PARCO))
c          TSOCF=GFUNCT/                 ! approximation 2 (Thevoz)
c    .      (PENDCL*SICON*(1.0D0-PARCO))
          ELSE                           ! Scheil's equation (fsolg=1.0)
           PARTX=1.0D0/(PARCO-1.0D0)     ! Exponent (see micphas.f)
C
           IF(TGAUST.LE.TESUPE) FINTF=1.0D0 ! form 1 (total)
           IF(TGAUST.GT.TESUPE.AND.TGAUST.LE.TESUP)
     .      FINTF=1.0D0-((TGAUST-TESUPF)/(TESUP-TESUPF))**PARTX
           IF(TGAUST.GT.TESUP) FINTF=0.0D0
           IF(FINTF.LT.FINTO) FINTF=FINTO       ! control (no remelting)
c          IF(TGAUST.LE.TESUPE) FINTF=1.0D0     ! form 2 (incremental)
c          IF(TGAUST.GT.TESUPE.AND.TGAUST.LE.TESUP)
c    .      FINTF=FINTF-PARTX*(((TGAUST-TESUPF)/
c    .             (TESUP-TESUPF))**(PARTX-1.0D0))/(TESUP-TESUPF)*DTEMPT
c          IF(TGAUST.GT.TESUP) FINTF=0.0D0
c          IF(FINTF.LT.FINTO) FINTF=FINTO ! control (no remelting)
C
           IF(TGAUST.LE.TESUPE) TSOCF=0.0D0
           IF(TGAUST.GT.TESUPE.AND.TGAUST.LE.TESUP)
     .      TSOCF=PARTX*(((TGAUST-TESUPF)/
     .                    (TESUP-TESUPF))**(PARTX-1.0D0))/(TESUP-TESUPF)
          IF(TGAUST.GT.TESUP) TSOCF=0.0D0
          ICORREC=1
         ENDIF                          ! fsolg.lt.1.0
        ENDIF                           ! gu.gt.0.0.and.fintd.lt.1.0....
       ENDIF                            ! imgro.eq.5
C
       IF(ICORREC.EQ.1) THEN            ! correction
        FINTF=FINTD+FINTDE
        FSOLD=FSOLG*FINTD
        IF(IAFLOJ.EQ.0) THEN
         FSOLGDE=FSOLG*FINTDE
         FSOLDE=FSOLGDE+FSOLGE
        ENDIF
        FSOLID=FSOLD+FSOLDE
       ENDIF                             ! icorrec.eq.1
C
       TSOCN=TSOCN*2.0D0/3.0D0*TWOPIT*RADGR**3.0D0*FINTF
       TSOCR=TSOCR*2.0D0*TWOPIT*DENNU*RADGR**2.0D0*FINTF
       TSOCF=TSOCF*2.0D0/3.0D0*TWOPIT*DENNU*RADGR**3.0D0
      ENDIF                      ! index1.eq.1
C
C**** EUTECTIC SOLIDIFICATION
C
      IF(INDEX2.EQ.1) THEN
       GUE=TESUPE-TGAUST         ! eutectic undercooling
       ICORREC=0
C
C**** NUCLEATION MODEL
C
       IF(DENNUE.GT.0.0D0.AND.DTEMPT.GT.0.0D0) INDEXE=1   ! end of nucl.
       IF(IMNUCE.EQ.1.OR.IMNUCE.EQ.2) THEN         ! instantaneous nucl.
        IF(GUE.GT.0.0D0.AND.DTEMPT.LT.0.0D0.AND.FSOLGE.LT.1.0D0
     .                 .AND.FSOLDE.LT.1.0D0.AND.FSOLID.LT.1.0D0
     .                 .AND.INDEXE.EQ.0) THEN
C
         NSUGUE=1                                    ! subincrementation
         IF(EXTGUE.GT.0.0D0) THEN
          IF(TEMAVE.GT.0.0D0) THEN
           NSUGUE=INT(GUE/(EXTGUE*TEMAVE))
          ELSE
           NSUGUE=INT(GUE/EXTGUE)
          ENDIF
          IF(NSUGUE.LT.1)       NSUGUE=1
          IF(NSUGUE.GT.NSUGUEM) NSUGUE=NSUGUEM
         ENDIF
         IF(NSUGUE.EQ.1) THEN
          AAAAA=((GUE-TEMAVE)/(2.0D0*TEMDEE))**2.0D0
          FNUCL=1.0D0/(DSQRT(TWOPIT))*DEXP(-AAAAA)
          DELNUE=-GRDEME*FNUCL*DTEMPT
          DENNUE=DENNUE+DELNUE*DTIMET
          IF(DENNUE.GT.GRDEME) DENNUE=GRDEME
c
c    ! same computation of N with error function => not working yet !!!!
c         BAAAA=(GUE-TEMAVE)/(2.0D0*TEMDEE)
c         DENNUEX=DERF(BAAAA)                           ! error function
c         DENNUEX=GRDEME*DENNUEX
c         IF(DENNUEX.GT.DENNUE) DENNUE=DENNUEX
c
         ELSE
          DTEMPX=DTEMPT/NSUGUE
          TGAUSX=TGAUAT
          DO I=1,NSUGUE
           TGAUSX=TGAUSX+DTEMPX*DTIMET
           GUEX=TESUPE-TGAUSX
           IF(GUEX.GT.0.0D0) THEN
            AAAAA=((GUEX-TEMAVE)/(2.0D0*TEMDEE))**2.0D0
            FNUCL=1.0D0/(DSQRT(TWOPIT))*DEXP(-AAAAA)
            DELNUE=-GRDEME*FNUCL*DTEMPX
            DENNUE=DENNUE+DELNUE*DTIMET
            IF(DENNUE.GT.GRDEME) DENNUE=GRDEME
           ENDIF
          ENDDO
         ENDIF
C
         FSOLGE=2.0D0/3.0D0*TWOPIT*DENNUE*RADGRE**3.0D0        ! control
         FSOLX=FSOLG                  ! only grain eutectic nuc.& growth
         IF(IMGROE.EQ.6) FSOLX=FSOLD  ! interd. & grain eutectic n. & g.
         FSOLIX=FSOLX+FSOLGE
         IF(FSOLIX.GT.1.0D0) THEN        ! assumption valid for small dt
          FSOLGE=1.0D0-FSOLX
          DENNUE=3.0D0*FSOLGE/(2.0D0*TWOPIT*RADGRE**3.0D0)
          IF(DENNUE.LT.DENNUEO) DENNUE=DENNUEO
          FNUCL=-(DENNUE-DENNUEO)/(DTIMET*GRDEME*DTEMPT)
         ENDIF
C
         TSOCNE=GRDEME*FNUCL
         ICORREC=1
        ENDIF                            ! gue.gt.0.and....
       ENDIF                             ! imnuce.eq.1.or....
C
       IF(ICORREC.EQ.1) THEN             ! correction
        IF(IAFLOJ.EQ.1) THEN             ! strict conditions
         FSOLDE=FSOLGDE+FSOLGE
         FSOLID=FSOLD+FSOLDE
        ENDIF                            ! iafloj.eq.1
       ENDIF                             ! icorrec.eq.1
C
C**** GROWTH MODEL
C
C     a) radius growth
C     b) internal fraction growth
C
       ICORREC=0                                                    ! a)
       IF(IMGROE.EQ.1.OR.IMGROE.EQ.2.OR.IMGROE.EQ.3.OR.
     .    IMGROE.EQ.4) THEN
        IF(GUE.GT.0.0D0.AND.FSOLGE.LT.1.0D0.AND.FSOLDE.LT.1.0D0
     .                 .AND.FSOLID.LT.1.0D0.AND.DENNUE.GT.0.0D0) THEN
         RADGRE=0.0D0
         FSOLGE=0.0D0
         ICORREC=1
        ENDIF                            ! gue.gt.0.0.and....
       ENDIF                             ! imgro.eq.1.or....
       IF(IMGROE.EQ.5.OR.IMGROE.EQ.6) THEN
        IF(GUE.GT.0.0D0.AND.FSOLGE.LT.1.0D0.AND.FSOLDE.LT.1.0D0
     .                 .AND.FSOLID.LT.1.0D0.AND.DENNUE.GT.0.0D0) THEN
         CALL LIQCOMT(DELRAE,PROPST,GUE,   ! computes grain velocity (v)
     .                NRAUNE,    6)
         RADGRE=RADGRE+DELRAE*DTIMET
C
         FSOLGE=2.0D0/3.0D0*TWOPIT*DENNUE*RADGRE**3.0D0        ! control
         FSOLX=FSOLG                  ! only grain eutectic nuc.& growth
         IF(IMGROE.EQ.6) FSOLX=FSOLD  ! interd. & grain eutectic n. & g.
         FSOLIX=FSOLX+FSOLGE
         IF(FSOLIX.GT.1.0D0) THEN        ! assumption valid for small dt
          FSOLGE=1.0D0-FSOLX
          RADGRE=(3.0D0*FSOLGE/(2.0D0*TWOPIT*DENNUE))**(1.0D0/3.0D0)
          IF(RADGRE.LT.RADGREO) RADGRE=RADGREO
C                                        ! no control on TSOCRE
         ENDIF
C
         CALL LIQCOMT(DEDELRAE,PROPST,GUE, ! computes dv/dT
     .                NRAUNE,    7)
         TSOCRE=DEDELRAE*DTIMET
         ICORREC=1
        ENDIF                            ! gue.gt.0.0.and....
       ENDIF                             ! imgroe.eq.5.or.imgroe.eq.6
C
       IF(ICORREC.EQ.1) THEN             ! correction
        IF(IAFLOJ.EQ.1) THEN             ! strict conditions
         FSOLDE=FSOLGDE+FSOLGE
         FSOLID=FSOLD+FSOLDE
        ENDIF                            ! iafloj.eq.1
       ENDIF                             ! icorrec.eq.1
C
       ICORREC=0                                                    ! b)
       IF(IMGROE.EQ.1.OR.IMGROE.EQ.5) THEN        ! internal fracion=0.0
        IF(GUE.GT.0.0D0.AND.FINTDE.LT.1.0D0.AND.FINTF.LT.1.0D0) THEN
         IF(FSOLG.GT.0.0D0) THEN         ! eutectic internal fraction
          FINTDE=0.0D0
          FINTX=FINTD+FINTDE                                   ! control
          IF(FINTX.GT.1.0D0) THEN        ! assumption valid for small dt
           FINTDE=1.0D0-FINTD
C                                        ! no control on TSOCFE
          ENDIF
          TSOCFE=0.0D0
          ICORREC=1
         ENDIF                           ! fsolg.lt.1.0
        ENDIF                            ! gue.gt.0.and....
       ENDIF                             ! imgroe.eq.1
       IF(IMGROE.EQ.2) THEN              ! linear model
        IF(GUE.GT.0.0D0.AND.FINTDE.LT.1.0D0.AND.FINTF.LT.1.0D0) THEN
         IF(FSOLG.GT.0.0D0) THEN         ! eutectic internal fraction
          IF(INDEX1.EQ.0) THEN           ! only eutectic solidification
           IF(TGAUST.LE.TEINF) FINTDE=1.0D0
           IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .      FINTDE=1.0D0-((TGAUST-TEINF)/(TESUP-TEINF))
           IF(TGAUST.GT.TESUP) FINTDE=0.0D0
C
           IF(TGAUST.LE.TEINF) TSOCFE=0.0D0
           IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .      TSOCFE=1.0D0/(TESUP-TEINF)
           IF(TGAUST.GT.TESUP) TSOCF=0.0D0
c          TSOCFE=DABS((FINTD-FINTDO)/(DTEMPT*DTIMET))          ! simpl.
          ELSE                           ! also dendritic solidification
           IF(TGAUST.LE.TEINF) FINTDE=1.0D0-FINTD       ! instead of 1.0
           IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUPE)
     .      FINTDE=(1.0D0-(TGAUST-TEINF)/(TESUPE-TEINF))*(1.0D0-FINTD)
           IF(TGAUST.GT.TESUPE) FINTDE=0.0D0
C
           IF(TGAUST.LE.TEINF) TSOCFE=0.0D0
           IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUPE)
     .      TSOCFE=1.0D0/(TESUPE-TEINF)*FINTD
           IF(TGAUST.GT.TESUPE) TSOCF=0.0D0
c          TSOCFE=DABS((FINTD-FINTDO)/(DTEMPT*DTIMET))          ! simpl.
          ENDIF
          IF(FINTDE.LT.FINTDEO) FINTDE=FINTDEO  ! control (no remelting)
          FINTX=FINTD+FINTDE                                   ! control
          IF(FINTX.GT.1.0D0) THEN        ! assumption valid for small dt
           FINTDE=1.0D0-FINTD
C                                        ! no control on TSOCFE
          ENDIF
          ICORREC=1
         ENDIF                           ! fsolg.lt.1.0
        ENDIF                            ! gue.gt.0.0.and.fintde.lt.1.0.
       ENDIF                             ! imgroe.eq.2
       IF(IMGROE.EQ.3) THEN              ! Scheil's eq.
        IF(GUE.GT.0.0D0.AND.FINTDE.LT.1.0D0.AND.FINTF.LT.1.0D0) THEN
         IF(FSOLG.GT.0.0D0) THEN         ! eutectic internal fraction
c         PARTX=1.0D0/(PARCO-1.0D0)      ! Exponent (see micphas.f) eq.
          PARTX=1.0D0/(DIFLI-1.0D0)      ! Exponent (see micphas.f) neq.
C
          IF(INDEX1.EQ.0) THEN           ! only eutectic solidification
           IF(TGAUST.LE.TEINF) FINTDE=1.0D0
           IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .      FINTDE=1.0D0-((TGAUST-TESUPF)/(TESUP-TESUPF))**PARTX
           IF(TGAUST.GT.TESUP) FINTDE=0.0D0
C
           IF(TGAUST.LE.TEINF) TSOCFE=0.0D0
           IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .      TSOCFE=PARTX*(((TGAUST-TESUPF)/
     .                    (TESUP-TESUPF))**(PARTX-1.0D0))/(TESUP-TESUPF)
           IF(TGAUST.GT.TESUP) TSOCFE=0.0D0
          ELSE                           ! also dendritic solidification
           IF(TGAUST.LE.TEINF) FINTDE=1.0D0-FINTD       ! instead of 1.0
           IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUPE)
     .      FINTDE=(1.0D0-((TGAUST-TESUPF)/(TESUPE-TESUPF))**PARTX)*
     .                                                     (1.0D0-FINTD)
           IF(TGAUST.GT.TESUPE) FINTDE=0.0D0
C
           IF(TGAUST.LE.TEINF) TSOCFE=0.0D0
           IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUPE)
     .      TSOCFE=PARTX*(((TGAUST-TESUPF)/
     .            (TESUPE-TESUPF))**(PARTX-1.0D0))/(TESUPE-TESUPF)*FINTD
           IF(TGAUST.GT.TESUPE) TSOCFE=0.0D0
          ENDIF
          IF(FINTDE.LT.FINTDEO) FINTDE=FINTDEO  ! control (no remelting)
          FINTX=FINTD+FINTDE                                   ! control
          IF(FINTX.GT.1.0D0) THEN        ! assumption valid for small dt
           FINTDE=1.0D0-FINTD
C                                        ! no control on TSOCFE
          ENDIF
          ICORREC=1
         ENDIF                        ! fsolg.lt.1.0
        ENDIF                         ! gue.gt.0.0.and.fintde.lt.1.0....
       ENDIF                          ! imgroe.eq.3
       IF(IMGROE.EQ.4) THEN           ! lever rule
        IF(GUE.GT.0.0D0.AND.FINTDE.LT.1.0D0.AND.FINTF.LT.1.0D0) THEN
         IF(FSOLG.GT.0.0D0) THEN         ! eutectic internal fraction
          TEMEL=TESUPF
          PARTX=((TESUP-TEMEL)/(TGAUST-TEMEL))  ! factor (see micphas.f)
C
          IF(INDEX1.EQ.0) THEN           ! only eutectic solidification
           IF(TGAUST.LE.TEINF) FINTDE=1.0D0
           IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .      FINTDE=1.0D0-((TGAUST-TEINF)/(TESUP-TEINF))*PARTX
           IF(TGAUST.GT.TESUP) FINTDE=0.0D0
C
           IF(TGAUST.LE.TEINF) TSOCFE=0.0D0
           IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .      TSOCFE=((TEMEL-TESUP)*(TEMEL-TEINF))/
     .                     ((TGAUST-TEMEL)*(TGAUST-TEMEL)*(TESUP-TEINF))
           IF(TGAUST.GT.TESUP) TSOCFE=0.0D0
c          TSOCFE=DABS((FINTDE-FINTDEO)/(DTEMPT*DTIMET))        ! simpl.
          ELSE                           ! also dendritic solidification
           PARTC=((TESUPE-TEMEL)/(TGAUST-TEMEL)) ! factor(see micphas.f)
           IF(TGAUST.LE.TEINF) FINTDE=1.0D0-FINTD       ! instead of 1.0
           IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUPE)
     .      FINTDE=(1.0D0-((TGAUST-TEINF)/(TESUPE-TEINF))*PARTC)*
     .                                                     (1.0D0-FINTD)
           IF(TGAUST.GT.TESUPE) FINTDE=0.0D0
C
           IF(TGAUST.LE.TEINF) TSOCFE=0.0D0
           IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUPE)
     .      TSOCFE=((TEMEL-TESUPE)*(TEMEL-TEINF))/
     .              ((TGAUST-TEMEL)*(TGAUST-TEMEL)*(TESUPE-TEINF))*FINTD
           IF(TGAUST.GT.TESUPE) TSOCFE=0.0D0
c          TSOCFE=DABS((FINTDE-FINTDEO)/(DTEMPT*DTIMET))        ! simpl.
          ENDIF
          IF(FINTDE.LT.FINTDEO) FINTDE=FINTDEO  ! control (no remelting)
          FINTX=FINTD+FINTDE                                   ! control
          IF(FINTX.GT.1.0D0) THEN        ! assumption valid for small dt
           FINTDE=1.0D0-FINTD
C                                        ! no control on TSOCFE
          ENDIF
          ICORREC=1
         ENDIF                        ! fsolg.lt.1.0
        ENDIF                         ! gue.gt.0.0.and.fintde.lt.1.0....
       ENDIF                          ! imgro.eq.4
       IF(IMGROE.EQ.6) THEN         ! interd. & grain eutectic n. & g.
        IF(GUE.GT.0.0D0.AND.FSOLD.LT.1.0D0.AND.FSOLGE.GT.0.0D0) THEN
         FSOLGX=FSOLGE              ! total obtained eutectic
         F3=1.0D0-FSOLD             ! total expected eutectic at the end
         F1=(1.0D0-FSOLG)/F3        ! proportional to grain
         F2=(FSOLG-FSOLD)/F3        ! proportional to interdendritic
C
         FSOLGE=FSOLGX*F1
         FSOLGDE=FSOLGX*F2
C
         IF(FSOLG.GT.0.0D0) FINTDE=FSOLGDE/FSOLG
         FINTX=FINTD+FINTDE                                    ! control
         IF(FINTX.GT.1.0D0) THEN         ! assumption valid for small dt
          FINTDE=1.0D0-FINTD
C                                        ! no control on TSOCFE
         ENDIF
         TSOCFE=0.0D0
         ICORREC=1
        ENDIF                            ! gue.gt.0.and....
       ENDIF                             ! imgroe.eq.6
C
       IF(ICORREC.EQ.1) THEN             ! correction
        FINTF=FINTD+FINTDE
        FSOLGDE=FSOLG*FINTDE
        FSOLDE=FSOLGDE+FSOLGE
        FSOLID=FSOLD+FSOLDE
       ENDIF                             ! icorrec.eq.1
C
       TSOCNE=TSOCNE*2.0D0/3.0D0*TWOPIT*RADGRE**3.0D0
       TSOCRE=TSOCRE*2.0D0*TWOPIT*DENNUE*RADGRE**2.0D0
       TSOCFE=TSOCFE*2.0D0/3.0D0*TWOPIT*DENNUE*RADGRE**3.0D0
      ENDIF                      ! index2.eq.1
C
      FSOLD=FSOLG*FINTD          ! always
      FSOLGDE=FSOLG*FINTDE
      FSOLDE=FSOLGDE+FSOLGE
      FSOLID=FSOLD+FSOLDE
C
C**** LIQUID FRACTION
C
      FLIQD=1.0D0-FSOLID                   ! liquid fraction
C
C**** DEFINES THE "MACROSCOPICAL" PROPERTIES (DENS., CAPACITY & CONDUC.)
C
c     BASMM=DENS                           ! no change
c     BASCC=ESPH
c     BASKK(1)=COND
C
C**** ESTABLISHES LATENT HEAT * PHASE-CHANGE FUNCTION (at time t & t+dt)
C
      IF(HENER.EQ.HENERE) THEN              ! classical
       TSOE1(IPLAT)=FLIQD                   ! f_pc at time t+dt
       TSOE2(IPLAT)=FLIQDO                  ! f_pc at time t
       TSOE2(IPLAT)=TSOE2(IPLAT)*HENER      ! L*f_pc at time t
       TSOE1(IPLAT)=TSOE1(IPLAT)*HENER      ! L*f_pc at time t+dt
      ELSE                                         ! improved
       TSOE2(IPLAT)=-HENER*FSOLDO-HENERE*FSOLDEO   ! L*f_pc at time t
       TSOE1(IPLAT)=-HENER*FSOLD -HENERE*FSOLDE    ! L*f_pc at time t+dt
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
C
C     Notes:
C
C     L is constant
C
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
      IF(HENER.EQ.HENERE) THEN              ! classical
       TSOC1(IPLAT)=TSOCN+TSOCR+TSOCF+TSOCNE+TSOCRE+TSOCFE
       TSOC1(IPLAT)=TSOC1(IPLAT)*HENER
      ELSE                                  ! improved
       TSOC1(IPLAT)=(TSOCN+TSOCR+TSOCF)*HENER+
     .              (TSOCNE+TSOCRE+TSOCFE)*HENERE
      ENDIF
      GO TO 10
C
    4 ICACOT=1
      IF(HENER.EQ.HENERE) THEN              ! classical
       TSOC1(IPLAT)=TSOCN+TSOCR+TSOCF+TSOCNE+TSOCRE+TSOCFE
       TSOC1(IPLAT)=TSOC1(IPLAT)*HENER
      ELSE                                  ! improved
       TSOC1(IPLAT)=(TSOCN+TSOCR+TSOCF)*HENER+
     .              (TSOCNE+TSOCRE+TSOCFE)*HENERE
      ENDIF
      GO TO 10
C
    5 ICACOT=1
      TSOC1(IPLAT)=0.0D0
      GO TO 10
C
   10 CONTINUE
C
C**** TRANSFERS MICROSCOPICAL VARIABLES TO "ALPHAM"
C
      ALPHAM(IN+1)=FLIQD
C
      ALPHAM(IN+2)=FSOLD
      ALPHAM(IN+3)=FSOLG
      ALPHAM(IN+4)=FINTD
      ALPHAM(IN+5)=DENNU
      ALPHAM(IN+6)=RADGR
C
      ALPHAM(IN+7)=FSOLDE
      ALPHAM(IN+8)=FSOLGE
      ALPHAM(IN+9)=FSOLGDE
      ALPHAM(IN+10)=FINTDE
      ALPHAM(IN+11)=DENNUE
      ALPHAM(IN+12)=RADGRE
C
      ALPHAM(IN+13)=FLOAT(INDEXD)
      ALPHAM(IN+14)=FLOAT(INDEXE)
C
C**** INCREMENTS "ALPHAM" INDEX
C
      INUPC=INUPC+14                       ! less than 19; see pointes.f
C
      RETURN
      END
