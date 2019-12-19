      SUBROUTINE INMICRO2(ALPHAM,INUPC,IPLAT)
C***********************************************************************
C
C**** THIS ROUTINE INITIALITES MICROSTRUCTURAL PARAMETERS FOR MODEL 2
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
      CA=VPLAT(IPLAT,6)                    ! carbon
      SI=VPLAT(IPLAT,7)                    ! silicon
      PH=VPLAT(IPLAT,8)                    ! phosporus
      CU=VPLAT(IPLAT,9)                    ! copper
      QM=VPLAT(IPLAT,10)                   ! manganese
      QG=VPLAT(IPLAT,11)                   ! magnesium
      CR=VPLAT(IPLAT,12)                   ! chromium
      AN=VPLAT(IPLAT,13)                   ! nickel
C
      DENS=VPLAT(IPLAT,14)
      ESPH=VPLAT(IPLAT,15)
      COND=VPLAT(IPLAT,16)
      HENER=VPLAT(IPLAT,17)
C
      IAUSM=INT(VPLAT(IPLAT,18))
      AUSTS=VPLAT(IPLAT,19)
      AUSTL=VPLAT(IPLAT,20)
      AUSLA=VPLAT(IPLAT,21)
C
      ANUCG=VPLAT(IPLAT,22)
      ANUCC=VPLAT(IPLAT,23)
      BGROG=VPLAT(IPLAT,24)
      BGROC=VPLAT(IPLAT,25)
C
      AK1LA=VPLAT(IPLAT,26)
      AK2LA=VPLAT(IPLAT,27)
      MAXNGL=INT(VPLAT(IPLAT,28))
C
      IFPCDT=INT(VPLAT(IPLAT,29))
      IAFLOJ=INT(VPLAT(IPLAT,30))
C
      IV=30                      ! number of VPLAT defined in input data
C
C**** INITIALITES MICROSTRUCTURAL PARAMETERS ("ALPHAM" & OTHERS)
C
C     Notes:
C
C     Relationships based on thermodynamics by:
C     -Heine R.: AFS Trans. A, 1986, vol. 94, pp. 391-402.
C     -Glover D., Bates C. and Monroe R.: AFS Trans., 1982, vol. 90,a
C                                         pp. 747-757.
C
C     included in:
C     -Fras E. and Lopez H.: A Generalized Theory of the Chilling
C      Tendency of Cast Iron, AFS Transactions, Vol. 132, pp. 355-363 
C      (1993).
C     -Goettsch D. and Dantzig J.: Modeling Microstructure Development
C      in Gray Cast Irons, Metallurgical and Materials Transactions A,
C      Vol. 25A, 1063-1079 (1994).
C     -Maijer D., Cockcroft S. and Patt W.: Mathematical Modeling of
C      Microstructural Development in Hypoeutectic Cast Iron,
C      Metallurgical and Materials Transactions A, Vol. 30A, 2147-2158
C      (1999).
C
C     Tecnometal relationship: useful for white hypoeutectic alloys
C                              (Rodrigo Bermudez's thesis - 2001)
C
C
C     Nomenclature:
C     TL:    liquidus temperature (hypo or hypereutectic)
C     TS:    solidus temperature (hypo or hypereutectic) => not used !
C     TEX:   temperature corresponding to the intersection of the
C            extrapolation of the liquid and solid temperature curves
C     TEG:   equilibrium graphite eutectic temperature
C     TEC:   "equilibrium" cementite eutectic temperature
C     CEQ:   equivalent carbon content
C     CEUT:  carbon concentration at eutectic point 
C     CGMAX: maximum carbon concentration in austenite
C     PK:    equilibrium asutenite distribution (partition) coefficient
C     REAG:  austenite volume/graphite volume (eutectic)
C
      TL= 1569.00D0- 97.3D0*(CA+0.25D0*SI)               ! hypoeutectic
c     TL= 1669.00D0-124.0D0*(CA+0.25D0*SI+0.5D0*PH)      ! other option
c     TL= 1636.00D0-113.0D0*(CA+0.25D0*SI+0.5D0*PH)      ! other option
c     TL= 1585.00D0- 97.3D0*(CA+0.25D0*SI)               ! Tecnometal
c     TL=  389.10D0*(CA+0.31D0*SI+0.33D0*PH-1.3D0)       ! hypereutectic
      TL1= 389.10D0*(CA+SI/3.0D0)-503.20D0               ! other option
      TL2= 389.10D0*(CA+0.31D0*SI)-505.80D0              ! other option
C
      TS= 1528.40D0-177.9D0*(CA+0.18D0*SI)
C
      TEX=1618.00D0-15.0332D0*SI
C
      TEG=1135.06D0+13.89D0*SI-2.05D0*SI*SI              ! hypo/hyper
c     TEG=1153.90D0+5.25D0*SI-14.88*PH                   ! other option
C
      TEC=1138.20D0-6.93D0*(SI+2.5D0*PH)-                ! hypo/hyper
     .    1.717D0*(SI+2.5D0*PH)*(SI+2.5D0*PH)
c     TEC=1138.20D0-6.93D0*(SI+2.5D0*PH)-                ! other option
c    .                               1.717D0*(SI+2.5D0*PH)*(SI+2.5D0*PH)
c     TEC=1134.90D0-11.02D0*SI-39.7D0*PH                 ! other option
c     TEC=1149.20D0-6.93D0*(SI+2.5D0*CU)-                ! Tecnometal
c    .                               1.717D0*(SI+2.5D0*PH)*(SI+2.5D0*PH)
C
      CEQ=CA+0.33D0*(SI+PH)                  ! equivalent carbon content
C
      IF(CEQ.GE.4.26D0.AND.CEQ.LT.4.40D0) TL=TL1
      IF(CEQ.GE.4.40D0.AND.CEQ.LT.4.60D0) TL=TL2
C
      CEUT=4.3D0-0.37D0*SI+0.02D0*SI*SI-0.5D0*PH
      CGMAX=2.2D0-0.26*SI-0.01D0*SI*SI
C
      PK=(2.1D0-0.2165D0*SI)/(4.26D0-0.3167D0*SI)        ! hypoeutectic
c     PK=CGMAX/CEUT                                      ! other option
C
      REAG=(100.0D0-CEUT)/(CEUT-CGMAX)*2300.0D0/7510.0D0
c     REAG=(26.26D0+0.087D0*SI)/(2.16D0-0.101D0*SI)      ! other option
      FGREUT=1.0D0/(1.0D0+REAG)
C
      VPLAT(IPLAT,IV+1)=TL
      VPLAT(IPLAT,IV+2)=TS
      VPLAT(IPLAT,IV+3)=TEX
      VPLAT(IPLAT,IV+4)=TEG
      VPLAT(IPLAT,IV+5)=TEC
      VPLAT(IPLAT,IV+6)=PK
      VPLAT(IPLAT,IV+7)=CEQ
      VPLAT(IPLAT,IV+8)=FGREUT
C
      IN=INUPC
      ALPHAM(IN+1)=1.0D0  ! liquid state is assumed as initial condition
      DO INU=2,14
       ALPHAM(IN+INU)=0.0D0
      ENDDO
C
C**** INCREMENTS "ALPHAM" INDEX
C
      INUPC=INUPC+14                   ! less than NBASES; see pointes.f
C
      RETURN
      END
