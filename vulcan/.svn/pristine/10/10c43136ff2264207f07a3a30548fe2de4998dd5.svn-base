      SUBROUTINE INMICRO4(ALPHAM,INUPC,IPLAT)
C***********************************************************************
C
C**** THIS ROUTINE INITIALITES MICROSTRUCTURAL PARAMETERS FOR MODEL 4
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
      CA =VPLAT(IPLAT, 6)          ! carbon
      SIO=VPLAT(IPLAT, 7)          ! silicon (initial)
      PH =VPLAT(IPLAT, 8)          ! phosporus
      CU =VPLAT(IPLAT, 9)          ! copper
      QMO=VPLAT(IPLAT,10)          ! manganese (initial)
      QG =VPLAT(IPLAT,11)          ! magnesium
      CR =VPLAT(IPLAT,12)          ! chromium
C
      AKSI=VPLAT(IPLAT,13)         ! silicon partition coeffic.
      AKQM=VPLAT(IPLAT,14)         ! manganese partition coeffic.
      IHELAC=INT(VPLAT(IPLAT,15))  ! T & C correl. (Heine, Lacaze, etc)
      IX0=0
      IF(IHELAC.EQ.4) THEN         ! Thermocalc
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
C**** INITIALITES MICROSTRUCTURAL VARIABLES ("ALPHAM")
C
      IN=INUPC
      I2=2+INNUM4T
      ALPHAM(IN+1)=1.0D0     ! assumed liquid state as initial condition
      ALPHAM(IN+2)=0.0D0     ! assumed zero initial austenite fraction
      ALPHAM(IN+3)=0.0D0     ! assumed zero initial graphite fraction
      ALPHAM(IN+4)=SIO       ! silicon content (initial)
      I5A=5
      IF(IMNMT.EQ.1) THEN
       ALPHAM(IN+5)=QMO      ! manganese content (initial)
       I5A=6
      ENDIF
      I5=I5A+NNUM4T-1        ! NNUM4T different grain densities
      DO INU=I5A,I5
       ALPHAM(IN+INU)=0.0D0
      ENDDO
      I6A=I5+1
      I6B=I6A+(1+INNUM4T)*NNUM4T-1
      DO INU=I6A,I6B         ! NNUM4T different grain radii
       ALPHAM(IN+INU)=0.0D0
      ENDDO
      IF(IWEUT.EQ.1) THEN    ! white fraction computation
       I7=I6B+1
       DO INU=I7,I7+2
        ALPHAM(IN+INU)=0.0D0
       ENDDO
      ENDIF
      I10=I6B+3*IWEUT+1
      ALPHAM(IN+I10)=0.0D0   ! number of families
      I11=I10+1
      IF(INUCAX.EQ.1)
     . ALPHAM(IN+I11)=0.0D0  ! index: nucleation up to Trec
      IF(INUCAX.EQ.2)
     . ALPHAM(IN+I11)=2.0D3  ! index: nucleation up to Tmin
C
C**** INCREMENTS "ALPHAM" INDEX
C
      INUPC=INUPC+4+IMNMT+I2*NNUM4T+   ! less than NBASES; see pointes.f
     .      3*IWEUT+2
C
      RETURN
      END
