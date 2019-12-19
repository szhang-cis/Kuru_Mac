      SUBROUTINE INMICRO5(ALPHAM,INUPC,IPLAT)
C***********************************************************************
C
C**** THIS ROUTINE INITIALITES MICROSTRUCTURAL PARAMETERS FOR MODEL 5
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
      QM =VPLAT(IPLAT,10)          ! manganese
      QG =VPLAT(IPLAT,11)          ! magnesium
C
      AKSI=VPLAT(IPLAT,12)         ! silicon partition coeffic.
      IHELAC=VPLAT(IPLAT,13)       ! T & C correlations (Heine & Lacaze)
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
      IAFLOJ=INT(VPLAT(IPLAT,27+IKAUX))  ! index for fraction correction
C
      IV=27+IKAUX                ! number of VPLAT defined in input data
C
C**** INITIALITES MICROSTRUCTURAL VARIABLES ("ALPHAM")
C
      IN=INUPC
      ALPHAM(IN+1)=1.0D0     ! assumed liquid state as initial condition
      ALPHAM(IN+2)=0.0D0     ! assumed zero initial austenite fraction
      ALPHAM(IN+3)=0.0D0     ! assumed zero initial graphite fraction
      ALPHAM(IN+4)=SIO       ! silicon content (initial)
      DO INU=5,5+(NNUM4T-1)  ! NNUM4T different grain densities
       ALPHAM(IN+INU)=0.0D0
      ENDDO
      DO INU=6+(NNUM4T-1),
     .       6+2*(NNUM4T-1)  ! NNUM4T different graphite nod. radii
       ALPHAM(IN+INU)=0.0D0
      ENDDO
      DO INU=7+2*(NNUM4T-1),
     .       7+3*(NNUM4T-1)  ! NNUM4T different aust. shell radii
       ALPHAM(IN+INU)=0.0D0
      ENDDO
      ALPHAM(IN+8+3*(NNUM4T-1))=0.0D0  ! graphite nucleation index
      IF(INUCAX.EQ.1)
     . ALPHAM(IN+9+3*(NNUM4T-1))=0.0D0 ! index: nucleation up to Trecal
      IF(INUCAX.EQ.2)
     . ALPHAM(IN+9+3*(NNUM4T-1))=2.0D3 ! index: nucleation up to Tmin
C
C**** INCREMENTS "ALPHAM" INDEX
C
      INUPC=INUPC+4+3*NNUM4T+2         ! less than NBASES; see pointes.f
C
      RETURN
      END
