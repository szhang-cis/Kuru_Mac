      SUBROUTINE INMICRO10(ALPHAM,INUPC,IPLAT)
C***********************************************************************
C
C**** THIS ROUTINE INITIALITES MICROSTRUCTURAL PARAMETERS FOR MODEL 10
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
      CTOK=VPLAT(IPLAT,6)          ! C to K conversion term
C
      OX =VPLAT(IPLAT, 7)          ! oxygen (initial)
      PH =VPLAT(IPLAT, 8)          ! phosphorus
      AS =VPLAT(IPLAT, 9)          ! arsenical
      SB =VPLAT(IPLAT,10)          ! antimony
      FE =VPLAT(IPLAT,11)          ! iron
      SE =VPLAT(IPLAT,12)          ! selenium
      TE =VPLAT(IPLAT,13)          ! telurium
      SU =VPLAT(IPLAT,14)          ! sulphur
      ZN =VPLAT(IPLAT,15)          ! zinc
      PB =VPLAT(IPLAT,16)          ! lead
C
      TLD=VPLAT(IPLAT,17)          ! liquidus temperature  
      EU =VPLAT(IPLAT,18)          ! eutectic temperature
      OP =VPLAT(IPLAT,19)          ! oxygen partition coefficient
C
      INUCMX=INT(VPLAT(IPLAT,20))
      IF(INUCMX.EQ.1) THEN         ! Turnbull & Fischer model
       ANUCA=VPLAT(IPLAT,21)       ! atoms per cubic meter
       ANUCB=VPLAT(IPLAT,22)       ! atoms in l surrounding surface of s
       ANUCC=VPLAT(IPLAT,23)       ! vibration frecuency
       ANUCD=VPLAT(IPLAT,24)       ! s-l interfacial energy
       ANUCF=VPLAT(IPLAT,25)       ! enthalpy of fusion
       ANUCG=VPLAT(IPLAT,26)       ! Boltzman constant
       ANUCH=VPLAT(IPLAT,27)       ! f(theta) parameter
       ICX=7
      ENDIF
      IF(INUCMX.EQ.2) THEN         ! Holloman & Turnbull model
       ANUCI=VPLAT(IPLAT,21)       ! B1 coefficient
       ANUCJ=VPLAT(IPLAT,22)       ! s-l interfacial energy
       ANUCL=VPLAT(IPLAT,23)       ! enthalpy of fusion
       ANUCM=VPLAT(IPLAT,24)       ! Boltzman constant
       ANUCN=VPLAT(IPLAT,25)       ! f(theta) parameter
       ICX=5
      ENDIF
      IF(INUCMX.EQ.3) THEN         ! Oldfield model
       ANUCQ=VPLAT(IPLAT,21)       ! A parameter
       ANUCR=VPLAT(IPLAT,22)       ! n coefficient
       ICX=2
      ENDIF
C
      IGROMX=INT(VPLAT(IPLAT,21+ICX))
      IF(IGROMX.EQ.1) THEN         ! Rappaz & Thevoz growth model 1
       GROCA=VPLAT(IPLAT,22+ICX)   ! diffusion coef. of Cu in liquid
       GROCB=VPLAT(IPLAT,23+ICX)   ! Gibbs-Thomson parameter for copper
       GROCC=VPLAT(IPLAT,24+ICX)   ! liquidus curve slope
       IGX=3
      ENDIF
      IF(IGROMX.EQ.2) THEN         ! Pero - Sanz growth model
       GROCF=VPLAT(IPLAT,22+ICX)   ! Z coefficient
       IGX=1
      ENDIF
      IF(IGROMX.EQ.3) THEN         ! Turnbull growth model 1
       GROCK=VPLAT(IPLAT,22+ICX)   ! enthalpy of fusion
       GROCL=VPLAT(IPLAT,23+ICX)   ! gases constant
       GROCM=VPLAT(IPLAT,24+ICX)   ! factor order 1
       GROCN=VPLAT(IPLAT,25+ICX)   ! sound speed
       IGX=4
      ENDIF
      IF(IGROMX.EQ.4) THEN         ! Turnbull growth model 2
       GROCP=VPLAT(IPLAT,22+ICX)   ! enthalpy of fusion
       GROCQ=VPLAT(IPLAT,23+ICX)   ! gases constant
       GROCR=VPLAT(IPLAT,24+ICX)   ! Bolztman constant
       GROCS=VPLAT(IPLAT,25+ICX)   ! atomic mass
       IGX=4
      ENDIF
      IF(IGROMX.EQ.5) THEN         ! Rappaz & Thevoz growth model 2
       GROCA=VPLAT(IPLAT,22+ICX)   ! diffusion coef. of Cu in liquid
       GROCB=VPLAT(IPLAT,23+ICX)   ! Gibbs-Thomson parameter for copper
       GROCC=VPLAT(IPLAT,24+ICX)   ! liquidus curve slope
       GROCD=VPLAT(IPLAT,25+ICX)   ! solute partition coefficient
       GROCE=VPLAT(IPLAT,26+ICX)   ! initial composition (oxygen)
       DIFOL=VPLAT(IPLAT,27+ICX)   ! diffusion coef. of O in liquid
       RTOTR=VPLAT(IPLAT,28+ICX)   ! maximum observed grain radius
       IGX=7
      ENDIF
C
      IFPCDT=INT(VPLAT(IPLAT,22+ICX+IGX))  ! index for temp. derivative
      IAFLOJ=INT(VPLAT(IPLAT,23+ICX+IGX))  ! index for fract. correction
C
      IV=23+ICX+IGX              ! number of VPLAT defined in input data
C
C**** INITIALITES MICROSTRUCTURAL VARIABLES ("ALPHAM")
C
      IN=INUPC
      ALPHAM(IN+1)=1.000D0   ! assumed liquid state as initial condition
      ALPHAM(IN+2)=0.0D0     ! assumed zero initial dendritic solid fraction
      ALPHAM(IN+3)=0.0D0     ! assumed zero initial eutectic solid fraction
      ALPHAM(IN+4)=0.0D0     ! nuclei density (initial)
      ALPHAM(IN+5)=0.0D0     ! radii of nuclei (initial)
      ALPHAM(IN+6)=0.0D0     ! nucleation index
      ALPHAM(IN+7)=0.0D0     ! dendritic internal solid fraction
      ALPHAM(IN+8)=0.0D0     ! eutectic internal solid fraction
      ALPHAM(IN+9)=0.0D0     ! radial globular fraction (FSOPC)
      ALPHAM(IN+10)=0.0174D0 ! OPE,sólo prueba
      ALPHAM(IN+11)=0.0D0    ! TSOCCP, sólo prueba
      ALPHAM(IN+12)=0.0D0    ! TSOCCE, sólo prueba
C
C**** INCREMENTS "ALPHAM" INDEX
C
      INUPC=INUPC+12                   ! less than NBASES; see pointes.f
C
      RETURN
      END
