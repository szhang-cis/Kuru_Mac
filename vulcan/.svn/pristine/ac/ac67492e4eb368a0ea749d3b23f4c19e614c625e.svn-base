      SUBROUTINE IDEPROS10(PROPST,IPLAT,IA1)
C***********************************************************************
C
C**** THIS ROUTINE ORDERS THE MICROSTRUCTURAL PROPERTIES OF MODEL
C     NUMBER 10 (IPCMO=10) OF RATE PHASE-CHANGE FORMULATIONS
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
      IA2=IA1+2         ! 2=ipcfo,ipcmo
C
      CTOK=PROPST(IA2+1)
C
      OX =  PROPST(IA2+2)
      PH =  PROPST(IA2+3)
      AS =  PROPST(IA2+4)
      SB =  PROPST(IA2+5)
      FE =  PROPST(IA2+6)
      SE =  PROPST(IA2+7)
      TE =  PROPST(IA2+8)
      SU =  PROPST(IA2+9)
      ZN =  PROPST(IA2+10)
      PB =  PROPST(IA2+11)
C
      TLD = PROPST(IA2+12) 
      EU =  PROPST(IA2+13)
      OP =  PROPST(IA2+14)
C
      INUCMX=INT(PROPST(IA2+15))
      IF(INUCMX.EQ.1) THEN       ! Turnbull & Fischer model
       ANUCA=PROPST(IA2+16)      ! atoms per cubic meter
       ANUCB=PROPST(IA2+17)      ! atoms in l surrounding surface of s
       ANUCC=PROPST(IA2+18)      ! vibration frecuency
       ANUCD=PROPST(IA2+19)      ! s-l interfacial energy
       ANUCF=PROPST(IA2+20)      ! enthalpy of fusion
       ANUCG=PROPST(IA2+21)      ! Boltzman constant
       ANUCH=PROPST(IA2+22)      ! f(theta) parameter
       ICX=7
      END IF
      IF(INUCMX.EQ.2) THEN       ! Holloman & Turnbull model
       ANUCI=PROPST(IA2+16)      ! B1 coefficient
       ANUCJ=PROPST(IA2+17)      ! s-l interfacial energy
       ANUCL=PROPST(IA2+18)      ! enthalpy of fusion
       ANUCM=PROPST(IA2+19)      ! Boltzman constant
       ANUCN=PROPST(IA2+20)      ! f(theta) parameter
       ICX=5
      END IF
      IF(INUCMX.EQ.3) THEN       ! Oldfield model
       ANUCQ=PROPST(IA2+16)      ! A parameter
       ANUCR=PROPST(IA2+17)      ! n exponent 
       ICX=2
      ENDIF
C
      IGROMX=INT(PROPST(IA2+16+ICX))
      IF(IGROMX.EQ.1) THEN       ! Rappaz & Thevoz growth model 1
       GROCA=PROPST(IA2+17+ICX)  ! diffusion coef. of Cu in liquid
       GROCB=PROPST(IA2+18+ICX)  ! Gibbs-Thomson parameter for copper
       GROCC=PROPST(IA2+19+ICX)  ! liquidus curve slope
       IGX=3
      ENDIF
      IF(IGROMX.EQ.2) THEN       ! Pero - Sanz growth model
       GROCF=PROPST(IA2+17+ICX)  ! Z coefficient
       IGX=1
      ENDIF
      IF(IGROMX.EQ.3) THEN       ! Turnbull growth model 1
       GROCK=PROPST(IA2+17+ICX)  ! enthalpy of fusion
       GROCL=PROPST(IA2+18+ICX)  ! gases constant
       GROCM=PROPST(IA2+19+ICX)  ! factor order 1
       GROCN=PROPST(IA2+20+ICX)  ! sound speed
       IGX=4
      ENDIF
      IF(IGROMX.EQ.4) THEN       ! Turnbull growth model 2
       GROCP=PROPST(IA2+17+ICX)  ! enthalpy of fusion
       GROCQ=PROPST(IA2+18+ICX)  ! gases constant
       GROCR=PROPST(IA2+19+ICX)  ! Boltzman constant
       GROCS=PROPST(IA2+20+ICX)  ! atomic mass
       IGX=4
      ENDIF
      IF(IGROMX.EQ.5) THEN       ! Rappaz & Thevoz growth model 2
       GROCA=PROPST(IA2+17+ICX)  ! diffusion coef. of Cu in liquid
       GROCB=PROPST(IA2+18+ICX)  ! Gibbs-Thomson parameter for copper
       GROCC=PROPST(IA2+19+ICX)  ! liquidus curve slope
       GROCD=PROPST(IA2+20+ICX)  ! solute partition coefficient
       GROCE=PROPST(IA2+21+ICX)  ! initial composition (oxygen)
       DIFOL=PROPST(IA2+22+ICX)  ! diffusion coef. of O in liquid
       RTOTR=PROPST(IA2+23+ICX)  ! maximum observed grain radius
       IGX=7
      ENDIF
C
      IFPCDT=INT(PROPST(IA2+17+ICX+IGX))
      IAFLOJ=INT(PROPST(IA2+18+ICX+IGX))
C
      VPLAT(IPLAT, 6)=CTOK
C
      VPLAT(IPLAT, 7)=OX
      VPLAT(IPLAT, 8)=PH
      VPLAT(IPLAT, 9)=AS
      VPLAT(IPLAT,10)=SB
      VPLAT(IPLAT,11)=FE
      VPLAT(IPLAT,12)=SE
      VPLAT(IPLAT,13)=TE
      VPLAT(IPLAT,14)=SU   
      VPLAT(IPLAT,15)=ZN
      VPLAT(IPLAT,16)=PB
C
      VPLAT(IPLAT,17)=TLD
      VPLAT(IPLAT,18)=EU
      VPLAT(IPLAT,19)=OP
C
      VPLAT(IPLAT,20)=FLOAT(INUCMX)
      IF(INUCMX.EQ.1) THEN       ! Turnbull & Fischer model
       VPLAT(IPLAT,21)=ANUCA
       VPLAT(IPLAT,22)=ANUCB
       VPLAT(IPLAT,23)=ANUCC
       VPLAT(IPLAT,24)=ANUCD
       VPLAT(IPLAT,25)=ANUCF
       VPLAT(IPLAT,26)=ANUCG
       VPLAT(IPLAT,27)=ANUCH
      END IF
      IF(INUCMX.EQ.2) THEN       ! Holloman & Turnbull model
       VPLAT(IPLAT,21)=ANUCI
       VPLAT(IPLAT,22)=ANUCJ
       VPLAT(IPLAT,23)=ANUCL
       VPLAT(IPLAT,24)=ANUCM
       VPLAT(IPLAT,25)=ANUCN
      ENDIF
      IF(INUCMX.EQ.3) THEN       ! Oldfield model
       VPLAT(IPLAT,21)=ANUCQ
       VPLAT(IPLAT,22)=ANUCR
      ENDIF
C
      VPLAT(IPLAT,21+ICX)=FLOAT(IGROMX)
      IF(IGROMX.EQ.1) THEN       ! Rappaz & Thevoz growth model 1
       VPLAT(IPLAT,22+ICX)=GROCA
       VPLAT(IPLAT,23+ICX)=GROCB
       VPLAT(IPLAT,24+ICX)=GROCC
      ENDIF
      IF(IGROMX.EQ.2) THEN       ! Pero - Sanz growth model
       VPLAT(IPLAT,22+ICX)=GROCF
      ENDIF
      IF(IGROMX.EQ.3) THEN       ! Turnbull growth model 1
       VPLAT(IPLAT,22+ICX)=GROCK
       VPLAT(IPLAT,23+ICX)=GROCL
       VPLAT(IPLAT,24+ICX)=GROCM
       VPLAT(IPLAT,25+ICX)=GROCN
      ENDIF
      IF(IGROMX.EQ.4) THEN       ! Turnbull growth model 2
       VPLAT(IPLAT,22+ICX)=GROCP
       VPLAT(IPLAT,23+ICX)=GROCQ
       VPLAT(IPLAT,24+ICX)=GROCR
       VPLAT(IPLAT,25+ICX)=GROCS
      ENDIF
      IF(IGROMX.EQ.5) THEN       ! Rappaz & Thevoz growth model 2
       VPLAT(IPLAT,22+ICX)=GROCA
       VPLAT(IPLAT,23+ICX)=GROCB
       VPLAT(IPLAT,24+ICX)=GROCC
       VPLAT(IPLAT,25+ICX)=GROCD
       VPLAT(IPLAT,26+ICX)=GROCE
       VPLAT(IPLAT,27+ICX)=DIFOL
       VPLAT(IPLAT,28+ICX)=RTOTR
      ENDIF
C
      VPLAT(IPLAT,22+ICX+IGX)=FLOAT(IFPCDT)
      VPLAT(IPLAT,23+ICX+IGX)=FLOAT(IAFLOJ)
C
      IMODE=18+ICX+IGX         ! imode=total number of prop. of model 10
C
      IA1=IA2+IMODE
C
      RETURN
      END
