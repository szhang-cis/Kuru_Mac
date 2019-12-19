      SUBROUTINE IDEIPN(PROPS)
C***********************************************************************
C
C**** THIS ROUTINE ORDER THE PROPERTIES FOR THE ISOTROPIC, PLASTIC &
C     NON TEMPERATURE-DEPENDENT MODEL
C
C***********************************************************************
C
C     VYOUN(1)=Young's Modulus
C
C     VPOIS(1)=Poisson's Ratio
C
C     VALPH(1)=Thermal Dilatation Coefficient
C
C     VCCER(1)=Hardening Function
C
C
C     ISOTT=indicator of isotropic hardening model
C
C     VCCOE(1)=Hardening Coefficient
C
C     VCCOB(1)=Coefficient b
C
C     VCCOQ(1)=Coefficient Q
C
C     VCCOEE(1)=Hardening Coefficient         VCCOEE(2)=Temperature
C
C     IPORO=indicator of porosity model (Gurson)
C
C     VPCOE(1)=Porosity Coefficient
C
C     VPCO1(1)=Q1 Coefficient
C
C     VPCO2(1)=Q2 Coefficient
C
C     TEREF=Reference Temperature
C
C     TEMPS=Solidus Temperature input in mechanical data file
C     TEMPL=Liquidus Temperature input in mechanical data file
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION PROPS(*)
C
      NCRIT=INT(PROPS(36))
      NCRIP=INT(PROPS(52))
      NASNA=INT(PROPS(47))
C
C**** YOUNG MODULUS
C
      I1=61
      VYOUN(1,1)=PROPS(I1)
C
C**** POISSON RATIO
C
      I2=I1+1
      VPOIS(1,1)=PROPS(I2)
C
C**** REFERENCE TEMPERATURE
C
      TEREF=PROPS(11)        ! see propipn.f
C
C**** THERMAL DILATATION
C
      I3=I2+1
      VALPH(1,1)=PROPS(I3)
C
C**** HARDENING FUNCTION
C
      I4=I3+1
      VCCER(1,1)=PROPS(I4)
C
C**** ISOTROPIC HARDENING INDICATOR
C
      I5=I4+1
      ISOTT=INT(PROPS(I5))
c     IF(ISOTT.EQ.0) THEN
       I6=I5
c     ENDIF
C
C**** HARDENING COEFFICIENT
C
      IF(ISOTT.EQ.1.OR.ISOTT.EQ.2.OR.ISOTT.EQ.4.OR.ISOTT.EQ.5.OR.
     .   ISOTT.EQ.8.OR.ISOTT.EQ.10.OR.ISOTT.EQ.11.OR.ISOTT.EQ.12.OR.
     .   ISOTT.EQ.13) THEN
       I6=I5+1
       VCCOE(1,1)=PROPS(I6)
      ENDIF
C
C**** B COEFFICIENT
C
      IF(ISOTT.EQ.3.OR.ISOTT.EQ.6.OR.ISOTT.EQ.7.OR.ISOTT.EQ.8.OR.
     .   ISOTT.EQ.10.OR.ISOTT.EQ.11.OR.ISOTT.EQ.12.OR.ISOTT.EQ.13) THEN
c      I6=I5+1
       I6=I6+1
       VCCOB(1,1)=PROPS(I6)
C
C**** Q COEFFICIENT
C
       I6A=I6+1
       VCCOQ(1,1)=PROPS(I6A)
       I6=I6A
      ENDIF
C
C**** HARDENING COEFFICIENT (Ep & T dependent)
C
      IF(ISOTT.EQ.9) THEN
       I6=I5+1
       NCCOEE=INT(PROPS(I6))
C
       DO ICCOEE=1,NCCOEE
        IA=I6-1+2*ICCOEE
        IB=IA+1
        VCCOEE(ICCOEE,1)=PROPS(IA)
        VCCOEE(ICCOEE,2)=PROPS(IB)
       ENDDO
       I6=I6+2*NCCOEE
      ENDIF
C
C**** KINEMATIC HARDENING INDICATOR
C
      I7=I6+1
      IKINE=INT(PROPS(I7))
      IF(IKINE.EQ.0) THEN
       I8=I7
      ENDIF
C
C**** HARDENING COEFFICIENT
C
      IF(IKINE.EQ.1.OR.IKINE.EQ.2) THEN
       I8=I7+1
       VKCOE(1,1)=PROPS(I8)
      ENDIF
C
C**** B COEFFICIENT
C
      IF(IKINE.EQ.3.OR.IKINE.EQ.4) THEN
       I8=I7+1
       VKCOB(1,1)=PROPS(I8)
C
C**** Q COEFFICIENT
C
       I8A=I8+1
       VKCOQ(1,1)=PROPS(I8A)
C
       I8=I8A
      ENDIF
C
      IVIFL=0
      IRECR=0
C
C**** POROSITY INDICATOR
C
      IPORO=0
      IF(NCRIT.EQ.41.OR.NCRIP.EQ.41) THEN
       I7=I8+1
       IPORO=INT(PROPS(I7))
C
C**** POROSITY COEFFICIENT
C
       IF(IPORO.EQ.1.OR.IPORO.EQ.2.OR.IPORO.EQ.3.OR.IPORO.EQ.4) THEN
        I8=I7+1
        VPCOE(1,1)=PROPS(I8)
       ENDIF
C
       IF(IPORO.EQ.4) THEN
        I8=I8+1
        VPCO1(1,1)=PROPS(I8)
C
        I8=I8+1
        VPCO2(1,1)=PROPS(I8)
       ENDIF
      ENDIF
C
C**** SOLIDUS & LIQUIDUS TEMPERATURES
C
      I9=I8+1
C
      TEMPS=PROPS(I9)
      VPLATM(1,8)=PROPS(I9)
      TEMPL=PROPS(I9+1)
      VPLATM(1,9)=PROPS(I9+1)
C
C**** DAMAGE INDICATOR
C
      I10=I9+1+1
      IDAMG=INT(PROPS(I10))
      I11=I10
C
C**** PARAMETER_A
C
      IF(IDAMG.EQ.1) THEN
       I11=I10+1
       VDAMA(1,1)=PROPS(I11)
      ENDIF
C
C**** PARAMETERS: EPSILON_D, EPSILON_R & CRITICAL DAMAGE
C
      IF(IDAMG.EQ.4.OR.IDAMG.EQ.5) THEN
       I11=I10+1
       VDAMA(1,1)=PROPS(I11)
       I11=I11+1
       VFRAC(1,1)=PROPS(I11)
       I11=I11+1
       VCOE9(1,1)=PROPS(I11)
       I11=I11+1
       VCOE10(1,1)=PROPS(I11)
       I11=I11+1
       VCOE11(1,1)=PROPS(I11)
      ENDIF
C
C**** PARAMETERS: S_PARAMETER, S_EXPONENT, B_EXPONENT, A_PARAMETERS &
C     CRITICAL DAMAGE
C
      IF(IDAMG.EQ.6) THEN
       I11=I10+1
       VDAMA(1,1)=PROPS(I11)
       I11=I11+1
       VFRAC(1,1)=PROPS(I11)
       I11=I11+1
       VCOE9(1,1)=PROPS(I11)
       I11=I11+1
       VCOE10(1,1)=PROPS(I11)
       I11=I11+1
       VCOE11(1,1)=PROPS(I11)
       I11=I11+1
       VCOE12(1,1)=PROPS(I11)
      ENDIF
C
C**** FREE ENERGY MODEL
C
      I12=I11+1
      IFREN=INT(PROPS(I12))
C
      RETURN
      END
