      SUBROUTINE IDEIEN(PROPS)
C***********************************************************************
C
C**** THIS ROUTINE ORDER THE PROPERTIES FOR THE ISOTROPIC, ELASTIC &
C     NON TEMPERATURE-DEPENDENT MODEL
C
C***********************************************************************
C
C     VYOUN(1)=Young's Modulus
C
C     VPOIS(1)=Poisson's Ratio
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
      I3=I2+1
      VCCER(1,1)=PROPS(I3)                     ! infinite value for C^th
C
      IVIFL=0
      ISOTT=0
      IKINE=0
      IRECR=0
      IPORO=0
      IDAMG=0
C
      I4=I3+1
      IFREN=INT(PROPS(I4))
C
      RETURN
      END
