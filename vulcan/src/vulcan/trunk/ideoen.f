      SUBROUTINE IDEOEN(PROPS)
C***********************************************************************
C
C**** THIS ROUTINE ORDER THE PROPERTIES FOR THE ORTHOTROPIC, ELASTIC &
C     NON TEMPERATURE-DEPENDENT MODEL
C
C***********************************************************************
C
C     VYOUN(1) =Young's modulus, component xx
C     VYOUNF(1)=Young's modulus, component yy
C     VYOUNA(1)=Young's modulus, component zz
C
C     VPOIS(1) =Poisson's ratio, component xy
C     VPOISF(1)=Poisson's ratio, component xz
C     VPOISA(1)=Poisson's ratio, component yz
C
C     VSHEA(1) =Shear modulus, component xy
C     VSHEAF(1)=Shear modulus, component xz
C     VSHEAA(1)=Shear modulus, component yz
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
C**** YOUNG MODULI
C
      I1=61
      VYOUN(1,1)=PROPS(I1)
      I1=I1+1
      VYOUNF(1,1)=PROPS(I1)
      I1=I1+1
      VYOUNA(1,1)=PROPS(I1)
C
C**** POISSON RATII
C
      I2=I1+1
      VPOIS(1,1)=PROPS(I2)
      I2=I2+1
      VPOISF(1,1)=PROPS(I2)
      I2=I2+1
      VPOISA(1,1)=PROPS(I2)
C
C**** SHEAR MODULI
C
      I3=I2+1
      VSHEA(1,1)=PROPS(I3)
      IF(NDIME.EQ.3) THEN
       I3=I3+1
       VSHEAF(1,1)=PROPS(I3)
       I3=I3+1
       VSHEAA(1,1)=PROPS(I3)
      ENDIF
C
      I4=I3+1
      VCCER(1,1)=PROPS(I4)                     ! infinite value for C^th
C
      IVIFL=0
      ISOTT=0
      IKINE=0
      IRECR=0
      IPORO=0
      IDAMG=0
C
      I5=I4+1
      IFREN=INT(PROPS(I5))
C
      RETURN
      END
