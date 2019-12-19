      SUBROUTINE IDEFA1(PROPS)
C***********************************************************************
C
C**** THIS ROUTINE ORDER THE FATIGUE PARAMETERS FOR MODEL=1
C
C***********************************************************************
C
C     VFATI(1)=limit number of cycles
C
C     VFATI(2)=A parameter
C
C     VFATI(3)=B parameter
C
C     VFATI(4)=C parameter
C
C     VFATI(5)=D parameter
C
C     VFATI(6)=minimum reversibility factor
C
C     VFATI(7)=maximum reversibility factor
C
C     VFATI(8)=E parameter
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
C**** FATIGUE PARAMETERS
C
      I1=42
      VFATI(1)=PROPS(I1)
C
      I1=I1+1
      VFATI(2)=PROPS(I1)
C
      I1=I1+1
      VFATI(3)=PROPS(I1)
C
      I1=I1+1
      VFATI(4)=PROPS(I1)
C
      I1=I1+1
      VFATI(5)=PROPS(I1)
C
      I1=I1+1
      VFATI(6)=PROPS(I1)
C
      I1=I1+1
      VFATI(7)=PROPS(I1)
C
      I1=I1+1
      VFATI(8)=PROPS(I1)
C
      RETURN
      END
