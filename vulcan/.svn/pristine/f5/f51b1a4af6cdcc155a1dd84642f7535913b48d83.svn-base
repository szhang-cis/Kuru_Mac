      SUBROUTINE INMI05T(PROPST,EHISTT)
C***********************************************************************
C
C**** THIS ROUTINE INITIALITES SOME MICROSTRUCTURAL PARAMETERS
C     ( FOR ELEMENT NO. 5 )
C
C***********************************************************************
C
C     Index of variables:
C
C     EHISTT(   1) = Density
C     EHISTT(   2) = Specific Heat coefficient
C     EHISTT(   3) = Isotropic Conductivity or Conduct. x (orthot. mat.)
C     EHISTT(   4) = Conductivity y (orthotropic material)
C     EHISTT(   5) = Conductivity z (orthotropic material)
C     EHISTT(3:11) = Conductivity for the fully anisotropc mat. (3D)
C     EHISTT(4+IX) = L*Phase-change function
C     EHISTT(5+IX) = L*Phase-change function rate
C     EHISTT(6+IX) = Initial density
C     EHISTT(7+IX) = Coupling coefficient
C     EHISTT(8+IX) = Temperature derivative of phase-change function
C
C     EHISTT(NHISTT-NHISTM+1:NHISTT)=array of microstructural variables
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES       ! thermal-microstructural
C
      INCLUDE 'nued_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION PROPST(*), EHISTT(NHISTT,*)
C
C**** LOOP OVER INTEGRATION POINTS
C
      DO IGAUST=1,NGAULT
C
C**** INITIALIZATION OF SOME PARAMETERS
C
       CALL INMICRO(PROPST,EHISTT(NHISTT-NHISTM+1,IGAUST))
C
      END DO ! IGAUST=1,NGAULT
C
      RETURN
      END
