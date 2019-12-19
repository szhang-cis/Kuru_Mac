      SUBROUTINE CONCOFT(BASKK,PROPST,TGAUST,PSEUDO)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE CONDUCTIVITY COEFFICIENT
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
C
C     Heat flux models:
C
C     IFREKT=1  input data: conductivity coefficient at spatial config.
C     IFREKT=2  Armero's model for heat flux
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION PROPST(*), BASKK(NDIMETO)
C
C**** COMPUTES THE CONDUCTIVITY COEFFICIENT
C
      DO I=1,NDIMETO
       IF(TGAUST.LE.VCOND(1,2,I)) BASKK(I)=VCOND(1,1,I)
       DO ICOND=2,NCOND
        I1=ICOND-1
        I2=ICOND
        IF(TGAUST.GT.VCOND(I1,2,I).AND.TGAUST.LE.VCOND(I2,2,I))
     .     BASKK(I)=(VCOND(I2,1,I)-VCOND(I1,1,I))/
     .              (VCOND(I2,2,I)-VCOND(I1,2,I))*
     .              (TGAUST-VCOND(I1,2,I))+VCOND(I1,1,I)
       ENDDO
       IF(TGAUST.GT.VCOND(NCOND,2,I)) BASKK(I)=VCOND(NCOND,1,I)
      ENDDO                 ! i=1,ndimeto
C
C**** DEALS WITH FILLING MATERIAL
C
      IF(IFILL.EQ.0) RETURN
C
C**** COMPUTES THE AIR CONDUCTIVITY COEFFICIENT
C
      DO I=1,NDIMETO
       IF(TGAUST.LE.VCONDFI(1,2,I)) BASKKFI=VCONDFI(1,1,I)
       DO ICOND=2,NCONDFI
        I1=ICOND-1
        I2=ICOND
        IF(TGAUST.GT.VCONDFI(I1,2,I).AND.TGAUST.LE.VCONDFI(I2,2,I))
     .     BASKKFI=(VCONDFI(I2,1,I)-VCONDFI(I1,1,I))/
     .             (VCONDFI(I2,2,I)-VCONDFI(I1,2,I))*
     .             (TGAUST-VCONDFI(I1,2,I))+VCONDFI(I1,1,I)
       ENDDO
       IF(TGAUST.GT.VCONDFI(NCONDFI,2,I)) BASKKFI=VCONDFI(NCONDFI,1,I)
C
C**** COMPUTES THE WEIGHTED CONDUCTIVITY COEFFICIENT
C
       BASKK(I)=PSEUDO*BASKK(I)+(1.0D0-PSEUDO)*BASKKFI
      ENDDO                 ! i=1,ndimeto
C
      RETURN
      END
