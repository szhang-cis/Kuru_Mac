      SUBROUTINE CONCOFS(BASKK,PROPSS,TGAUSS,PSEUDOS)
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
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      DIMENSION PROPSS(*), BASKK(NDIMETOS)
C
C**** COMPUTES THE CONDUCTIVITY COEFFICIENT
C
      DO I=1,NDIMETOS
       IF(TGAUSS.LE.VCONDS(1,2,I)) BASKK(I)=VCONDS(1,1,I)
       DO ICOND=2,NCONDS
        I1=ICOND-1
        I2=ICOND
        IF(TGAUSS.GT.VCONDS(I1,2,I).AND.TGAUSS.LE.VCONDS(I2,2,I))
     .     BASKK(I)=(VCONDS(I2,1,I)-VCONDS(I1,1,I))/
     .              (VCONDS(I2,2,I)-VCONDS(I1,2,I))*
     .              (TGAUSS-VCONDS(I1,2,I))+VCONDS(I1,1,I)
       ENDDO
       IF(TGAUSS.GT.VCONDS(NCONDS,2,I)) BASKK(I)=VCONDS(NCONDS,1,I)
      ENDDO                 ! i=1,ndimetos
C
C**** DEALS WITH FILLING MATERIAL
C
      IF(IFILLS.EQ.0) RETURN

      call runends('ifills=1 not implemented - concofs')

C
C**** COMPUTES THE AIR CONDUCTIVITY COEFFICIENT
C
c     DO I=1,NDIMETO
c      IF(TGAUST.LE.VCONDFI(1,2,I)) BASKKFI=VCONDFI(1,1,I)
c      DO ICOND=2,NCONDFI
c       I1=ICOND-1
c       I2=ICOND
c       IF(TGAUST.GT.VCONDFI(I1,2,I).AND.TGAUST.LE.VCONDFI(I2,2,I))
c    .     BASKKFI=(VCONDFI(I2,1,I)-VCONDFI(I1,1,I))/
c    .             (VCONDFI(I2,2,I)-VCONDFI(I1,2,I))*
c    .             (TGAUST-VCONDFI(I1,2,I))+VCONDFI(I1,1,I)
c      ENDDO
c      IF(TGAUST.GT.VCONDFI(NCONDFI,2,I)) BASKKFI=VCONDFI(NCONDFI,1,I)
C
C**** COMPUTES THE WEIGHTED CONDUCTIVITY COEFFICIENT
C
c      BASKK(I)=PSEUDO*BASKK(I)+(1.0-PSEUDO)*BASKKFI
c     ENDDO                 ! i=1,ndimeto
C
      RETURN
      END
