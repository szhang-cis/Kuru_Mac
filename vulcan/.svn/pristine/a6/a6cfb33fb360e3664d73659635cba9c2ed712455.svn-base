      SUBROUTINE SHEART(TEMPG,PROPS,
     .                  SHEAR,SHEARF,SHEARA)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE SHEAR MODULUS COEFFICIENT (TEMPERATURE
C     DEPENDENT)
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION PROPS(*)
C
      ISOTR=INT(PROPS(1))   ! 0=isotropic; 1=orthotropic
      IPEP1=100
      IF(ISOTR.EQ.1) IPEP1=200
C
      IF(IPEP1.EQ.100) RETURN               ! isotropic
C
#ifndef restricted
      IF(IPEP1.EQ.200) THEN                 ! orthotropic
       IPEP2=INT(PROPS(2))
       IF(IPEP2.EQ.40) RETURN               ! hyperelastic
C
       IPEP3=INT(PROPS(3))
       IF(IPEP3.EQ.1) THEN                  ! no temp.-dependent
        SHEAR=VSHEA(1,1)
        IF(NDIME.EQ.3) THEN
         SHEARF=VSHEAF(1,1)
         SHEARA=VSHEAA(1,1)
        ENDIF
        RETURN
       ENDIF
C
       IPEP4=INT(PROPS(4))
       IF(IPEP4.EQ.0) THEN                  ! standard
        IF(TEMPG.LE.VSHEA(1,2)) SHEAR=VSHEA(1,1)
        DO ISHEA=2,NSHEA
         I1=ISHEA-1
         I2=ISHEA
         IF(TEMPG.GE.VSHEA(I1,2).AND.TEMPG.LE.VSHEA(I2,2))
     .    SHEAR=(VSHEA(I2,1)-VSHEA(I1,1))/(VSHEA(I2,2)-VSHEA(I1,2))*
     .          (TEMPG-VSHEA(I1,2))+VSHEA(I1,1)
        ENDDO
        IF(TEMPG.GE.VSHEA(NSHEA,2)) SHEAR=VSHEA(NSHEA,1)
C
        DSHEA=0.0D0                         ! temperature derivative
C
        IF(NDIME.EQ.3) THEN
         IF(TEMPG.LE.VSHEAF(1,2)) SHEARF=VSHEAF(1,1)
         DO ISHEA=2,NSHEAF
          I1=ISHEA-1
          I2=ISHEA
          IF(TEMPG.GE.VSHEAF(I1,2).AND.TEMPG.LE.VSHEAF(I2,2))
     .     SHEARF=(VSHEAF(I2,1)-VSHEAF(I1,1))/
     .            (VSHEAF(I2,2)-VSHEAF(I1,2))*
     .            (TEMPG-VSHEAF(I1,2))+VSHEAF(I1,1)
         ENDDO
         IF(TEMPG.GE.VSHEAF(NSHEAF,2)) SHEARF=VSHEAF(NSHEAF,1)
C
         DSHEAF=0.0D0                       ! temperature derivative
C
         IF(TEMPG.LE.VSHEAA(1,2)) SHEARA=VSHEAA(1,1)
         DO ISHEA=2,NSHEAA
          I1=ISHEA-1
          I2=ISHEA
          IF(TEMPG.GE.VSHEAA(I1,2).AND.TEMPG.LE.VSHEAA(I2,2))
     .     SHEARA=(VSHEAA(I2,1)-VSHEAA(I1,1))/
     .            (VSHEAA(I2,2)-VSHEAA(I1,2))*
     .            (TEMPG-VSHEAA(I1,2))+VSHEAA(I1,1)
         ENDDO
         IF(TEMPG.GE.VSHEAA(NSHEAA,2)) SHEARA=VSHEAA(NSHEAA,1)
C
         DSHEAA=0.0D0                       ! temperature derivative
        ENDIF
       ENDIF                ! ipep4.eq.0
      ENDIF                 ! ipep1.eq.200
#endif
C
      RETURN
      END
