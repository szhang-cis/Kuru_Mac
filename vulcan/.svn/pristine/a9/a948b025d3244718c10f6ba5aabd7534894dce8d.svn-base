      SUBROUTINE POISMT(TEMPG,PROPS,
     .                  POISM,POISMF,POISMA)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THERMAL POISSON RATIO COEFFICIENT
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
      IF(IPEP1.EQ.100) THEN                 ! isotropic
       IPEP2=INT(PROPS(2))
       IF(IPEP2.EQ.40) RETURN               ! hyperelastic
C
       IPEP3=INT(PROPS(3))                  ! no temp.-dependent
       IF(IPEP3.EQ.1) THEN
        POISM=VPOIS(1,1)
        RETURN
       ENDIF
C
#ifndef restricted
       IPEP4=INT(PROPS(4))
       IF(IPEP4.EQ.0.OR.IPEP4.EQ.2) THEN    ! standard or green sand
        IF(TEMPG.LE.VPOIS(1,2)) POISM=VPOIS(1,1)
        DO IPOIS=2,NPOIS
         I1=IPOIS-1
         I2=IPOIS
         IF(TEMPG.GE.VPOIS(I1,2).AND.TEMPG.LE.VPOIS(I2,2)) 
     .    POISM=(VPOIS(I2,1)-VPOIS(I1,1))/(VPOIS(I2,2)-VPOIS(I1,2))*
     .          (TEMPG-VPOIS(I1,2))+VPOIS(I1,1)
        ENDDO
        IF(TEMPG.GE.VPOIS(NPOIS,2)) POISM=VPOIS(NPOIS,1)
C
        DPOIS=0.0D0                         ! temperature derivative
       ENDIF                ! ipep4.eq.0
C
       IF(IPEP4.EQ.1.OR.IPEP4.EQ.3) THEN    ! non-standard: SGCI or DP
        IF(TEMPG.LE.VPOISF(1,2)) POISMF=VPOISF(1,1)
        DO IPOIS=2,NPOISF
         I1=IPOIS-1
         I2=IPOIS
         IF(TEMPG.GE.VPOISF(I1,2).AND.TEMPG.LE.VPOISF(I2,2))
     .    POISMF=(VPOISF(I2,1)-VPOISF(I1,1))/
     .           (VPOISF(I2,2)-VPOISF(I1,2))*
     .           (TEMPG-VPOISF(I1,2))+VPOISF(I1,1)
        ENDDO
        IF(TEMPG.GE.VPOISF(NPOISF,2)) POISMF=VPOISF(NPOISF,1)
C
        DPOISF=0.0D0                        ! temperature derivative
C
        IF(TEMPG.LE.VPOISA(1,2)) POISMA=VPOISA(1,1)
        DO IPOIS=2,NPOISA
         I1=IPOIS-1
         I2=IPOIS
         IF(TEMPG.GE.VPOISA(I1,2).AND.TEMPG.LE.VPOISA(I2,2))
     .    POISMA=(VPOISA(I2,1)-VPOISA(I1,1))/
     .           (VPOISA(I2,2)-VPOISA(I1,2))*
     .           (TEMPG-VPOISA(I1,2))+VPOISA(I1,1)
        ENDDO
        IF(TEMPG.GE.VPOISA(NPOISA,2)) POISMA=VPOISA(NPOISA,1)
C
        DPOISA=0.0D0                        ! temperature derivative
       ENDIF                ! ipep4.eq.1
      ENDIF                 ! ipep1.eq.100
C
      IF(IPEP1.EQ.200) THEN                 ! orthotropic
       IPEP2=INT(PROPS(2))
       IF(IPEP2.EQ.40) RETURN               ! hyperelastic
C
       IPEP3=INT(PROPS(3))
       IF(IPEP3.EQ.1) THEN                  ! no temp.-dependent
        POISM=VPOIS(1,1)
        POISMF=VPOISF(1,1)
        POISMA=VPOISA(1,1)
        RETURN
       ENDIF
C
       IPEP4=INT(PROPS(4))
       IF(IPEP4.EQ.0) THEN                  ! standard
        IF(TEMPG.LE.VPOIS(1,2)) POISM=VPOIS(1,1)
        DO IPOIS=2,NPOIS
         I1=IPOIS-1
         I2=IPOIS
         IF(TEMPG.GE.VPOIS(I1,2).AND.TEMPG.LE.VPOIS(I2,2)) 
     .    POISM=(VPOIS(I2,1)-VPOIS(I1,1))/(VPOIS(I2,2)-VPOIS(I1,2))*
     .          (TEMPG-VPOIS(I1,2))+VPOIS(I1,1)
        ENDDO
        IF(TEMPG.GE.VPOIS(NPOIS,2)) POISM=VPOIS(NPOIS,1)
C
        DPOIS=0.0D0                         ! temperature derivative
C
        IF(TEMPG.LE.VPOISF(1,2)) POISMF=VPOISF(1,1)
        DO IPOIS=2,NPOISF
         I1=IPOIS-1
         I2=IPOIS
         IF(TEMPG.GE.VPOISF(I1,2).AND.TEMPG.LE.VPOISF(I2,2))
     .    POISMF=(VPOISF(I2,1)-VPOISF(I1,1))/
     .           (VPOISF(I2,2)-VPOISF(I1,2))*
     .           (TEMPG-VPOISF(I1,2))+VPOISF(I1,1)
        ENDDO
        IF(TEMPG.GE.VPOISF(NPOISF,2)) POISMF=VPOISF(NPOISF,1)
C
        DPOISF=0.0D0                        ! temperature derivative
C
        IF(TEMPG.LE.VPOISA(1,2)) POISMA=VPOISA(1,1)
        DO IPOIS=2,NPOISA
         I1=IPOIS-1
         I2=IPOIS
         IF(TEMPG.GE.VPOISA(I1,2).AND.TEMPG.LE.VPOISA(I2,2))
     .    POISMA=(VPOISA(I2,1)-VPOISA(I1,1))/
     .           (VPOISA(I2,2)-VPOISA(I1,2))*
     .           (TEMPG-VPOISA(I1,2))+VPOISA(I1,1)
        ENDDO
        IF(TEMPG.GE.VPOISA(NPOISA,2)) POISMA=VPOISA(NPOISA,1)
C
        DPOISA=0.0D0                        ! temperature derivative
       ENDIF                ! ipep4.eq.0
#endif
      ENDIF                 ! ipep1.eq.200
C
      RETURN
      END
