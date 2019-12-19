      SUBROUTINE ALPHAT(TEMPG,PROPS,
     .                  ALPHA,ALPHAF,ALPHAA)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THERMAL EXPANSION COEFFICIENT (TEMPERATURE
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
      IF(IPEP1.EQ.100) THEN                 ! isotropic
       IPEP2=INT(PROPS(2))
       IF(IPEP2.EQ.40) RETURN               ! hyperelastic
C
       IPEP3=INT(PROPS(3))
       IF(IPEP3.EQ.1) THEN
        ALPHA=VALPH(1,1)
        RETURN
       ENDIF
C
       IPEP4=INT(PROPS(4))
       IF(IPEP4.EQ.0.OR.IPEP4.EQ.2) THEN    ! standard or green sand
        NALPU=INT(VALPU(1,2))
        IF(NALPU.EQ.0) THEN                 ! secant
         IF(TEMPG.LE.VALPH(1,2)) ALPHA=VALPH(1,1)
         DO IALPH=2,NALPH
          I1=IALPH-1
          I2=IALPH
          IF(TEMPG.GE.VALPH(I1,2).AND.TEMPG.LE.VALPH(I2,2)) 
     .       ALPHA=(VALPH(I2,1)-VALPH(I1,1))/(VALPH(I2,2)-VALPH(I1,2))*
     .             (TEMPG-VALPH(I1,2))+VALPH(I1,1)
         ENDDO
         IF(TEMPG.GE.VALPH(NALPH,2)) ALPHA=VALPH(NALPH,1)
C
         DALPH=0.0D0                        ! temperature derivative
        ELSE                                ! secant derived from tang.
         IF(TEMPG.LE.VALPH(1,2)) ALPHA=VALPU(1,1)
         DO IALPH=2,NALPH
          I1=IALPH-1
          I2=IALPH
          IF(TEMPG.GE.VALPH(I1,2).AND.TEMPG.LE.VALPH(I2,2))
     .       ALPHA=(VALPU(I2,1)-VALPU(I1,1))/(VALPH(I2,2)-VALPH(I1,2))*
     .             (TEMPG-VALPH(I1,2))+VALPU(I1,1)
         ENDDO
         IF(TEMPG.GE.VALPH(NALPH,2)) ALPHA=VALPU(NALPH,1)
C
         DALPH=0.0D0                        ! temperature derivative
        ENDIF               ! nalpu.eq.0
       ENDIF                ! ipep4.eq.0
C
       IF(IPEP4.EQ.1.OR.IPEP4.EQ.3) THEN    ! non-standard: SGCI or DP
        NALPUF=INT(VALPUF(1,2))
        IF(NALPUF.EQ.0) THEN                ! secant
         IF(TEMPG.LE.VALPHF(1,2)) ALPHAF=VALPHF(1,1)
         DO IALPH=2,NALPHF
          I1=IALPH-1
          I2=IALPH
          IF(TEMPG.GE.VALPHF(I1,2).AND.TEMPG.LE.VALPHF(I2,2))
     .     ALPHAF=(VALPHF(I2,1)-VALPHF(I1,1))/
     .            (VALPHF(I2,2)-VALPHF(I1,2))*
     .            (TEMPG-VALPHF(I1,2))+VALPHF(I1,1)
         ENDDO
         IF(TEMPG.GE.VALPHF(NALPHF,2)) ALPHAF=VALPHF(NALPHF,1)
C
         DALPHF=0.0D0                       ! temperature derivative
        ELSE                                ! secant derived from tang.
         IF(TEMPG.LE.VALPHF(1,2)) ALPHAF=VALPUF(1,1)
         DO IALPH=2,NALPHF
          I1=IALPH-1
          I2=IALPH
          IF(TEMPG.GE.VALPHF(I1,2).AND.TEMPG.LE.VALPHF(I2,2))
     .     ALPHAF=(VALPUF(I2,1)-VALPUF(I1,1))/
     .            (VALPHF(I2,2)-VALPHF(I1,2))*
     .            (TEMPG-VALPHF(I1,2))+VALPUF(I1,1)
         ENDDO
         IF(TEMPG.GE.VALPHF(NALPHF,2)) ALPHAF=VALPUF(NALPHF,1)
C
         DALPHF=0.0D0                       ! temperature derivative
        ENDIF               ! nalpuf.eq.0
C
        NALPUA=INT(VALPUA(1,2))
        IF(NALPUA.EQ.0) THEN                ! secant
         IF(TEMPG.LE.VALPHA(1,2)) ALPHAA=VALPHA(1,1)
         DO IALPH=2,NALPHA
          I1=IALPH-1
          I2=IALPH
          IF(TEMPG.GE.VALPHA(I1,2).AND.TEMPG.LE.VALPHA(I2,2))
     .     ALPHAA=(VALPHA(I2,1)-VALPHA(I1,1))/
     .            (VALPHA(I2,2)-VALPHA(I1,2))*
     .            (TEMPG-VALPHA(I1,2))+VALPHA(I1,1)
         ENDDO
         IF(TEMPG.GE.VALPHA(NALPHA,2)) ALPHAA=VALPHA(NALPHA,1)
C
         DALPHA=0.0D0                       ! temperature derivative
        ELSE                                ! secant derived from tang.
         IF(TEMPG.LE.VALPHA(1,2)) ALPHAA=VALPUA(1,1)
         DO IALPH=2,NALPHA
          I1=IALPH-1
          I2=IALPH
          IF(TEMPG.GE.VALPHA(I1,2).AND.TEMPG.LE.VALPHA(I2,2))
     .     ALPHAA=(VALPUA(I2,1)-VALPUA(I1,1))/
     .            (VALPHA(I2,2)-VALPHA(I1,2))*
     .            (TEMPG-VALPHA(I1,2))+VALPUA(I1,1)
         ENDDO
         IF(TEMPG.GE.VALPHA(NALPHA,2)) ALPHAA=VALPUA(NALPHA,1)
         DALPHA=0.0D0                       ! temperature derivative
        ENDIF               ! nalpua.eq.0
       ENDIF                ! ipep4.eq.1
      ENDIF                 ! ipep1.eq.100
C
      IF(IPEP1.EQ.200) THEN                 ! orthotropic
       IPEP2=INT(PROPS(2))
       IF(IPEP2.EQ.40) RETURN               ! hyperelastic
C
       IPEP3=INT(PROPS(3))
       IF(IPEP3.EQ.1) THEN                  ! no temp.-dependent
        ALPHA=VALPH(1,1)
        ALPHAF=VALPHF(1,1)
        ALPHAA=VALPHA(1,1)
        RETURN
       ENDIF
C
       IPEP4=INT(PROPS(4))
       IF(IPEP4.EQ.0) THEN                  ! standard
C
        NALPU=INT(VALPU(1,2))
        IF(NALPU.EQ.0) THEN                 ! secant
         IF(TEMPG.LE.VALPH(1,2)) ALPHA=VALPH(1,1)
         DO IALPH=2,NALPH
          I1=IALPH-1
          I2=IALPH
          IF(TEMPG.GE.VALPH(I1,2).AND.TEMPG.LE.VALPH(I2,2))
     .       ALPHA=(VALPH(I2,1)-VALPH(I1,1))/(VALPH(I2,2)-VALPH(I1,2))*
     .             (TEMPG-VALPH(I1,2))+VALPH(I1,1)
         ENDDO
         IF(TEMPG.GE.VALPH(NALPH,2)) ALPHA=VALPH(NALPH,1)
C
         DALPH=0.0D0                        ! temperature derivative
C
         IF(TEMPG.LE.VALPHF(1,2)) ALPHAF=VALPHF(1,1)
         DO IALPH=2,NALPHF
          I1=IALPH-1
          I2=IALPH
          IF(TEMPG.GE.VALPHF(I1,2).AND.TEMPG.LE.VALPHF(I2,2))
     .     ALPHAF=(VALPHF(I2,1)-VALPHF(I1,1))/
     .            (VALPHF(I2,2)-VALPHF(I1,2))*
     .            (TEMPG-VALPHF(I1,2))+VALPHF(I1,1)
         ENDDO
         IF(TEMPG.GE.VALPHF(NALPHF,2)) ALPHAF=VALPHF(NALPHF,1)
C
         DALPHF=0.0D0                       ! temperature derivative
C
         IF(TEMPG.LE.VALPHA(1,2)) ALPHAA=VALPHA(1,1)
         DO IALPH=2,NALPHA
          I1=IALPH-1
          I2=IALPH
          IF(TEMPG.GE.VALPHA(I1,2).AND.TEMPG.LE.VALPHA(I2,2))
     .     ALPHAA=(VALPHA(I2,1)-VALPHA(I1,1))/
     .            (VALPHA(I2,2)-VALPHA(I1,2))*
     .            (TEMPG-VALPHA(I1,2))+VALPHA(I1,1)
         ENDDO
         IF(TEMPG.GE.VALPHA(NALPHA,2)) ALPHAA=VALPHA(NALPHA,1)
C
         DALPHA=0.0D0                       ! temperature derivative
        ENDIF               ! nalpu.eq.0
       ENDIF                ! ipep4.eq.0
      ENDIF                 ! ipep1.eq.200
C
      RETURN
      END
