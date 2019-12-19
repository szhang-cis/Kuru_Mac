      SUBROUTINE CCOEFET(CCOEFE,DCCOEFE,
     .                   TEMPG,EFFPL,PROPS)
c    .                   CCOEFFE,CCOEFAE)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE PLASTIC OR VISCOPLASTIC HARDENING
C     COEFFICIENTS (TEMPERATURE DEPENDENT)
C     
C     ISOTT defines the isotropic hardening model considered:
C
C     MODEL 9 Closed form: C = f(\barEp,T), where f is an explicit
C             function of the effective plastic strain and the
C             temperature
C
C     These models can be applied to plastic or viscoplastic analysis
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
C**** INITIALIZATION
C
      CCOEFE=0.0
      DCCOEE=0.0
C
      IPEP4=INT(PROPS(4))
C
      IPEP3=INT(PROPS(3))
      IF(IPEP3.EQ.1) THEN
       IF(ISOTT.EQ.0) RETURN
       IF(ISOTT.GE.1.AND.ISOTT.LE.8) RETURN
       IF(ISOTT.EQ.9) THEN
        IF(EFFPL.LE.VCCOEE(1,2)) THEN
         CCOEFE=VCCOEE(1,1)
         DCCOEFE=0.0D0
        ENDIF
        DO ICCOEE=2,NCCOEE
         I1=ICCOEE-1
         I2=ICCOEE
         IF(EFFPL.GE.VCCOEE(I1,2).AND.EFFPL.LE.VCCOEE(I2,2)) THEN
          CCOEFE=(VCCOEE(I2,1)-VCCOEE(I1,1))/
     .           (VCCOEE(I2,2)-VCCOEE(I1,2))*
     .           (EFFPL-VCCOEE(I1,2))+VCCOEE(I1,1)
          DCCOEFE=(VCCOEE(I2,1)-VCCOEE(I1,1))/
     .            (VCCOEE(I2,2)-VCCOEE(I1,2))
         ENDIF
        ENDDO
        IF(EFFPL.GE.VCCOEE(NCCOEE,2)) THEN
         CCOEFE=VCCOEE(NCCOEE,1)
         DCCOEFE=0.0D0
        ENDIF
        DCCOEFES=(VCCOEE(NCCOEE,1)-VCCOEE(1,1))/       ! secant
     .           (VCCOEE(NCCOEE,2)-VCCOEE(1,2))
c       IF(DCCOEFE.LT.DCCOEFES) DCCOEFE=DCCOEFES       ! possible crit.
        IF(DCCOEFE.eq.0.0d0) DCCOEFE=DCCOEFES          ! possible crit.
C
C**** TEMPERATURE DERIVATIVE OF CCOEF
C
        DCCOEE=0.0D0
       ENDIF
       RETURN
      ENDIF
C
c     IF(IPEP4.EQ.0.OR.IPEP4.EQ.2) THEN
c      IF(ISOTT.EQ.0) RETURN
c      IF(ISOTT.EQ.0) RETURN
c      IF(ISOTT.GE.1.AND.ISOTT.LE.8) RETURN
c      IF(ISOTT.EQ.9) THEN
c
c       IF(TEMPG.LE.VCCOE(1,2)) CCOEF=VCCOE(1,1)
c        DO ICCOE=2,NCCOE
c        I1=ICCOE-1
c        I2=ICCOE
c        IF(TEMPG.GE.VCCOE(I1,2).AND.TEMPG.LE.VCCOE(I2,2)) 
c    .    CCOEF=(VCCOE(I2,1)-VCCOE(I1,1))/(VCCOE(I2,2)-VCCOE(I1,2))*
c    .          (TEMPG-VCCOE(I1,2))+VCCOE(I1,1)
c       ENDDO
c       IF(TEMPG.GE.VCCOE(NCCOE,2)) CCOEF=VCCOE(NCCOE,1)
C
C**** TEMPERATURE DERIVATIVE OF CCOEF
C
c       DCCOE=0.0
c      ENDIF
c      RETURN
c     ENDIF             ! ipep4.eq.0
C
c     IF(IPEP4.EQ.1) THEN
c      IF(ISOTT.EQ.0) RETURN
c      IF(ISOTT.EQ.0) RETURN
c      IF(ISOTT.GE.1.AND.ISOTT.LE.8) RETURN
c      IF(ISOTT.EQ.9) THEN
C
c       IF(TEMPG.LE.VCCOEF(1,2)) CCOEFF=VCCOEF(1,1)
c       DO ICCOE=2,NCCOEF
c        I1=ICCOE-1
c        I2=ICCOE
c        IF(TEMPG.GE.VCCOEF(I1,2).AND.TEMPG.LE.VCCOEF(I2,2))
c    .   CCOEFF=(VCCOEF(I2,1)-VCCOEF(I1,1))/(VCCOEF(I2,2)-VCCOEF(I1,2))*
c    .          (TEMPG-VCCOEF(I1,2))+VCCOEF(I1,1)
c       ENDDO
c       IF(TEMPG.GE.VCCOEF(NCCOEF,2)) CCOEFF=VCCOEF(NCCOEF,1)
C
C**** TEMPERATURE DERIVATIVE OF CCOEFF
C
c       DCCOEF=0.0
C
c       IF(TEMPG.LE.VCCOEA(1,2)) CCOEFA=VCCOEA(1,1)
c       DO ICCOE=2,NCCOEA
c        I1=ICCOE-1
c        I2=ICCOE
c        IF(TEMPG.GE.VCCOEA(I1,2).AND.TEMPG.LE.VCCOEA(I2,2))
c    .   CCOEFA=(VCCOEA(I2,1)-VCCOEA(I1,1))/(VCCOEA(I2,2)-VCCOEA(I1,2))*
c    .          (TEMPG-VCCOEA(I1,2))+VCCOEA(I1,1)
c       ENDDO
c       IF(TEMPG.GE.VCCOEA(NCCOEA,2)) CCOEFA=VCCOEA(NCCOEA,1)
C
C**** TEMPERATURE DERIVATIVE OF CCOEFA
C
c       DCCOEA=0.0
c      ENDIF
c      RETURN
c     ENDIF             ! ipep4.eq.1
C
      END
