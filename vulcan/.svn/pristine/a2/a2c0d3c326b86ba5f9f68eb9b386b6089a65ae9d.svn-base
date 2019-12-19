      SUBROUTINE CCOEFT(CCOEF,CCOEB,CCOEQ,
     .                  TEMPG,PROPS,
     .                  CCOEFF,CCOEFA,CCOEBF,CCOEBA,CCOEQF,CCOEQA)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE PLASTIC OR VISCOPLASTIC HARDENING
C     COEFFICIENTS (TEMPERATURE DEPENDENT)
C     
C     ISOTT defines the isotropic hardening model considered:
C
C     MODEL 1 (work hardening; i.e., thesis):
C             CCOEF is the hardening coefficient (h_Cp in thesis)
C             Evolution equation: \dot Cp = \dot lambda h_Cp \sigma : R
C
C     MODEL 2 (linear strain hardening):
C             CCOEF is the hardening coefficient (h_Cp in thesis)
C             Evolution equation: \dot Cp = h_Cp dot \barEp where
C             \barEp=effec. pl. str.
C
C     MODEL 3 (V2 version):
C             CCOEB & CCOEQ: b & Q parameters of V2 viscoplastic model
C             Closed form: Cp = Q (1-exp(-b \barEp)) where
C             \barEp=effec. pl. str.
C
C     MODEL 4 (linear version of MODEL 3)
C             CCOEF is the hardening coefficient (secant)
C             Closed form: Cp = h_Cp \barEp
C
C     MODEL 5 (linear strain hardening):
C             CCOEF is the hardening coefficient (h_Cp in thesis)
C             Evolution equation: \dot Cp = h_Cp dot \barEp where
C             \barEp=effec. pl. str. computed without factor 2/3
C
C     MODEL 6 (non-linear strain hardening; non-linear version of
C             MODEL 2):
C             CCOEB & CCOEQ: n & Q parameters
C             Closed form: Cp = A \barEp^n where \barEp=effec. pl. str.
C
C     MODEL 7 (non-linear strain hardening; similar to MODEL=6):
C             CCOEB & CCOEQ: n & Q parameters
C             Closed form: Cp = A (\barEp_o+\barEp)^n
C             where \barEp=effec. pl. str.
C
C     MODEL 8 (model 3 + model 4):
C             CCOEF, CCOEB & CCOEQ: h_Cp, b & Q parameters
C             Closed form: Cp = h_Cp \barEp + Q (1-exp(-b \barEp)) where
C             \barEp=effec. pl. str.
C
C     MODEL 9 (see ccoefet.f)
C
C     MODEL 10 Johnson & Cook hardening model (strain, strain rate &
C              temperature-dependent)
C              CCOEF, CCOEB & CCOEQ: H, n & A parameters
C              Closed form: Cp = C = (Cth + A \barEp^n)*
C                                                  (1 + H ln \dot\barEp)
C              where \barEp=effec. pl. str.
C                    \dot\barEp=effec. pl. str. rate
C              The temperature-dependency is included in Cth & A
C              H & n are normally temperature-independent
C
C     MODEL 11 Idem MODEL 10 with the strain hardening term of MODEL 7
C
C     MODEL 12 Idem MODEL 7 with the consideration of a critical
C              effective plastic strain in order to model the effect
C              of plastification without hardening which is, e.g.,
C              a typical phenomenon of low carbon steels (Luders band
C              formation)
C              CCOEF, CCOEB & CCOEQ: \barEp_c, n & A parameters
C              Closed form: Cp = C = A (\barEp_o+<\barEp-\barEp_c>)^n
C              where \barEp=effec. pl. str., \barEp_c=critical effec.
C              pl. str. and A*\barEp_o^n=Cth
C
C     MODEL 13 Idem MODEL 6 with the consideration of a critical
C              effective plastic strain in order to model the effect
C              of plastification without hardening (Idem MODEL 12)
C              CCOEF, CCOEB & CCOEQ: \barEp_c, n & A parameters
C              Closed form: Cp = A (<\barEp-\barEp_c>)^n
C              where \barEp=effec. pl. str., \barEp_c=critical effec.
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
      CCOEF=0.0D0
      DCCOE=0.0D0
C
      CCOEB=0.0D0
      DCCOB=0.0D0
      CCOEQ=0.0D0
      DCCOQ=0.0D0
C
      IPEP4=INT(PROPS(4))
C
      IPEP3=INT(PROPS(3))
      IF(IPEP3.EQ.1) THEN
       IF(ISOTT.EQ.0) RETURN
c      GO TO (21,21,22,21,21,22,22) ISOTT
       IF(ISOTT.EQ.1.OR.ISOTT.EQ.2.OR.ISOTT.EQ.4.OR.ISOTT.EQ.5.OR.
     .    ISOTT.EQ.8.OR.ISOTT.EQ.10.OR.ISOTT.EQ.11.OR.ISOTT.EQ.12.OR.
     .    ISOTT.EQ.13) THEN
        CCOEF=VCCOE(1,1)
       ENDIF
       IF(ISOTT.EQ.3.OR.ISOTT.EQ.6.OR.ISOTT.EQ.7.OR.ISOTT.EQ.8.OR.
     .    ISOTT.EQ.10.OR.ISOTT.EQ.11.OR.ISOTT.EQ.12.OR.ISOTT.EQ.13) THEN
        CCOEB=VCCOB(1,1)
        CCOEQ=VCCOQ(1,1)
       ENDIF
       RETURN
      ENDIF
C
#ifndef restricted
      IF(IPEP4.EQ.0.OR.IPEP4.EQ.2) THEN    ! standard or green sand
       IF(ISOTT.EQ.0) RETURN
c      GO TO (1,1,2,1,1,2,2) ISOTT
       IF(ISOTT.EQ.1.OR.ISOTT.EQ.2.OR.ISOTT.EQ.4.OR.ISOTT.EQ.5.OR.
     .    ISOTT.EQ.8.OR.ISOTT.EQ.10.OR.ISOTT.EQ.11.OR.ISOTT.EQ.12.OR.
     .    ISOTT.EQ.13) THEN
C
        IF(TEMPG.LE.VCCOE(1,2)) CCOEF=VCCOE(1,1)
         DO ICCOE=2,NCCOE
         I1=ICCOE-1
         I2=ICCOE
         IF(TEMPG.GE.VCCOE(I1,2).AND.TEMPG.LE.VCCOE(I2,2)) 
     .    CCOEF=(VCCOE(I2,1)-VCCOE(I1,1))/(VCCOE(I2,2)-VCCOE(I1,2))*
     .          (TEMPG-VCCOE(I1,2))+VCCOE(I1,1)
        ENDDO
        IF(TEMPG.GE.VCCOE(NCCOE,2)) CCOEF=VCCOE(NCCOE,1)
C
C**** TEMPERATURE DERIVATIVE OF CCOEF
C
        DCCOE=0.0D0
       ENDIF
C
       IF(ISOTT.EQ.3.OR.ISOTT.EQ.6.OR.ISOTT.EQ.7.OR.ISOTT.EQ.8.OR.
     .    ISOTT.EQ.10.OR.ISOTT.EQ.11.OR.ISOTT.EQ.12.OR.ISOTT.EQ.13) THEN
        IF(TEMPG.LE.VCCOB(1,2)) CCOEB=VCCOB(1,1)
        DO ICCOB=2,NCCOB
         I1=ICCOB-1
         I2=ICCOB
         IF(TEMPG.GE.VCCOB(I1,2).AND.TEMPG.LE.VCCOB(I2,2)) 
     .    CCOEB=(VCCOB(I2,1)-VCCOB(I1,1))/(VCCOB(I2,2)-VCCOB(I1,2))*
     .          (TEMPG-VCCOB(I1,2))+VCCOB(I1,1)
        ENDDO
        IF(TEMPG.GE.VCCOB(NCCOB,2)) CCOEB=VCCOB(NCCOB,1)
C
        IF(TEMPG.LE.VCCOQ(1,2)) CCOEQ=VCCOQ(1,1)
        DO ICCOQ=2,NCCOQ
         I1=ICCOQ-1
         I2=ICCOQ
         IF(TEMPG.GE.VCCOQ(I1,2).AND.TEMPG.LE.VCCOQ(I2,2)) 
     .    CCOEQ=(VCCOQ(I2,1)-VCCOQ(I1,1))/(VCCOQ(I2,2)-VCCOQ(I1,2))*
     .          (TEMPG-VCCOQ(I1,2))+VCCOQ(I1,1)
        ENDDO
        IF(TEMPG.GE.VCCOQ(NCCOQ,2)) CCOEQ=VCCOQ(NCCOQ,1)
C
C**** TEMPERATURE DERIVATIVE OF CCOEB & CCOEQ
C
        DCCOB=0.0D0
        DCCOQ=0.0D0
       ENDIF
       RETURN
      ENDIF             ! ipep4.eq.0
C
      IF(IPEP4.EQ.1.OR.IPEP4.EQ.3) THEN    ! non-standard: SGCI or DP
       IF(ISOTT.EQ.0) RETURN
c      GO TO (11,11,12,11,11,12,12) ISOTT
       IF(ISOTT.EQ.1.OR.ISOTT.EQ.2.OR.ISOTT.EQ.4.OR.ISOTT.EQ.5.OR.
     .    ISOTT.EQ.8.OR.ISOTT.EQ.10.OR.ISOTT.EQ.11.OR.ISOTT.EQ.12.OR.
     .    ISOTT.EQ.13) THEN
C
        IF(TEMPG.LE.VCCOEF(1,2)) CCOEFF=VCCOEF(1,1)
        DO ICCOE=2,NCCOEF
         I1=ICCOE-1
         I2=ICCOE
         IF(TEMPG.GE.VCCOEF(I1,2).AND.TEMPG.LE.VCCOEF(I2,2))
     .   CCOEFF=(VCCOEF(I2,1)-VCCOEF(I1,1))/(VCCOEF(I2,2)-VCCOEF(I1,2))*
     .          (TEMPG-VCCOEF(I1,2))+VCCOEF(I1,1)
        ENDDO
        IF(TEMPG.GE.VCCOEF(NCCOEF,2)) CCOEFF=VCCOEF(NCCOEF,1)
C
C**** TEMPERATURE DERIVATIVE OF CCOEFF
C
        DCCOEF=0.0D0
C
        IF(TEMPG.LE.VCCOEA(1,2)) CCOEFA=VCCOEA(1,1)
        DO ICCOE=2,NCCOEA
         I1=ICCOE-1
         I2=ICCOE
         IF(TEMPG.GE.VCCOEA(I1,2).AND.TEMPG.LE.VCCOEA(I2,2))
     .   CCOEFA=(VCCOEA(I2,1)-VCCOEA(I1,1))/(VCCOEA(I2,2)-VCCOEA(I1,2))*
     .          (TEMPG-VCCOEA(I1,2))+VCCOEA(I1,1)
        ENDDO
        IF(TEMPG.GE.VCCOEA(NCCOEA,2)) CCOEFA=VCCOEA(NCCOEA,1)
C
C**** TEMPERATURE DERIVATIVE OF CCOEFA
C
        DCCOEA=0.0D0
       ENDIF
C
       IF(ISOTT.EQ.3.OR.ISOTT.EQ.6.OR.ISOTT.EQ.7.OR.ISOTT.EQ.8.OR.
     .    ISOTT.EQ.10.OR.ISOTT.EQ.11.OR.ISOTT.EQ.12.OR.ISOTT.EQ.13) THEN
        IF(TEMPG.LE.VCCOBF(1,2)) CCOEBF=VCCOBF(1,1)
        DO ICCOB=2,NCCOBF
         I1=ICCOB-1
         I2=ICCOB
         IF(TEMPG.GE.VCCOBF(I1,2).AND.TEMPG.LE.VCCOBF(I2,2))
     .   CCOEBF=(VCCOBF(I2,1)-VCCOBF(I1,1))/(VCCOBF(I2,2)-VCCOBF(I1,2))*
     .          (TEMPG-VCCOBF(I1,2))+VCCOBF(I1,1)
        ENDDO
        IF(TEMPG.GE.VCCOBF(NCCOBF,2)) CCOEBF=VCCOBF(NCCOBF,1)
C
        IF(TEMPG.LE.VCCOQF(1,2)) CCOEQF=VCCOQF(1,1)
        DO ICCOQ=2,NCCOQF
         I1=ICCOQ-1
         I2=ICCOQ
         IF(TEMPG.GE.VCCOQF(I1,2).AND.TEMPG.LE.VCCOQF(I2,2))
     .   CCOEQF=(VCCOQF(I2,1)-VCCOQF(I1,1))/(VCCOQF(I2,2)-VCCOQF(I1,2))*
     .          (TEMPG-VCCOQF(I1,2))+VCCOQF(I1,1)
        ENDDO
        IF(TEMPG.GE.VCCOQF(NCCOQF,2)) CCOEQF=VCCOQF(NCCOQF,1)
C
C**** TEMPERATURE DERIVATIVE OF CCOEBF & CCOEQF
C
        DCCOBF=0.0D0
        DCCOQF=0.0D0
C
        IF(TEMPG.LE.VCCOBA(1,2)) CCOEBA=VCCOBA(1,1)
        DO ICCOB=2,NCCOBA
         I1=ICCOB-1
         I2=ICCOB
         IF(TEMPG.GE.VCCOBA(I1,2).AND.TEMPG.LE.VCCOBA(I2,2))
     .   CCOEBA=(VCCOBA(I2,1)-VCCOBA(I1,1))/(VCCOBA(I2,2)-VCCOBA(I1,2))*
     .          (TEMPG-VCCOBA(I1,2))+VCCOBA(I1,1)
        ENDDO
        IF(TEMPG.GE.VCCOBA(NCCOBA,2)) CCOEBA=VCCOBA(NCCOBA,1)
C
        IF(TEMPG.LE.VCCOQA(1,2)) CCOEQA=VCCOQA(1,1)
        DO ICCOQ=2,NCCOQA
         I1=ICCOQ-1
         I2=ICCOQ
         IF(TEMPG.GE.VCCOQA(I1,2).AND.TEMPG.LE.VCCOQA(I2,2))
     .   CCOEQA=(VCCOQA(I2,1)-VCCOQA(I1,1))/(VCCOQA(I2,2)-VCCOQA(I1,2))*
     .          (TEMPG-VCCOQA(I1,2))+VCCOQA(I1,1)
        ENDDO
        IF(TEMPG.GE.VCCOQA(NCCOQA,2)) CCOEQA=VCCOQA(NCCOQA,1)
C
C**** TEMPERATURE DERIVATIVE OF CCOEBA & CCOEQA
C
        DCCOBA=0.0D0
        DCCOQA=0.0D0
       ENDIF
       RETURN
      ENDIF             ! ipep4.eq.1
C
#endif
      RETURN
      END
