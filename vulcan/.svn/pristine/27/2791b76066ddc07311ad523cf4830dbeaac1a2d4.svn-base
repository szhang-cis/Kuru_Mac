      SUBROUTINE CKOEFT(CKOEF,CKOEB,CKOEQ,
     .                  TEMPG,PROPS,
     .                  CKOEFF,CKOEFA,CKOEBF,CKOEBA,CKOEQF,CKOEQA)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE PLASTIC OR VISCOPLASTIC KINEMATIC
C     HARDENING COEFFICIENTS (TEMPERATURE DEPENDENT)
C     
C     IKINE defines the kinematic hardening model considered:
C
C     to be revised (see plvari.f) !!!!!!!!
C
C     MODEL 1 (work hardening; i.e., thesis):
C             CKOEF is the hardening coefficient (h_Cp in thesis)
C             Evolution equation: \dot Cp = \dot lambda h_Cp \sigma : R
C
C     MODEL 2 (linear strain hardening):
C             CKOEF is the hardening coefficient (h_Cp in thesis)
C             Evolution equation: \dot Cp = h_Cp dot \barEp where
C             \barEp=effec. pl. str.
C
C     MODEL 3 (V2 version):
C             CKOEB & CKOEQ: b & Q parameters of V2 viscoplastic model
C             Closed form: Cp = Q (1-exp(-b \barEp)) where
C             \barEp=effec. pl. str.
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
      CKOEF=0.0
      DKCOE=0.0
C
      CKOEB=0.0
      DKCOB=0.0
      CKOEQ=0.0
      DKCOQ=0.0
C
      IPEP4=INT(PROPS(4))
C
      IPEP3=INT(PROPS(3))
      IF(IPEP3.EQ.1) THEN                  ! non-temperature dep. models
       IF(IKINE.EQ.0) RETURN
       GO TO (21,21,22,22) IKINE
   21  CKOEF=VKCOE(1,1)
       RETURN
   22  CKOEB=VKCOB(1,1)
       CKOEQ=VKCOQ(1,1)
       RETURN
      ENDIF
C
#ifndef restricted
      IF(IPEP4.EQ.0.OR.IPEP4.EQ.2) THEN    ! standard or green sand
       IF(IKINE.EQ.0) RETURN
       GO TO (1,1,2,2) IKINE
C
    1  IF(TEMPG.LE.VKCOE(1,2)) CKOEF=VKCOE(1,1)
       DO IKCOE=2,NKCOE
        I1=IKCOE-1
        I2=IKCOE
        IF(TEMPG.GE.VKCOE(I1,2).AND.TEMPG.LE.VKCOE(I2,2)) 
     .   CKOEF=(VKCOE(I2,1)-VKCOE(I1,1))/(VKCOE(I2,2)-VKCOE(I1,2))*
     .         (TEMPG-VKCOE(I1,2))+VKCOE(I1,1)
       ENDDO
       IF(TEMPG.GE.VKCOE(NKCOE,2)) CKOEF=VKCOE(NKCOE,1)
C
C**** TEMPERATURE DERIVATIVE OF CKOEF
C
       DKCOE=0.0
       RETURN
C
    2  IF(TEMPG.LE.VKCOB(1,2)) CKOEB=VKCOB(1,1)
       DO IKCOB=2,NKCOB
        I1=IKCOB-1
        I2=IKCOB
        IF(TEMPG.GE.VKCOB(I1,2).AND.TEMPG.LE.VKCOB(I2,2)) 
     .   CKOEB=(VKCOB(I2,1)-VKCOB(I1,1))/(VKCOB(I2,2)-VKCOB(I1,2))*
     .         (TEMPG-VKCOB(I1,2))+VKCOB(I1,1)
       ENDDO
       IF(TEMPG.GE.VKCOB(NKCOB,2)) CKOEB=VKCOB(NKCOB,1)
C
       IF(TEMPG.LE.VKCOQ(1,2)) CKOEQ=VKCOQ(1,1)
       DO IKCOQ=2,NKCOQ
        I1=IKCOQ-1
        I2=IKCOQ
        IF(TEMPG.GE.VKCOQ(I1,2).AND.TEMPG.LE.VKCOQ(I2,2)) 
     .   CKOEQ=(VKCOQ(I2,1)-VKCOQ(I1,1))/(VKCOQ(I2,2)-VKCOQ(I1,2))*
     .         (TEMPG-VKCOQ(I1,2))+VKCOQ(I1,1)
       ENDDO
       IF(TEMPG.GE.VKCOQ(NKCOQ,2)) CKOEQ=VKCOQ(NKCOQ,1)
C
C**** TEMPERATURE DERIVATIVE OF CKOEB & CKOEQ
C
       DKCOB=0.0
       DKCOQ=0.0
       RETURN
      ENDIF             ! ipep4.eq.0
C
      IF(IPEP4.EQ.1.OR.IPEP4.EQ.3) THEN    ! non-standard: SGCI or DP
       IF(IKINE.EQ.0) RETURN
       GO TO (11,11,12,12) IKINE
C
   11  IF(TEMPG.LE.VKCOEF(1,2)) CKOEFF=VKCOEF(1,1)
       DO IKCOE=2,NKCOEF
        I1=IKCOE-1
        I2=IKCOE
        IF(TEMPG.GE.VKCOEF(I1,2).AND.TEMPG.LE.VKCOEF(I2,2))
     .   CKOEFF=(VKCOEF(I2,1)-VKCOEF(I1,1))/(VKCOEF(I2,2)-VKCOEF(I1,2))*
     .          (TEMPG-VKCOEF(I1,2))+VKCOEF(I1,1)
       ENDDO
       IF(TEMPG.GE.VKCOEF(NKCOEF,2)) CKOEFF=VKCOEF(NKCOEF,1)
C
C**** TEMPERATURE DERIVATIVE OF CKOEFF
C
       DKCOEF=0.0
C
       IF(TEMPG.LE.VKCOEA(1,2)) CKOEFA=VKCOEA(1,1)
       DO IKCOE=2,NKCOEA
        I1=IKCOE-1
        I2=IKCOE
        IF(TEMPG.GE.VKCOEA(I1,2).AND.TEMPG.LE.VKCOEA(I2,2))
     .   CKOEFA=(VKCOEA(I2,1)-VKCOEA(I1,1))/(VKCOEA(I2,2)-VKCOEA(I1,2))*
     .          (TEMPG-VKCOEA(I1,2))+VKCOEA(I1,1)
       ENDDO
       IF(TEMPG.GE.VKCOEA(NKCOEA,2)) CKOEFA=VKCOEA(NKCOEA,1)
C
C**** TEMPERATURE DERIVATIVE OF CKOEFA
C
       DKCOEA=0.0
       RETURN
C
   12  IF(TEMPG.LE.VKCOBF(1,2)) CKOEBF=VKCOBF(1,1)
       DO IKCOB=2,NKCOBF
        I1=IKCOB-1
        I2=IKCOB
        IF(TEMPG.GE.VKCOBF(I1,2).AND.TEMPG.LE.VKCOBF(I2,2))
     .   CKOEBF=(VKCOBF(I2,1)-VKCOBF(I1,1))/(VKCOBF(I2,2)-VKCOBF(I1,2))*
     .          (TEMPG-VKCOBF(I1,2))+VKCOBF(I1,1)
       ENDDO
       IF(TEMPG.GE.VKCOBF(NKCOBF,2)) CKOEBF=VKCOBF(NKCOBF,1)
C
       IF(TEMPG.LE.VKCOQF(1,2)) CKOEQF=VKCOQF(1,1)
       DO IKCOQ=2,NKCOQF
        I1=IKCOQ-1
        I2=IKCOQ
        IF(TEMPG.GE.VKCOQF(I1,2).AND.TEMPG.LE.VKCOQF(I2,2))
     .   CKOEQF=(VKCOQF(I2,1)-VKCOQF(I1,1))/(VKCOQF(I2,2)-VKCOQF(I1,2))*
     .          (TEMPG-VKCOQF(I1,2))+VKCOQF(I1,1)
       ENDDO
       IF(TEMPG.GE.VKCOQF(NKCOQF,2)) CKOEQF=VKCOQF(NKCOQF,1)
C
C**** TEMPERATURE DERIVATIVE OF CKOEBF & CKOEQF
C
       DKCOBF=0.0
       DKCOQF=0.0
C
       IF(TEMPG.LE.VKCOBA(1,2)) CKOEBA=VKCOBA(1,1)
       DO IKCOB=2,NKCOBA
        I1=IKCOB-1
        I2=IKCOB
        IF(TEMPG.GE.VKCOBA(I1,2).AND.TEMPG.LE.VKCOBA(I2,2))
     .   CKOEBA=(VKCOBA(I2,1)-VKCOBA(I1,1))/(VKCOBA(I2,2)-VKCOBA(I1,2))*
     .          (TEMPG-VKCOBA(I1,2))+VKCOBA(I1,1)
       ENDDO
       IF(TEMPG.GE.VKCOBA(NKCOBA,2)) CKOEBA=VKCOBA(NKCOBA,1)
C
       IF(TEMPG.LE.VKCOQA(1,2)) CKOEQA=VKCOQA(1,1)
       DO IKCOQ=2,NKCOQA
        I1=IKCOQ-1
        I2=IKCOQ
        IF(TEMPG.GE.VKCOQA(I1,2).AND.TEMPG.LE.VKCOQA(I2,2))
     .   CKOEQA=(VKCOQA(I2,1)-VKCOQA(I1,1))/(VKCOQA(I2,2)-VKCOQA(I1,2))*
     .          (TEMPG-VKCOQA(I1,2))+VKCOQA(I1,1)
       ENDDO
       IF(TEMPG.GE.VKCOQA(NKCOQA,2)) CKOEQA=VKCOQA(NKCOQA,1)
C
C**** TEMPERATURE DERIVATIVE OF CKOEBA & CKOEQA
C
       DKCOBA=0.0
       DKCOQA=0.0
       RETURN
      ENDIF             ! ipep4.eq.1
C
#endif
      RETURN
      END
