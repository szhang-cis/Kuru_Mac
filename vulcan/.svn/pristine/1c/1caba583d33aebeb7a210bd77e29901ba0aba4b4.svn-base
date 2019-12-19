      SUBROUTINE CCEROT(TEMPG,PROPS,DTEMG,
     .                  CCERO,DCCER,CCEROF,CCEROA)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE ELASTIC LIMIT OR TOTAL HARDENING
C     FUNCTION (TEMPERATURE DEPENDENT FOR COUPLED PROBLEMS)
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
      IPEP3=INT(PROPS(3))
      IF(IPEP3.EQ.1) THEN
       CCERO=VCCER(1,1)
       RETURN
      ENDIF
C
#ifndef restricted
      IPEP4=INT(PROPS(4))
      IF(IPEP4.EQ.0.OR.IPEP4.EQ.2) THEN    ! standard or green sand
       TEMPX=TEMPG
       IF(TEMPX.LE.VCCER(1,2)) CCERX=VCCER(1,1)
       DO ICCER=2,NCCER
        I1=ICCER-1
        I2=ICCER
        IF(TEMPX.GE.VCCER(I1,2).AND.TEMPX.LE.VCCER(I2,2)) 
     .    CCERX=(VCCER(I2,1)-VCCER(I1,1))/(VCCER(I2,2)-VCCER(I1,2))*
     .          (TEMPX-VCCER(I1,2))+VCCER(I1,1)
       ENDDO
       IF(TEMPX.GE.VCCER(NCCER,2)) CCERX=VCCER(NCCER,1)
       CCERO=CCERX
C
C**** TEMPERATURE DERIVATIVE OF CCERO
C
       DCCER=0.0D0
       IF(ITERME.GE.0) THEN     ! coupled problems
        DTEM1=DABS(DTEMG)
        TOLTT=1.0D-06
        IF(DTEM1.GT.TOLTT) THEN
         CCER2=CCERO
         TEMPX=TEMPG-DTEMG
         IF(TEMPX.LE.VCCER(1,2)) CCERX=VCCER(1,1)
         DO ICCER=2,NCCER
          I1=ICCER-1
          I2=ICCER
          IF(TEMPX.GE.VCCER(I1,2).AND.TEMPX.LE.VCCER(I2,2)) 
     .      CCERX=(VCCER(I2,1)-VCCER(I1,1))/(VCCER(I2,2)-VCCER(I1,2))*
     .            (TEMPX-VCCER(I1,2))+VCCER(I1,1)
         ENDDO
         IF(TEMPX.GE.VCCER(NCCER,2)) CCERX=VCCER(NCCER,1)
         CCER1=CCERX
         DCCER=(CCER2-CCER1)/DTEMG
        ENDIF
       ENDIF
      ENDIF                   ! ipep4.eq.0
C
      IF(IPEP4.EQ.1.OR.IPEP4.EQ.3) THEN    ! non-standard: SGCI or DP
       TEMPX=TEMPG
       IF(TEMPX.LE.VCCERF(1,2)) CCERXF=VCCERF(1,1)
       DO ICCER=2,NCCERF
        I1=ICCER-1
        I2=ICCER
        IF(TEMPX.GE.VCCERF(I1,2).AND.TEMPX.LE.VCCERF(I2,2))
     .   CCERXF=(VCCERF(I2,1)-VCCERF(I1,1))/(VCCERF(I2,2)-VCCERF(I1,2))*
     .          (TEMPX-VCCERF(I1,2))+VCCERF(I1,1)
       ENDDO
       IF(TEMPX.GE.VCCERF(NCCERF,2)) CCERXF=VCCERF(NCCERF,1)
       CCEROF=CCERXF
C
       IF(TEMPX.LE.VCCERA(1,2)) CCERXA=VCCERA(1,1)
       DO ICCER=2,NCCERA
        I1=ICCER-1
        I2=ICCER
        IF(TEMPX.GE.VCCERA(I1,2).AND.TEMPX.LE.VCCERA(I2,2))
     .   CCERXA=(VCCERA(I2,1)-VCCERA(I1,1))/(VCCERA(I2,2)-VCCERA(I1,2))*
     .          (TEMPX-VCCERA(I1,2))+VCCERA(I1,1)
       ENDDO
       IF(TEMPX.GE.VCCERA(NCCERA,2)) CCERXA=VCCERA(NCCERA,1)
       CCEROA=CCERXA
      ENDIF                   ! ipep4.eq.1
#endif
C
      RETURN
      END
