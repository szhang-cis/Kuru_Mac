      SUBROUTINE EXPONT(EXPON,TEMPG,PROPS,
     .                  EXPONF,EXPONA)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE EXPONENT (TEMPERATURE DEPENDENT)
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
      IPEP4=INT(PROPS(4))
C
      IPEP3=INT(PROPS(3))
      IF(IPEP3.EQ.1) THEN
       EXPON=VEXPO(1,1)
       RETURN
      ENDIF
#ifndef restricted
C
      IF(IPEP4.EQ.0.OR.IPEP4.EQ.2) THEN    ! standard or green sand
       IF(TEMPG.LE.VEXPO(1,2)) EXPON=VEXPO(1,1)
       DO IEXPO=2,NEXPO
        I1=IEXPO-1
        I2=IEXPO
        IF(TEMPG.GE.VEXPO(I1,2).AND.TEMPG.LE.VEXPO(I2,2)) 
     .    EXPON=(VEXPO(I2,1)-VEXPO(I1,1))/(VEXPO(I2,2)-VEXPO(I1,2))*
     .          (TEMPG-VEXPO(I1,2))+VEXPO(I1,1)
       ENDDO
       IF(TEMPG.GE.VEXPO(NEXPO,2)) EXPON=VEXPO(NEXPO,1)
C
C**** TEMPERATURE DERIVATIVE OF EXPON
C
       DEXPO=0.0
      ENDIF               ! ipep4.eq.0
C
      IF(IPEP4.EQ.1.OR.IPEP4.EQ.3) THEN    ! non-standard: SGCI or DP
       IF(TEMPG.LE.VEXPOF(1,2)) EXPONF=VEXPOF(1,1)
       DO IEXPO=2,NEXPOF
        I1=IEXPO-1
        I2=IEXPO
        IF(TEMPG.GE.VEXPOF(I1,2).AND.TEMPG.LE.VEXPOF(I2,2))
     .   EXPONF=(VEXPOF(I2,1)-VEXPOF(I1,1))/(VEXPOF(I2,2)-VEXPOF(I1,2))*
     .          (TEMPG-VEXPOF(I1,2))+VEXPOF(I1,1)
       ENDDO
       IF(TEMPG.GE.VEXPOF(NEXPOF,2)) EXPONF=VEXPOF(NEXPOF,1)
C
C**** TEMPERATURE DERIVATIVE OF EXPONF
C
       DEXPOF=0.0
C
       IF(TEMPG.LE.VEXPOA(1,2)) EXPONA=VEXPOA(1,1)
       DO IEXPO=2,NEXPOA
        I1=IEXPO-1
        I2=IEXPO
        IF(TEMPG.GE.VEXPOA(I1,2).AND.TEMPG.LE.VEXPOA(I2,2))
     .   EXPONA=(VEXPOA(I2,1)-VEXPOA(I1,1))/(VEXPOA(I2,2)-VEXPOA(I1,2))*
     .          (TEMPG-VEXPOA(I1,2))+VEXPOA(I1,1)
       ENDDO
       IF(TEMPG.GE.VEXPOA(NEXPOA,2)) EXPONA=VEXPOA(NEXPOA,1)
C
C**** TEMPERATURE DERIVATIVE OF EXPONA
C
       DEXPOA=0.0
      ENDIF               ! ipep4.eq.1
C
#endif
      RETURN
      END
