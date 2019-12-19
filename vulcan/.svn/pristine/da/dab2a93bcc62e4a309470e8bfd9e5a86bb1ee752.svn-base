      SUBROUTINE VISCOT(VISCO,TEMPG,PROPS,
     .                  VISCOF,VISCOA)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE VISCOSITY COEFFICIENT (TEMPERATURE
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
      IPEP4=INT(PROPS(4))
C
      IPEP3=INT(PROPS(3))
      IF(IPEP3.EQ.1) THEN
       VISCO=VVISC(1,1)
       RETURN
      ENDIF
C
#ifndef restricted
      IF(IPEP4.EQ.0.OR.IPEP4.EQ.2) THEN    ! standard or green sand
       IF(TEMPG.LE.VVISC(1,2)) VISCO=VVISC(1,1)
       DO IVISC=2,NVISC
        I1=IVISC-1
        I2=IVISC
        IF(TEMPG.GE.VVISC(I1,2).AND.TEMPG.LE.VVISC(I2,2)) 
     .    VISCO=(VVISC(I2,1)-VVISC(I1,1))/(VVISC(I2,2)-VVISC(I1,2))*
     .          (TEMPG-VVISC(I1,2))+VVISC(I1,1)
       ENDDO
       IF(TEMPG.GE.VVISC(NVISC,2)) VISCO=VVISC(NVISC,1)
C
C**** TEMPERATURE DERIVATIVE OF VISCO
C
       DVISC=0.0
      ENDIF             ! ipep4.eq.0
C
      IF(IPEP4.EQ.1.OR.IPEP4.EQ.3) THEN    ! non-standard: SGCI or DP
       IF(TEMPG.LE.VVISCF(1,2)) VISCOF=VVISCF(1,1)
       DO IVISC=2,NVISCF
        I1=IVISC-1
        I2=IVISC
        IF(TEMPG.GE.VVISCF(I1,2).AND.TEMPG.LE.VVISCF(I2,2))
     .   VISCOF=(VVISCF(I2,1)-VVISCF(I1,1))/(VVISCF(I2,2)-VVISCF(I1,2))*
     .          (TEMPG-VVISCF(I1,2))+VVISCF(I1,1)
       ENDDO
       IF(TEMPG.GE.VVISCF(NVISCF,2)) VISCOF=VVISCF(NVISCF,1)
C
C**** TEMPERATURE DERIVATIVE OF VISCOF
C
       DVISCF=0.0
C
       IF(TEMPG.LE.VVISCA(1,2)) VISCOA=VVISCA(1,1)
       DO IVISC=2,NVISCA
        I1=IVISC-1
        I2=IVISC
        IF(TEMPG.GE.VVISCA(I1,2).AND.TEMPG.LE.VVISCA(I2,2))
     .   VISCOA=(VVISCA(I2,1)-VVISCA(I1,1))/(VVISCA(I2,2)-VVISCA(I1,2))*
     .          (TEMPG-VVISCA(I1,2))+VVISCA(I1,1)
       ENDDO
       IF(TEMPG.GE.VVISCA(NVISCA,2)) VISCOA=VVISCA(NVISCA,1)
C
C**** TEMPERATURE DERIVATIVE OF VISCOA
C
       DVISCA=0.0
      ENDIF             ! ipep4.eq.1
C
#endif
      RETURN
      END
