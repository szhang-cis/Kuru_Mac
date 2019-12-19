      SUBROUTINE CPOEFT(CPOEF,CPOE1,CPOE2,TEMPG,PROPS)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE POROSITY COEFFICIENTS (TEMP.-DEPENDENT)
C     
C     IPORO defines the porosity model considered:
C
C     MODEL 1 Gurson model: only nucleation
C
C     MODEL 2 Gurson model: only growth
C
C     MODEL 3 Gurson model: nucleation & growth
C
C     MODEL 4 Gurson-Tvergaard model: nucleation & growth
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
      CPOEF=0.0
      DPCOE=0.0
C
      IPEP4=INT(PROPS(4))
C
      IPEP3=INT(PROPS(3))
      IF(IPEP3.EQ.1) THEN
       IF(IPORO.EQ.0) RETURN
       CPOEF=VPCOE(1,1)
       IF(IPORO.EQ.4) THEN
        CPOE1=VPCO1(1,1)
        CPOE2=VPCO2(1,1)
       ENDIF
       RETURN
      ENDIF
#ifndef restricted
C
      IF(IPEP4.EQ.0.OR.IPEP4.EQ.2) THEN
       IF(IPORO.EQ.0) RETURN
C
       IF(TEMPG.LE.VPCOE(1,2)) CPOEF=VPCOE(1,1)
       DO IPCOE=2,NPCOE
        I1=IPCOE-1
        I2=IPCOE
        IF(TEMPG.GE.VPCOE(I1,2).AND.TEMPG.LE.VPCOE(I2,2)) 
     .   CPOEF=(VPCOE(I2,1)-VPCOE(I1,1))/(VPCOE(I2,2)-VPCOE(I1,2))*
     .         (TEMPG-VPCOE(I1,2))+VPCOE(I1,1)
       ENDDO
       IF(TEMPG.GE.VPCOE(NPCOE,2)) CPOEF=VPCOE(NPCOE,1)
C
C**** TEMPERATURE DERIVATIVE OF CPOEF
C
       DPCOE=0.0
C
       IF(IPORO.EQ.4) THEN
        IF(TEMPG.LE.VPCO1(1,2)) CPOE1=VPCO1(1,1)
        DO IPCO1=2,NPCO1
         I1=IPCO1-1
         I2=IPCO1
         IF(TEMPG.GE.VPCO1(I1,2).AND.TEMPG.LE.VPCO1(I2,2))
     .    CPOE1=(VPCO1(I2,1)-VPCO1(I1,1))/(VPCO1(I2,2)-VPCO1(I1,2))*
     .          (TEMPG-VPCO1(I1,2))+VPCO1(I1,1)
        ENDDO
        IF(TEMPG.GE.VPCO1(NPCO1,2)) CPOE1=VPCO1(NPCO1,1)
C
C**** TEMPERATURE DERIVATIVE OF CPOE1
C
        DPCO1=0.0
C
        IF(TEMPG.LE.VPCO2(1,2)) CPOE2=VPCO2(1,1)
        DO IPCO2=2,NPCO2
         I1=IPCO2-1
         I2=IPCO2
         IF(TEMPG.GE.VPCO2(I1,2).AND.TEMPG.LE.VPCO2(I2,2))
     .    CPOE2=(VPCO2(I2,1)-VPCO2(I1,1))/(VPCO2(I2,2)-VPCO2(I1,2))*
     .          (TEMPG-VPCO2(I1,2))+VPCO2(I1,1)
        ENDDO
        IF(TEMPG.GE.VPCO2(NPCO2,2)) CPOE2=VPCO2(NPCO2,1)
C
C**** TEMPERATURE DERIVATIVE OF CPOE2
C
        DPCO2=0.0
C
       ENDIF            ! iporo.eq.4
       RETURN
      ENDIF             ! ipep4.eq.0
C
      IF(IPEP4.EQ.1) THEN
       IF(IPORO.EQ.0) RETURN
C
       IF(TEMPG.LE.VPCOEF(1,2)) CPOEFF=VPCOEF(1,1)
       DO IPCOE=2,NPCOEF
        I1=IPCOE-1
        I2=IPCOE
        IF(TEMPG.GE.VPCOEF(I1,2).AND.TEMPG.LE.VPCOEF(I2,2))
     .   CPOEFF=(VPCOEF(I2,1)-VPCOEF(I1,1))/(VPCOEF(I2,2)-VPCOEF(I1,2))*
     .          (TEMPG-VPCOEF(I1,2))+VPCOEF(I1,1)
       ENDDO
       IF(TEMPG.GE.VPCOEF(NPCOEF,2)) CPOEFF=VPCOEF(NPCOEF,1)
C
C**** TEMPERATURE DERIVATIVE OF CPOEFF
C
       DPCOEF=0.0
C
       IF(TEMPG.LE.VPCOEA(1,2)) CPOEFA=VPCOEA(1,1)
       DO IPCOE=2,NPCOEA
        I1=IPCOE-1
        I2=IPCOE
        IF(TEMPG.GE.VPCOEA(I1,2).AND.TEMPG.LE.VPCOEA(I2,2))
     .   CPOEFA=(VPCOEA(I2,1)-VPCOEA(I1,1))/(VPCOEA(I2,2)-VPCOEA(I1,2))*
     .          (TEMPG-VPCOEA(I1,2))+VPCOEA(I1,1)
       ENDDO
       IF(TEMPG.GE.VPCOEA(NPCOEA,2)) CPOEFA=VPCOEA(NPCOEA,1)
C
C**** TEMPERATURE DERIVATIVE OF CPOEFA
C
       DPCOEA=0.0
C
       IF(IPORO.EQ.4) THEN
        CALL RUNEND('CPOEFT: MODEL 4 NOT IMPLEMENTED')
       ENDIF            ! iporo.eq.4
       RETURN
      ENDIF             ! ipep4.eq.1
C
#endif
      RETURN
      END
