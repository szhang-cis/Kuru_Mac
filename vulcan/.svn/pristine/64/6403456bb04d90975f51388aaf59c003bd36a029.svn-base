      SUBROUTINE CCOBTT(CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,CCOB7,CCOB8,
     .                  CCOB9,CCOB10,CCOB11,CCOB12,CLANK,TEMPG,PROPS,
     .                  CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES SOME COEFFICIENTS (TEMPERATURE-DEPENDENT)
C     
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION PROPS(*)
      DIMENSION CLAN1(2),       CLAN2(2),       CLAN3(2),
     .          CLAN4(2),       CLAN5(2),       CLAN6(2)
C
      NCRIT=INT(PROPS(36))
      NCRIP=INT(PROPS(52))
      NASNA=INT(PROPS(47))
C
      ISOTR=INT(PROPS(1))   ! 0=isotropic; 1=orthotropic
      IPEP1=100
      IF(ISOTR.EQ.1) IPEP1=200
C
      IF(IPEP1.EQ.100) THEN               ! isotropic
       IPEP2=INT(PROPS(2))
       IF(IPEP2.EQ.40) RETURN             ! hyperelastic
C
       IPEP3=INT(PROPS(3))
C
       IF(IPEP3.EQ.1) THEN                ! no temp.-dependent
        IF(NCRIT.EQ.39) THEN
         CCOB1=VCOE1(1,1)
         CCOB2=VCOE2(1,1)
         CCOB3=VCOE3(1,1)
         CCOB4=VCOE4(1,1)
         CCOB5=VCOE5(1,1)
        ENDIF
        IF(NCRIT.EQ.42) THEN
         CLANK=VLANK(1,1)
        ENDIF
        IF(NCRIP.EQ.40) THEN
         CCOB6=VCOE6(1,1)
         CCOB7=VCOE7(1,1)
         CCOB8=VCOE8(1,1)
        ENDIF
        IF(IDAMG.EQ.3.OR.IDAMG.EQ.4.OR.IDAMG.EQ.5.OR.IDAMG.EQ.6) THEN
         CCOB9=VCOE9(1,1)
         CCOB10=VCOE10(1,1)
         CCOB11=VCOE11(1,1)
        ENDIF
        IF(IDAMG.EQ.6) THEN
         CCOB12=VCOE12(1,1)
        ENDIF
        RETURN
       ENDIF
#ifndef restricted
C
       IPEP4=INT(PROPS(4))
       IF(IPEP4.EQ.0) THEN                ! standard
C
C**** COEFFICIENT 1
C
        IF(NCRIT.EQ.39) THEN
         IF(TEMPG.LE.VCOE1(1,2)) CCOB1=VCOE1(1,1)
         DO ICOB1=2,NCOE1
          I1=ICOB1-1
          I2=ICOB1
          IF(TEMPG.GE.VCOE1(I1,2).AND.TEMPG.LE.VCOE1(I2,2)) 
     .     CCOB1=(VCOE1(I2,1)-VCOE1(I1,1))/(VCOE1(I2,2)-VCOE1(I1,2))*
     .           (TEMPG-VCOE1(I1,2))+VCOE1(I1,1)
         ENDDO
         IF(TEMPG.GE.VCOE1(NCOE1,2)) CCOB1=VCOE1(NCOE1,1)
C
         DCOE1=0.0D0                      ! temperature derivative
        ENDIF
C
C**** COEFFICIENT 2
C
        IF(NCRIT.EQ.39) THEN
         IF(TEMPG.LE.VCOE2(1,2)) CCOB2=VCOE2(1,1)
         DO ICOB2=2,NCOE2
          I1=ICOB2-1
          I2=ICOB2
          IF(TEMPG.GE.VCOE2(I1,2).AND.TEMPG.LE.VCOE2(I2,2)) 
     .     CCOB2=(VCOE2(I2,1)-VCOE2(I1,1))/(VCOE2(I2,2)-VCOE2(I1,2))*
     .           (TEMPG-VCOE2(I1,2))+VCOE2(I1,1)
         ENDDO
         IF(TEMPG.GE.VCOE2(NCOE2,2)) CCOB2=VCOE2(NCOE2,1)
C
         DCOE2=0.0D0                      ! temperature derivative
        ENDIF
C
C**** COEFFICIENT 3
C
        IF(NCRIT.EQ.39) THEN
         IF(TEMPG.LE.VCOE3(1,2)) CCOB3=VCOE3(1,1)
         DO ICOB3=2,NCOE3
          I1=ICOB3-1
          I2=ICOB3
          IF(TEMPG.GE.VCOE3(I1,2).AND.TEMPG.LE.VCOE3(I2,2)) 
     .     CCOB3=(VCOE3(I2,1)-VCOE3(I1,1))/(VCOE3(I2,2)-VCOE3(I1,2))*
     .           (TEMPG-VCOE3(I1,2))+VCOE3(I1,1)
         ENDDO
         IF(TEMPG.GE.VCOE3(NCOE3,2)) CCOB3=VCOE3(NCOE3,1)
C
         DCOE3=0.0D0                      ! temperature derivative
        ENDIF
C
C**** COEFFICIENT 4
C
        IF(NCRIT.EQ.39) THEN
         IF(TEMPG.LE.VCOE4(1,2)) CCOB4=VCOE4(1,1)
         DO ICOB4=2,NCOE4
          I1=ICOB4-1
          I2=ICOB4
          IF(TEMPG.GE.VCOE4(I1,2).AND.TEMPG.LE.VCOE4(I2,2)) 
     .     CCOB4=(VCOE4(I2,1)-VCOE4(I1,1))/(VCOE4(I2,2)-VCOE4(I1,2))*
     .           (TEMPG-VCOE4(I1,2))+VCOE4(I1,1)
         ENDDO
         IF(TEMPG.GE.VCOE4(NCOE4,2)) CCOB4=VCOE4(NCOE4,1)
C
         DCOE4=0.0D0                      ! temperature derivative
        ENDIF
C
C**** COEFFICIENT 5
C
        IF(NCRIT.EQ.39) THEN
         IF(TEMPG.LE.VCOE5(1,2)) CCOB5=VCOE5(1,1)
         DO ICOB5=2,NCOE5
          I1=ICOB5-1
          I2=ICOB5
          IF(TEMPG.GE.VCOE5(I1,2).AND.TEMPG.LE.VCOE5(I2,2)) 
     .     CCOB5=(VCOE5(I2,1)-VCOE5(I1,1))/(VCOE5(I2,2)-VCOE5(I1,2))*
     .           (TEMPG-VCOE5(I1,2))+VCOE5(I1,1)
         ENDDO
         IF(TEMPG.GE.VCOE5(NCOE5,2)) CCOB5=VCOE5(NCOE5,1)
C
         DCOE5=0.0D0                      ! temperature derivative
        ENDIF
C
C**** AVERAGE LANKFORD COEFFICIENT
C
        IF(NCRIT.EQ.42) THEN
         IF(TEMPG.LE.VLANK(1,2)) CLANK=VLANK(1,1)
         DO ILANK=2,NLANK
          I1=ILANK-1
          I2=ILANK
          IF(TEMPG.GE.VLANK(I1,2).AND.TEMPG.LE.VLANK(I2,2))
     .     CLANK=(VLANK(I2,1)-VLANK(I1,1))/(VLANK(I2,2)-VLANK(I1,2))*
     .           (TEMPG-VLANK(I1,2))+VLANK(I1,1)
         ENDDO
         IF(TEMPG.GE.VLANK(NLANK,2)) CLANK=VLANK(NLANK,1)
C
         DLANK=0.0D0                      ! temperature derivative
        ENDIF
C
C**** COEFFICIENT 6
C
        IF(NCRIP.EQ.40) THEN
         IF(TEMPG.LE.VCOE6(1,2)) CCOB6=VCOE6(1,1)
         DO ICOB6=2,NCOE6
          I1=ICOB6-1
          I2=ICOB6
          IF(TEMPG.GE.VCOE6(I1,2).AND.TEMPG.LE.VCOE6(I2,2)) 
     .     CCOB6=(VCOE6(I2,1)-VCOE6(I1,1))/(VCOE6(I2,2)-VCOE6(I1,2))*
     .           (TEMPG-VCOE6(I1,2))+VCOE6(I1,1)
         ENDDO
         IF(TEMPG.GE.VCOE6(NCOE6,2)) CCOB6=VCOE6(NCOE6,1)
C
         DCOE6=0.0D0                      ! temperature derivative
        ENDIF
C
C**** COEFFICIENT 7
C
        IF(NCRIP.EQ.40) THEN
         IF(TEMPG.LE.VCOE7(1,2)) CCOB7=VCOE7(1,1)
         DO ICOB7=2,NCOE7
          I1=ICOB7-1
          I2=ICOB7
          IF(TEMPG.GE.VCOE7(I1,2).AND.TEMPG.LE.VCOE7(I2,2)) 
     .     CCOB7=(VCOE7(I2,1)-VCOE7(I1,1))/(VCOE7(I2,2)-VCOE7(I1,2))*
     .           (TEMPG-VCOE7(I1,2))+VCOE7(I1,1)
         ENDDO
         IF(TEMPG.GE.VCOE7(NCOE7,2)) CCOB7=VCOE7(NCOE7,1)
C
         DCOE7=0.0D0                      ! temperature derivative
        ENDIF
C
C**** COEFFICIENT 8
C
        IF(NCRIP.EQ.40) THEN
         IF(TEMPG.LE.VCOE8(1,2)) CCOB8=VCOE8(1,1)
         DO ICOB8=2,NCOE8
          I1=ICOB8-1
          I2=ICOB8
          IF(TEMPG.GE.VCOE8(I1,2).AND.TEMPG.LE.VCOE8(I2,2)) 
     .     CCOB8=(VCOE8(I2,1)-VCOE8(I1,1))/(VCOE8(I2,2)-VCOE8(I1,2))*
     .           (TEMPG-VCOE8(I1,2))+VCOE8(I1,1)
         ENDDO
         IF(TEMPG.GE.VCOE8(NCOE8,2)) CCOB8=VCOE8(NCOE8,1)
C
         DCOE8=0.0D0                      ! temperature derivative
        ENDIF
C
C**** COEFFICIENTS 9, 10 & 11
C
        IF(IDAMG.EQ.3.OR.IDAMG.EQ.4.OR.IDAMG.EQ.5.OR.IDAMG.EQ.6) THEN
         IF(TEMPG.LE.VCOE9(1,2)) CCOB9=VCOE9(1,1)
         DO ICOB9=2,NCOE9
          I1=ICOB9-1
          I2=ICOB9
          IF(TEMPG.GE.VCOE9(I1,2).AND.TEMPG.LE.VCOE9(I2,2)) 
     .     CCOB9=(VCOE9(I2,1)-VCOE9(I1,1))/(VCOE9(I2,2)-VCOE9(I1,2))*
     .           (TEMPG-VCOE9(I1,2))+VCOE9(I1,1)
         ENDDO
         IF(TEMPG.GE.VCOE9(NCOE9,2)) CCOB9=VCOE9(NCOE9,1)
C
         DCOE9=0.0D0                      ! temperature derivative
C
         IF(TEMPG.LE.VCOE10(1,2)) CCOB10=VCOE10(1,1)
         DO ICOB10=2,NCOE10
          I1=ICOB10-1
          I2=ICOB10
          IF(TEMPG.GE.VCOE10(I1,2).AND.TEMPG.LE.VCOE10(I2,2))
     .     CCOB10=(VCOE10(I2,1)-VCOE10(I1,1))/
     .            (VCOE10(I2,2)-VCOE10(I1,2))*
     .           (TEMPG-VCOE10(I1,2))+VCOE10(I1,1)
         ENDDO
         IF(TEMPG.GE.VCOE10(NCOE10,2)) CCOB10=VCOE10(NCOE10,1)
C
         DCOE10=0.0D0                     ! temperature derivative
C
         IF(TEMPG.LE.VCOE11(1,2)) CCOB11=VCOE11(1,1)
         DO ICOB11=2,NCOE11
          I1=ICOB11-1
          I2=ICOB11
          IF(TEMPG.GE.VCOE11(I1,2).AND.TEMPG.LE.VCOE11(I2,2))
     .     CCOB11=(VCOE11(I2,1)-VCOE11(I1,1))/
     .            (VCOE11(I2,2)-VCOE11(I1,2))*
     .           (TEMPG-VCOE11(I1,2))+VCOE11(I1,1)
         ENDDO
         IF(TEMPG.GE.VCOE11(NCOE11,2)) CCOB11=VCOE11(NCOE11,1)
C
         DCOE11=0.0D0                     ! temperature derivative
        ENDIF
C
C**** COEFFICIENT 12
C
        IF(IDAMG.EQ.6) THEN
         IF(TEMPG.LE.VCOE12(1,2)) CCOB12=VCOE12(1,1)
         DO ICOB12=2,NCOE12
          I1=ICOB12-1
          I2=ICOB12
          IF(TEMPG.GE.VCOE12(I1,2).AND.TEMPG.LE.VCOE12(I2,2))
     .     CCOB12=(VCOE12(I2,1)-VCOE12(I1,1))/
     .            (VCOE12(I2,2)-VCOE12(I1,2))*
     .            (TEMPG-VCOE9(I1,2))+VCOE9(I1,1)
         ENDDO
         IF(TEMPG.GE.VCOE12(NCOE12,2)) CCOB12=VCOE12(NCOE12,1)
C
         DCOE12=0.0D0                     ! temperature derivative
        ENDIF
       ENDIF                ! ipep4.eq.0
      ENDIF                 ! ipep1.eq.100
C
      IF(IPEP1.EQ.200) THEN               ! orthotropic
       IPEP2=INT(PROPS(2))
       IF(IPEP2.EQ.40) RETURN             ! hyperelastic
C
       IPEP3=INT(PROPS(3))
       IF(IPEP3.EQ.1) THEN                ! no temp.-dependent

        IF(NCRIT.EQ.42) THEN              ! to be revised & completed !
         CLANK=VLANK(1,1)
        ENDIF
        RETURN

       ENDIF
C
       IPEP4=INT(PROPS(4))
       IF(IPEP4.EQ.0) THEN                ! standard
C
C**** HILL, STRENGTH OR LANKFORD COEFFICIENTS FOR YIELD FUNCTION
C
        IF(NCRIT.EQ.42) THEN
         NC42T=INT(PROPS(45))
         IF(NC42T.EQ.1) THEN              ! GENERAL MODEL
          IF(TEMPG.LE.VLAN1(1,2)) CLAN1(1)=VLAN1(1,1)
          DO ILANK=2,NLAN1(1)
           I1=ILANK-1
           I2=ILANK
           IF(TEMPG.GE.VLAN1(I1,2).AND.TEMPG.LE.VLAN1(I2,2))
     .      CLAN1(1)=(VLAN1(I2,1)-VLAN1(I1,1))/
     .               (VLAN1(I2,2)-VLAN1(I1,2))*
     .               (TEMPG-VLAN1(I1,2))+VLAN1(I1,1)
          ENDDO
          IF(TEMPG.GE.VLAN1(NLAN1(1),2)) CLAN1(1)=VLAN1(NLAN1(1),1)
C
          DLAN1=0.0D0                     ! temperature derivative
C
          IF(TEMPG.LE.VLAN2(1,2)) CLAN2(1)=VLAN2(1,1)
          DO ILANK=2,NLAN2(1)
           I1=ILANK-1
           I2=ILANK
           IF(TEMPG.GE.VLAN2(I1,2).AND.TEMPG.LE.VLAN2(I2,2))
     .      CLAN2(1)=(VLAN2(I2,1)-VLAN2(I1,1))/
     .               (VLAN2(I2,2)-VLAN2(I1,2))*
     .               (TEMPG-VLAN2(I1,2))+VLAN2(I1,1)
          ENDDO
          IF(TEMPG.GE.VLAN2(NLAN2(1),2)) CLAN2(1)=VLAN2(NLAN2(1),1)
C
          DLAN2=0.0D0                     ! temperature derivative
C
          IF(TEMPG.LE.VLAN3(1,2)) CLAN3(1)=VLAN3(1,1)
          DO ILANK=2,NLAN3(1)
           I1=ILANK-1
           I2=ILANK
           IF(TEMPG.GE.VLAN3(I1,2).AND.TEMPG.LE.VLAN3(I2,2))
     .      CLAN3(1)=(VLAN3(I2,1)-VLAN3(I1,1))/
     .               (VLAN3(I2,2)-VLAN3(I1,2))*
     .               (TEMPG-VLAN3(I1,2))+VLAN3(I1,1)
          ENDDO
          IF(TEMPG.GE.VLAN3(NLAN3(1),2)) CLAN3(1)=VLAN3(NLAN3(1),1)
C
          DLAN3=0.0D0                     ! temperature derivative
C
          IF(TEMPG.LE.VLAN4(1,2)) CLAN4(1)=VLAN4(1,1)
          DO ILANK=2,NLAN4(1)
           I1=ILANK-1
           I2=ILANK
           IF(TEMPG.GE.VLAN4(I1,2).AND.TEMPG.LE.VLAN4(I2,2))
     .      CLAN4(1)=(VLAN4(I2,1)-VLAN4(I1,1))/
     .               (VLAN4(I2,2)-VLAN4(I1,2))*
     .               (TEMPG-VLAN4(I1,2))+VLAN4(I1,1)
          ENDDO
          IF(TEMPG.GE.VLAN4(NLAN4(1),2)) CLAN4(1)=VLAN4(NLAN4(1),1)
C
          DLAN4=0.0D0                     ! temperature derivative
C
          IF(NTYPE.EQ.4) THEN
           IF(TEMPG.LE.VLAN5(1,2)) CLAN5(1)=VLAN5(1,1)
           DO ILANK=2,NLAN5(1)
            I1=ILANK-1
            I2=ILANK
            IF(TEMPG.GE.VLAN5(I1,2).AND.TEMPG.LE.VLAN5(I2,2))
     .       CLAN5(1)=(VLAN5(I2,1)-VLAN5(I1,1))/
     .                (VLAN5(I2,2)-VLAN5(I1,2))*
     .                (TEMPG-VLAN5(I1,2))+VLAN5(I1,1)
           ENDDO
           IF(TEMPG.GE.VLAN5(NLAN5(1),2)) CLAN5(1)=VLAN5(NLAN5(1),1)
C
           DLAN5=0.0D0                    ! temperature derivative
C
           IF(TEMPG.LE.VLAN6(1,2)) CLAN6(1)=VLAN6(1,1)
           DO ILANK=2,NLAN6(1)
            I1=ILANK-1
            I2=ILANK
            IF(TEMPG.GE.VLAN6(I1,2).AND.TEMPG.LE.VLAN6(I2,2))
     .       CLAN6(1)=(VLAN6(I2,1)-VLAN6(I1,1))/
     .                (VLAN6(I2,2)-VLAN6(I1,2))*
     .                (TEMPG-VLAN6(I1,2))+VLAN6(I1,1)
           ENDDO
           IF(TEMPG.GE.VLAN6(NLAN6(1),2)) CLAN6(1)=VLAN6(NLAN6(1),1)
C
           DLAN6=0.0D0                    ! temperature derivative
          ENDIF
         ENDIF              ! nc42t.eq.1
         IF(NC42T.EQ.2.OR.                ! PLANE STRESS & S-BASED MODEL
     .      NC42T.EQ.3) THEN              ! PLANE STRESS & R-BASED MODEL
          IF(TEMPG.LE.VLAN1(1,2)) CLAN1(1)=VLAN1(1,1)
          DO ILANK=2,NLAN1(1)
           I1=ILANK-1
           I2=ILANK
           IF(TEMPG.GE.VLAN1(I1,2).AND.TEMPG.LE.VLAN1(I2,2))
     .      CLAN1(1)=(VLAN1(I2,1)-VLAN1(I1,1))/
     .               (VLAN1(I2,2)-VLAN1(I1,2))*
     .               (TEMPG-VLAN1(I1,2))+VLAN1(I1,1)
          ENDDO
          IF(TEMPG.GE.VLAN1(NLAN1(1),2)) CLAN1(1)=VLAN1(NLAN1(1),1)
C
          DLAN1=0.0D0                     ! temperature derivative
C
          IF(NTYPE.EQ.3) THEN
           CLAN2(1)=CLAN1(1)
           CLAN3(1)=CLAN1(1)
          ELSE
           IF(TEMPG.LE.VLAN2(1,2)) CLAN2(1)=VLAN2(1,1)
           DO ILANK=2,NLAN2(1)
            I1=ILANK-1
            I2=ILANK
            IF(TEMPG.GE.VLAN2(I1,2).AND.TEMPG.LE.VLAN2(I2,2))
     .       CLAN2(1)=(VLAN2(I2,1)-VLAN2(I1,1))/
     .                (VLAN2(I2,2)-VLAN2(I1,2))*
     .                (TEMPG-VLAN2(I1,2))+VLAN2(I1,1)
           ENDDO
           IF(TEMPG.GE.VLAN2(NLAN2(1),2)) CLAN2(1)=VLAN2(NLAN2(1),1)
C
           DLAN2=0.0D0                    ! temperature derivative
C
           IF(TEMPG.LE.VLAN3(1,2)) CLAN3(1)=VLAN3(1,1)
           DO ILANK=2,NLAN3(1)
            I1=ILANK-1
            I2=ILANK
            IF(TEMPG.GE.VLAN3(I1,2).AND.TEMPG.LE.VLAN3(I2,2))
     .       CLAN3(1)=(VLAN3(I2,1)-VLAN3(I1,1))/
     .                (VLAN3(I2,2)-VLAN3(I1,2))*
     .                (TEMPG-VLAN3(I1,2))+VLAN3(I1,1)
           ENDDO
           IF(TEMPG.GE.VLAN3(NLAN3(1),2)) CLAN3(1)=VLAN3(NLAN3(1),1)
C
           DLAN3=0.0D0                    ! temperature derivative
          ENDIF
         ENDIF              ! nc42t.eq.2,3
        ENDIF               ! ncrit.eq.42
C
C**** HILL OR LANKFORD COEFFICIENTS FOR FLOW POTENTIAL
C
        IF(NCRIP.EQ.42) THEN
         NC42P=INT(PROPS(46))
         IF(NC42P.EQ.1) THEN              ! GENERAL MODEL
          IF(NASNA.EQ.0) THEN             ! associate plasticity
           CLAN1(2)=CLAN1(1)
           CLAN2(2)=CLAN2(1)
           CLAN3(2)=CLAN3(1)
           CLAN4(2)=CLAN4(1)
           CLAN5(2)=CLAN5(1)
           CLAN6(2)=CLAN6(1)
          ELSE                            ! non-associate plasticity
           IF(TEMPG.LE.VLAN1(1,4)) CLAN1(2)=VLAN1(1,3)
           DO ILANK=2,NLAN1(2)
            I1=ILANK-1
            I2=ILANK
            IF(TEMPG.GE.VLAN1(I1,4).AND.TEMPG.LE.VLAN1(I2,4))
     .       CLAN1(2)=(VLAN1(I2,3)-VLAN1(I1,3))/
     .                (VLAN1(I2,4)-VLAN1(I1,4))*
     .                (TEMPG-VLAN1(I1,4))+VLAN1(I1,3)
           ENDDO
           IF(TEMPG.GE.VLAN1(NLAN1(2),4)) CLAN1(2)=VLAN1(NLAN1(2),3)
C
           DLAN1=0.0D0                    ! temperature derivative
C
           IF(TEMPG.LE.VLAN2(1,4)) CLAN2(2)=VLAN2(1,3)
           DO ILANK=2,NLAN2(2)
            I1=ILANK-1
            I2=ILANK
            IF(TEMPG.GE.VLAN2(I1,4).AND.TEMPG.LE.VLAN2(I2,4))
     .       CLAN2(2)=(VLAN2(I2,3)-VLAN2(I1,3))/
     .                (VLAN2(I2,4)-VLAN2(I1,4))*
     .                (TEMPG-VLAN2(I1,4))+VLAN2(I1,3)
           ENDDO
           IF(TEMPG.GE.VLAN2(NLAN1(2),4)) CLAN2(2)=VLAN2(NLAN2(2),3)
C
           DLAN2=0.0D0                    ! temperature derivative
C
           IF(TEMPG.LE.VLAN3(1,4)) CLAN3(2)=VLAN3(1,3)
           DO ILANK=2,NLAN3(2)
            I1=ILANK-1
            I2=ILANK
            IF(TEMPG.GE.VLAN3(I1,4).AND.TEMPG.LE.VLAN3(I2,4))
     .       CLAN3(2)=(VLAN3(I2,3)-VLAN3(I1,3))/
     .                (VLAN3(I2,4)-VLAN3(I1,4))*
     .                (TEMPG-VLAN3(I1,4))+VLAN3(I1,3)
           ENDDO
           IF(TEMPG.GE.VLAN3(NLAN3(2),4)) CLAN3(2)=VLAN3(NLAN3(2),3)
C
           DLAN3=0.0D0                    ! temperature derivative
C
           IF(TEMPG.LE.VLAN4(1,4)) CLAN4(2)=VLAN4(1,3)
           DO ILANK=2,NLAN4(2)
            I1=ILANK-1
            I2=ILANK
            IF(TEMPG.GE.VLAN4(I1,4).AND.TEMPG.LE.VLAN4(I2,4))
     .       CLAN4(2)=(VLAN4(I2,3)-VLAN4(I1,3))/
     .                (VLAN4(I2,4)-VLAN4(I1,4))*
     .                (TEMPG-VLAN4(I1,4))+VLAN4(I1,3)
           ENDDO
           IF(TEMPG.GE.VLAN4(NLAN4(2),4)) CLAN4(2)=VLAN4(NLAN4(2),3)
C
           DLAN4=0.0D0                    ! temperature derivative
C
           IF(NTYPE.EQ.4) THEN
            IF(TEMPG.LE.VLAN5(1,4)) CLAN5(2)=VLAN5(1,3)
            DO ILANK=2,NLAN5(2)
             I1=ILANK-1
             I2=ILANK
             IF(TEMPG.GE.VLAN5(I1,4).AND.TEMPG.LE.VLAN5(I2,4))
     .        CLAN5(2)=(VLAN5(I2,3)-VLAN5(I1,3))/
     .                 (VLAN5(I2,4)-VLAN5(I1,4))*
     .                 (TEMPG-VLAN5(I1,4))+VLAN5(I1,3)
            ENDDO
            IF(TEMPG.GE.VLAN5(NLAN5(2),4)) CLAN5(2)=VLAN5(NLAN5(2),3)
C
            DLAN5=0.0D0                   ! temperature derivative
C
            IF(TEMPG.LE.VLAN6(1,4)) CLAN6(2)=VLAN6(1,3)
            DO ILANK=2,NLAN6(2)
             I1=ILANK-1
             I2=ILANK
             IF(TEMPG.GE.VLAN6(I1,4).AND.TEMPG.LE.VLAN6(I2,4))
     .        CLAN6(2)=(VLAN6(I2,3)-VLAN6(I1,3))/
     .                 (VLAN6(I2,4)-VLAN6(I1,4))*
     .                 (TEMPG-VLAN6(I1,4))+VLAN6(I1,3)
            ENDDO
            IF(TEMPG.GE.VLAN6(NLAN6(2),4)) CLAN6(2)=VLAN6(NLAN6(2),3)
C
            DLAN6=0.0D0                   ! temperature derivative
           ENDIF
          ENDIF             ! nasna.eq.0
         ENDIF              ! nc42p.eq.1
         IF(NC42P.EQ.2.OR.                ! PLANE STRESS & S-BASED MODEL
     .      NC42P.EQ.3) THEN              ! PLANE STRESS & R-BASED MODEL
          IF(NASNA.EQ.0) THEN             ! associate plasticity
           CLAN1(2)=CLAN1(1)
           CLAN2(2)=CLAN2(1)
           CLAN3(2)=CLAN3(1)
          ELSE                            ! non-associate plasticity
           IF(TEMPG.LE.VLAN1(1,4)) CLAN1(2)=VLAN1(1,3)
           DO ILANK=2,NLAN1(2)
            I1=ILANK-1
            I2=ILANK
            IF(TEMPG.GE.VLAN1(I1,4).AND.TEMPG.LE.VLAN1(I2,4))
     .       CLAN1(2)=(VLAN1(I2,3)-VLAN1(I1,3))/
     .                (VLAN1(I2,4)-VLAN1(I1,4))*
     .                (TEMPG-VLAN1(I1,4))+VLAN1(I1,3)
           ENDDO
           IF(TEMPG.GE.VLAN1(NLAN1(2),4)) CLAN1(2)=VLAN1(NLAN1(2),3)
C
           DLAN1=0.0D0                    ! temperature derivative
C
           IF(NTYPE.EQ.3) THEN
            CLAN2(2)=CLAN1(2)
            CLAN3(2)=CLAN1(2)
           ELSE
            IF(TEMPG.LE.VLAN2(1,4)) CLAN2(2)=VLAN2(1,3)
            DO ILANK=2,NLAN2(2)
             I1=ILANK-1
             I2=ILANK
             IF(TEMPG.GE.VLAN2(I1,4).AND.TEMPG.LE.VLAN2(I2,4))
     .        CLAN2(2)=(VLAN2(I2,3)-VLAN2(I1,3))/
     .                 (VLAN2(I2,4)-VLAN2(I1,4))*
     .                 (TEMPG-VLAN2(I1,4))+VLAN2(I1,3)
            ENDDO
            IF(TEMPG.GE.VLAN2(NLAN2(2),4)) CLAN2(2)=VLAN2(NLAN2(2),3)
C
            DLAN2=0.0D0                   ! temperature derivative
C
            IF(TEMPG.LE.VLAN3(1,4)) CLAN3(2)=VLAN3(1,3)
            DO ILANK=2,NLAN3(2)
             I1=ILANK-1
             I2=ILANK
             IF(TEMPG.GE.VLAN3(I1,4).AND.TEMPG.LE.VLAN3(I2,4))
     .        CLAN3(2)=(VLAN3(I2,3)-VLAN3(I1,3))/
     .                 (VLAN3(I2,4)-VLAN3(I1,4))*
     .                 (TEMPG-VLAN3(I1,4))+VLAN3(I1,3)
            ENDDO
            IF(TEMPG.GE.VLAN3(NLAN3(2),4)) CLAN3(2)=VLAN3(NLAN3(2),3)
C
            DLAN3=0.0D0                   ! temperature derivative
           ENDIF
          ENDIF             ! nasna.eq.0
         ENDIF              ! nc42p.eq.3
        ENDIF               ! ncrip.eq.42
       ENDIF                ! ipep4.eq.0
#endif
      ENDIF                 ! ipep1.eq.200
C
      RETURN
      END
