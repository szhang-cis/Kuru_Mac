      SUBROUTINE YOUNGT(TEMPG,PROPS,
     .                  YOUNG,YOUNGF,YOUNGA)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE YOUNG MODULUS COEFFICIENT (TEMPERATURE
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
       IF(IPEP3.EQ.1) THEN                  ! no temp.-dependent
        YOUNG=VYOUN(1,1)
        RETURN
       ENDIF
C
#ifndef restricted
       IPEP4=INT(PROPS(4))
       IF(IPEP4.EQ.0.OR.IPEP4.EQ.2) THEN    ! standard or green sand
        IF(TEMPG.LE.VYOUN(1,2)) YOUNG=VYOUN(1,1)
        DO IYOUN=2,NYOUN
         I1=IYOUN-1
         I2=IYOUN
         IF(TEMPG.GE.VYOUN(I1,2).AND.TEMPG.LE.VYOUN(I2,2)) 
     .    YOUNG=(VYOUN(I2,1)-VYOUN(I1,1))/(VYOUN(I2,2)-VYOUN(I1,2))*
     .          (TEMPG-VYOUN(I1,2))+VYOUN(I1,1)
        ENDDO
        IF(TEMPG.GE.VYOUN(NYOUN,2)) YOUNG=VYOUN(NYOUN,1)
C
        DYOUN=0.0D0                         ! temperature derivative
       ENDIF                ! ipep4.eq.0
C
       IF(IPEP4.EQ.1.OR.IPEP4.EQ.3) THEN    ! non-standard: SGCI or DP
        IF(TEMPG.LE.VYOUNF(1,2)) YOUNGF=VYOUNF(1,1)
        DO IYOUN=2,NYOUNF
         I1=IYOUN-1
         I2=IYOUN
         IF(TEMPG.GE.VYOUNF(I1,2).AND.TEMPG.LE.VYOUNF(I2,2))
     .    YOUNGF=(VYOUNF(I2,1)-VYOUNF(I1,1))/
     .           (VYOUNF(I2,2)-VYOUNF(I1,2))*
     .           (TEMPG-VYOUNF(I1,2))+VYOUNF(I1,1)
        ENDDO
        IF(TEMPG.GE.VYOUNF(NYOUNF,2)) YOUNGF=VYOUNF(NYOUNF,1)
C
        DYOUNF=0.0D0                        ! temperature derivative
C
        IF(TEMPG.LE.VYOUNA(1,2)) YOUNGA=VYOUNA(1,1)
        DO IYOUN=2,NYOUNA
         I1=IYOUN-1
         I2=IYOUN
         IF(TEMPG.GE.VYOUNA(I1,2).AND.TEMPG.LE.VYOUNA(I2,2))
     .    YOUNGA=(VYOUNA(I2,1)-VYOUNA(I1,1))/
     .           (VYOUNA(I2,2)-VYOUNA(I1,2))*
     .           (TEMPG-VYOUNA(I1,2))+VYOUNA(I1,1)
        ENDDO
        IF(TEMPG.GE.VYOUNA(NYOUNA,2)) YOUNGA=VYOUNA(NYOUNA,1)
C
        DYOUNA=0.0D0                        ! temperature derivative
       ENDIF                ! ipep4.eq.1
      ENDIF                 ! ipep1.eq.100
C
      IF(IPEP1.EQ.200) THEN                 ! orthotropic
       IPEP2=INT(PROPS(2))
       IF(IPEP2.EQ.40) RETURN               ! hyperelastic
C
       IPEP3=INT(PROPS(3))
       IF(IPEP3.EQ.1) THEN                  ! no temp.-dependent
        YOUNG=VYOUN(1,1)
        YOUNGF=VYOUNF(1,1)
        YOUNGA=VYOUNA(1,1)
        RETURN
       ENDIF
C
       IPEP4=INT(PROPS(4))
       IF(IPEP4.EQ.0) THEN                  ! standard
        IF(TEMPG.LE.VYOUN(1,2)) YOUNG=VYOUN(1,1)
        DO IYOUN=2,NYOUN
         I1=IYOUN-1
         I2=IYOUN
         IF(TEMPG.GE.VYOUN(I1,2).AND.TEMPG.LE.VYOUN(I2,2))
     .    YOUNG=(VYOUN(I2,1)-VYOUN(I1,1))/(VYOUN(I2,2)-VYOUN(I1,2))*
     .          (TEMPG-VYOUN(I1,2))+VYOUN(I1,1)
        ENDDO
        IF(TEMPG.GE.VYOUN(NYOUN,2)) YOUNG=VYOUN(NYOUN,1)
C
        DYOUN=0.0D0                         ! temperature derivative
C
        IF(TEMPG.LE.VYOUNF(1,2)) YOUNGF=VYOUNF(1,1)
        DO IYOUN=2,NYOUNF
         I1=IYOUN-1
         I2=IYOUN
         IF(TEMPG.GE.VYOUNF(I1,2).AND.TEMPG.LE.VYOUNF(I2,2))
     .    YOUNGF=(VYOUNF(I2,1)-VYOUNF(I1,1))/
     .           (VYOUNF(I2,2)-VYOUNF(I1,2))*
     .           (TEMPG-VYOUNF(I1,2))+VYOUNF(I1,1)
        ENDDO
        IF(TEMPG.GE.VYOUNF(NYOUNF,2)) YOUNGF=VYOUNF(NYOUNF,1)
C
        DYOUNF=0.0D0                        ! temperature derivative
C
        IF(TEMPG.LE.VYOUNA(1,2)) YOUNGA=VYOUNA(1,1)
        DO IYOUN=2,NYOUNA
         I1=IYOUN-1
         I2=IYOUN
         IF(TEMPG.GE.VYOUNA(I1,2).AND.TEMPG.LE.VYOUNA(I2,2))
     .    YOUNGA=(VYOUNA(I2,1)-VYOUNA(I1,1))/
     .           (VYOUNA(I2,2)-VYOUNA(I1,2))*
     .           (TEMPG-VYOUNA(I1,2))+VYOUNA(I1,1)
        ENDDO
        IF(TEMPG.GE.VYOUNA(NYOUNA,2)) YOUNGA=VYOUNA(NYOUNA,1)
C
        DYOUNA=0.0D0                        ! temperature derivative
       ENDIF                ! ipep4.eq.0
#endif
      ENDIF                 ! ipep1.eq.200
C
      RETURN
      END
