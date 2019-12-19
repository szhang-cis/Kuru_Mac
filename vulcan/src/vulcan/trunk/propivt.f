      SUBROUTINE PROPIVT(ITAPE,PROPS)
C***********************************************************************
C
C**** THIS ROUTINE READS THE MATERIAL PROPERTIES FOR THE ISOTROPIC, 
C     VISCOPLASTIC & TEMPERATURE-DEPENDENT SOLID MODEL WITH DAMAGE
C     (ELEMENT NO. 30)
C
C**** OUTPUT (see pointe.f)
C
C     KPLA1: isotropic hardening
C     KPLA2: kinematic hardening
C     KPLA3: damage
C     KPLA4: porosity (Gurson model)
C     KPLA5: shrinkage
C     KPLA6: fatigue
C     KPLA7: strain-rate-dependent hardening models
C     KPLA8: recrystallization
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'auxl_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'inpo_om.f'
C
      DIMENSION PROPS(*)
C
C**** ESTABLISH SOME PARAMETERS
C
      PROPS( 1)= 0.0D0          ! isotropic
      PROPS( 2)=30.0D0          ! viscoplastic
      PROPS( 3)= 2.0D0          ! temperature-dependent
      PROPS(59)= 3.0D0          ! solidification
C
      PROPS(37)= 0.0D0          ! non prescribed strains (see incstn.f)
      PROPS(19)= 0.0D0
C
      PROPS(21)= 0.0D0          ! no external damping
      PROPS(22)= 0.0D0
C
C**** READ & WRITE MATERIAL PROPERTIES
C
      NPRIN=0
C
C**** LOOK FOR 'DENSITY' CARD
C
      CALL LISTEN('PROPIVT',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'DENSI') THEN
       NLINE=INT(PARAM(1))
       IF(NLINE.GT.20)
     .  CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN DENSITY')
       NPRIN=1
       DO LINEA=1,NLINE
        CALL LISTEN('PROPIVT',NPRIN,ITAPE)
        PROPS(12)=PARAM(1)
       ENDDO
       WRITE(LURES,801) PROPS(12)
       WRITE(LURES,899)
      ELSE
       IF(KDYNA.EQ.1) THEN
        CALL RUNEND('PROPIVT: DENSITY CARD NOT FOUND')
       ELSE
        PROPS(12)=0.0D0
        CALL RUNMEN('WARNING: ZERO DENSITY VALUE IS ASSUMED')
        GO TO 1
       ENDIF
      ENDIF
C
C**** LOOK FOR 'YOUNG_MODULUS' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIVT',NPRIN,ITAPE)
    1 IF(WORDS(1).EQ.'YOUNG') THEN
       NPOI1=61
       NPOI2=NPOI1
       PROPS(NPOI1)=PARAM(1)
       NLINE=INT(PARAM(1))
       IF(NLINE.GT.20)
     .  CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN YOUNG MODULUS')
       NPRIN=1
       DO LINEA=1,NLINE
        CALL LISTEN('PROPIVT',NPRIN,ITAPE)
        DO I=1,2
         NPOI2=NPOI2+1
         PROPS(NPOI2)=PARAM(I)
        ENDDO
       ENDDO
       WRITE(LURES,803)
       DO IPROP=NPOI1+1,NPOI2,2
        WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
       ENDDO
       WRITE(LURES,899)
      ELSE
       CALL RUNEND('PROPIVT: YOUNG_MODULUS CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'POISSON_RATIO' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIVT',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'POISS') THEN
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI1)=PARAM(1)
       NLINE=INT(PARAM(1))
       IF(NLINE.GT.20)
     .  CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN POISSON_RATIO')
       NPRIN=1
       DO LINEA=1,NLINE
        CALL LISTEN('PROPIVT',NPRIN,ITAPE)
        DO I=1,2
         NPOI2=NPOI2+1
         PROPS(NPOI2)=PARAM(I)
        ENDDO
       ENDDO
       WRITE(LURES,805)
       DO IPROP=NPOI1+1,NPOI2,2
        WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
       ENDDO
       WRITE(LURES,899)
      ELSE
       CALL RUNEND('PROPIVT: POISSON_RATIO CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'REFERENCE_TEMPERATURE' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIVT',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'REFER') THEN
       PROPS(11)=PARAM(1)
      ELSE
       CALL RUNEND('PROPIVT: REFERENCE TEMPERATURE CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'THERMAL_DILATATION' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIVT',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'THERM') THEN
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI1)=PARAM(1)
       NLINE=INT(PARAM(1))
       IF(NLINE.GT.20)
     .  CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN THERMAL DIL.')
C
       NPOI2=NPOI2+1
       PROPS(NPOI2)=0.0D0                          ! nalpu
       IF(WORDS(2).EQ.'TANGE') PROPS(NPOI2)=1.0D0  ! nalpu
C
       NPRIN=1
       DO LINEA=1,NLINE
        CALL LISTEN('PROPIVT',NPRIN,ITAPE)
        DO I=1,2
         NPOI2=NPOI2+1
         PROPS(NPOI2)=PARAM(I)
        ENDDO
       ENDDO
       WRITE(LURES,806)
       DO IPROP=NPOI1+2,NPOI2,2
        WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
       ENDDO
       WRITE(LURES,899)
      ELSE
       CALL RUNEND('PROPIVT: THERMAL_DILAT. CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'TH_HD_FUNCTION' CARD (THERMAL HARDENING)
C
      NPRIN=0
      CALL LISTEN('PROPIVT',NPRIN,ITAPE)
      IF((WORDS(1).EQ.'HD_FU').OR.(WORDS(1).EQ.'TH_HD')) THEN
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI1)=PARAM(1)
       NLINE=INT(PARAM(1))
       IF(NLINE.GT.20)
     .  CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN TH. HD. FUNC.')
       NPRIN=1
       DO LINEA=1,NLINE
        CALL LISTEN('PROPIVT',NPRIN,ITAPE)
        DO I=1,2
         NPOI2=NPOI2+1
         PROPS(NPOI2)=PARAM(I)
        ENDDO
       ENDDO
       WRITE(LURES,807)
       DO IPROP=NPOI1+1,NPOI2,2
        WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
       ENDDO
       WRITE(LURES,899)
      ELSE
       CALL RUNEND('PROPIVT: THERMAL HARDEN. CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'FLOW_PARAMETERS' CARD (PARAMETERS OF VISCOUS LAW)
C
      NPRIN=0
      CALL LISTEN('PROPIVT',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'FLOW_') THEN
       IF(WORDS(2).EQ.'MODEL') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)      ! ivifl
        IVIFL=INT(PARAM(1))
       ELSE                        ! default: flow model=1
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=1.0D0
        IVIFL=1
       ENDIF
      ELSE
       CALL RUNEND('PROPIVT: FLOW_PARAMETERS CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'VISCOSITY' CARD OR SMA COEFFICIENT
C
      IF(IVIFL.EQ.1.OR.IVIFL.EQ.2.OR.IVIFL.EQ.3.OR.IVIFL.EQ.6) THEN
       NPRIN=0
       CALL LISTEN('PROPIVT',NPRIN,ITAPE)
       IF(WORDS(1).EQ.'VISCO') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN VISCOSITY')
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         DO I=1,2
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(I)
         ENDDO
        ENDDO
        WRITE(LURES,808)
        DO IPROP=NPOI1+1,NPOI2,2
         WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
        ENDDO
        WRITE(LURES,899)
       ELSE
        CALL RUNEND('PROPIVT: VISCOSITY CARD NOT FOUND')
       ENDIF
      ENDIF           ! ivifl.eq.1.or.ivifl.eq.2.or.ivifl.eq.3.or. ...
C
C**** LOOK FOR ZABARAS' COEFFICIENTS (A, B & C) OR SMA COEFFICIENTS
C
      IF(IVIFL.EQ.4.OR.IVIFL.EQ.6) THEN
C
C**** COEFFICIENT A
C
       NPRIN=0
       CALL LISTEN('PROPIVT',NPRIN,ITAPE)
       IF(WORDS(1).EQ.'A_COE') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN A_COEF')
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         DO I=1,2
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(I)
         ENDDO
        ENDDO
        WRITE(LURES,823)
        DO IPROP=NPOI1+1,NPOI2,2
         WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
        ENDDO
        WRITE(LURES,899)
       ELSE
        CALL RUNEND('PROPIVT: A_COEFFICIENT CARD NOT FOUND')
       ENDIF
C
C**** COEFFICIENT B
C
       NPRIN=0
       CALL LISTEN('PROPIVT',NPRIN,ITAPE)
       IF(WORDS(1).EQ.'B_COE') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN B_COEF')
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         DO I=1,2
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(I)
         ENDDO
        ENDDO
        WRITE(LURES,824)
        DO IPROP=NPOI1+1,NPOI2,2
         WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
        ENDDO
        WRITE(LURES,899)
       ELSE
        CALL RUNEND('PROPIVT: B_COEFFICIENT CARD NOT FOUND')
       ENDIF
C
C**** COEFFICIENT C
C
       NPRIN=0
       CALL LISTEN('PROPIVT',NPRIN,ITAPE)
       IF(WORDS(1).EQ.'C_COE') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN C_COEF')
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         DO I=1,2
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(I)
         ENDDO
        ENDDO
        WRITE(LURES,825)
        DO IPROP=NPOI1+1,NPOI2,2
         WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
        ENDDO
        WRITE(LURES,899)
       ELSE
        CALL RUNEND('PROPIVT: C_COEFFICIENT CARD NOT FOUND')
       ENDIF
C
       CALL RUNMEN('WARNING: ZABARAS LAW USES ZERO THERMAL HD FUNCTION &
     . ZERO HD COEFFICIENTS')
C
      ENDIF           ! ivifl.eq.4.or.ivifl.eq.6
C
C**** LOOK FOR 'EXPONENT' CARD (EXPONENT OF VISCOUS LAW) OR SMA COEF.
C
      IF((IVIFL.GE.1.AND.IVIFL.LE.4).OR.IVIFL.EQ.6) THEN
       NPRIN=0
       CALL LISTEN('PROPIVT',NPRIN,ITAPE)
       IF(WORDS(1).EQ.'EXPON') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN EXPONENT')
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         DO I=1,2
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(I)
         ENDDO
        ENDDO
        WRITE(LURES,809)
        DO IPROP=NPOI1+1,NPOI2,2
         WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
        ENDDO
        WRITE(LURES,899)
       ELSE
        CALL RUNEND('PROPIVT: EXPONENT CARD NOT FOUND')
       ENDIF
      ENDIF           ! (ivifl.ge.1.or.ivifl.le.4).or. ...
C
C**** LOOK FOR BUSSO'S COEFFICIENTS
C
      IF(IVIFL.EQ.5) THEN
       NPRIN=1
       CALL LISTEN('PROPIVT',NPRIN,ITAPE)
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI2)=PARAM(1)          ! Eo_dot
       NPOI2=NPOI2+1
       PROPS(NPOI2)=PARAM(2)          ! w
       NPOI2=NPOI2+1
       PROPS(NPOI2)=PARAM(3)          ! Qac
       NPOI2=NPOI2+1
       PROPS(NPOI2)=PARAM(4)          ! R (universal gas constant)
       NPOI2=NPOI2+1
       PROPS(NPOI2)=PARAM(5)          ! chi
       NPOI2=NPOI2+1
       PROPS(NPOI2)=PARAM(6)          ! n1
       WRITE(LURES,851) PROPS(NPOI1)
       WRITE(LURES,852) PROPS(NPOI1+1)
       WRITE(LURES,853) PROPS(NPOI1+2)
       WRITE(LURES,854) PROPS(NPOI1+3)
       WRITE(LURES,855) PROPS(NPOI1+4)
       WRITE(LURES,856) PROPS(NPOI1+5)
      ENDIF           ! ivifl.eq.5
C
      IF(IVIFL.EQ.6) THEN
       KAUXX=1
       IF(KAUXX.GT.KPLA9) KPLA9=1
      ENDIF           ! ivifl.eq.6
C
C**** LOOK FOR 'ISOTROPIC_HARDENING' CARD
C
      NPRIN=0
      ISOTT=0
      CALL LISTEN('PROPIVT',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'ISOTR') THEN
       KAUXX=1
       IF(KAUXX.GT.KPLA1) KPLA1=1
       IF(WORDS(2).EQ.'MODEL') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        IIAUX=INT(PARAM(1))
        ISOTT=IIAUX
        IF((IIAUX.LE.0.OR.IIAUX.GT.8).AND.(IIAUX.LE.9.OR.IIAUX.GT.13))
     .   CALL RUNEND('ERROR: WRONG NUMBER OF MODEL FOR ISOT. HARDENING')
        IF(IIAUX.EQ.10.OR.IIAUX.EQ.11) THEN
         IF(KDYNA.EQ.0)
     .    CALL RUNMEN('WARNING: ISOT. HARD. MODEL=10 & STATIC ANALYSIS')
         KAUXX=1
         IF(KAUXX.GT.KPLA7) KPLA7=1
        ENDIF
C
C**** LOOK FOR 'HD_COEFFICIENT' CARD (HARDENING COEFFICIENT)
C
        IF(IIAUX.EQ.1.OR.IIAUX.EQ.2.OR.IIAUX.EQ.4.OR.IIAUX.EQ.5.OR.
     .     IIAUX.EQ.8.OR.IIAUX.EQ.10.OR.IIAUX.EQ.11.OR.IIAUX.EQ.12.OR.
     .     IIAUX.EQ.13) THEN
         NPRIN=0
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         IF((WORDS(1).EQ.'HD_MO').OR.(WORDS(1).EQ.'HD_CO')) THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          NLINE=INT(PARAM(1))
          IF(NLINE.GT.20.OR.NLINE.LT.0)
     .     CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN HARD. COEF.')
          NPRIN=1
          DO LINEA=1,NLINE
           CALL LISTEN('PROPIVT',NPRIN,ITAPE)
           DO I=1,2
            NPOI2=NPOI2+1
            PROPS(NPOI2)=PARAM(I)
           ENDDO
          ENDDO
          WRITE(LURES,810)
          DO IPROP=NPOI1+1,NPOI2,2
           WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
          ENDDO
          WRITE(LURES,899)
         ELSE
          CALL RUNEND('ERROR: HARDEN. COEFF. CARD NOT FOUND')
         ENDIF
        ENDIF
C
C**** LOOK FOR 'B_COEFFICIENT' CARD (B COEFFICIENT)
C
        IF(IIAUX.EQ.3.OR.IIAUX.EQ.6.OR.IIAUX.EQ.7.OR.IIAUX.EQ.8.OR.
     .     IIAUX.EQ.10.OR.IIAUX.EQ.11.OR.IIAUX.EQ.12.OR.
     .     IIAUX.EQ.13) THEN
         NPRIN=0
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'B_COE') THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          NLINE=INT(PARAM(1))
          IF(NLINE.GT.20.OR.NLINE.LT.0)
     .     CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN B COEF.')
          NPRIN=1
          DO LINEA=1,NLINE
           CALL LISTEN('PROPIVT',NPRIN,ITAPE)
           DO I=1,2
            NPOI2=NPOI2+1
            PROPS(NPOI2)=PARAM(I)
           ENDDO
          ENDDO
          WRITE(LURES,811)
          DO IPROP=NPOI1+1,NPOI2,2
           WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
          ENDDO
          WRITE(LURES,899)
         ELSE
          CALL RUNEND('ERROR: B COEFFICIENT. CARD NOT FOUND')
         ENDIF
C
C**** LOOK FOR 'Q_COEFFICIENT' CARD (Q COEFFICIENT)
C
         NPRIN=0
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'Q_COE') THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          NLINE=INT(PARAM(1))
          IF(NLINE.GT.20.OR.NLINE.LT.0)
     .     CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN Q COEF.')
          NPRIN=1
          DO LINEA=1,NLINE
           CALL LISTEN('PROPIVT',NPRIN,ITAPE)
           DO I=1,2
            NPOI2=NPOI2+1
            PROPS(NPOI2)=PARAM(I)
           ENDDO
          ENDDO
          WRITE(LURES,812)
          DO IPROP=NPOI1+1,NPOI2,2
           WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
          ENDDO
          WRITE(LURES,899)
         ELSE
          CALL RUNEND('ERROR: Q COEFFICIENT CARD NOT FOUND')
         ENDIF
        ENDIF
       ELSE
        CALL RUNEND('ERROR: NO MODEL FOR ISOTROPIC HARDENING')
       ENDIF
      ELSE
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI2)=0.0D0
       GO TO 2
      ENDIF
C
C**** LOOK FOR 'KINEMATIC_HARDENING' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIVT',NPRIN,ITAPE)
    2 IF(WORDS(1).EQ.'KINEM') THEN
       KAUXX=1
       IF(KAUXX.GT.KPLA2) KPLA2=1
       IF(WORDS(2).EQ.'MODEL') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        IIAUX=INT(PARAM(1))
        IF(IIAUX.LE.0.OR.IIAUX.GT.4)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF MODEL FOR ISOT. HARDENING')
C
C**** LOOK FOR 'HD_COEFFICIENT' CARD (HARDENING COEFFICIENT)
C
        IF(IIAUX.EQ.1.OR.IIAUX.EQ.2) THEN
         NPRIN=0
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         IF((WORDS(1).EQ.'HD_MO').OR.(WORDS(1).EQ.'HD_CO')) THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          NLINE=INT(PARAM(1))
          IF(NLINE.GT.20.OR.NLINE.LT.0)
     .     CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN HARD. COEF.')
          NPRIN=1
          DO LINEA=1,NLINE
           CALL LISTEN('PROPIVT',NPRIN,ITAPE)
           DO I=1,2
            NPOI2=NPOI2+1
            PROPS(NPOI2)=PARAM(I)
           ENDDO
          ENDDO
          WRITE(LURES,810)
          DO IPROP=NPOI1+1,NPOI2,2
           WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
          ENDDO
          WRITE(LURES,899)
         ELSE
          CALL RUNEND('ERROR: HARDEN. COEFF. CARD NOT FOUND')
         ENDIF
        ENDIF
C
C**** LOOK FOR 'B_COEFFICIENT' CARD (B COEFFICIENT)
C
        IF(IIAUX.EQ.3.OR.IIAUX.EQ.4) THEN
         NPRIN=0
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'B_COE') THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          NLINE=INT(PARAM(1))
          IF(NLINE.GT.20.OR.NLINE.LT.0)
     .     CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN B COEF.')
          NPRIN=1
          DO LINEA=1,NLINE
           CALL LISTEN('PROPIVT',NPRIN,ITAPE)
           DO I=1,2
            NPOI2=NPOI2+1
            PROPS(NPOI2)=PARAM(I)
           ENDDO
          ENDDO
          WRITE(LURES,811)
          DO IPROP=NPOI1+1,NPOI2,2
           WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
          ENDDO
          WRITE(LURES,899)
         ELSE
          CALL RUNEND('ERROR: B COEFFICIENT. CARD NOT FOUND')
         ENDIF
C
C**** LOOK FOR 'Q_COEFFICIENT' CARD (Q COEFFICIENT)
C
         NPRIN=0
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'Q_COE') THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          NLINE=INT(PARAM(1))
          IF(NLINE.GT.20.OR.NLINE.LT.0)
     .     CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN Q COEF.')
          NPRIN=1
          DO LINEA=1,NLINE
           CALL LISTEN('PROPIVT',NPRIN,ITAPE)
           DO I=1,2
            NPOI2=NPOI2+1
            PROPS(NPOI2)=PARAM(I)
           ENDDO
          ENDDO
          WRITE(LURES,812)
          DO IPROP=NPOI1+1,NPOI2,2
           WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
          ENDDO
          WRITE(LURES,899)
         ELSE
          CALL RUNEND('ERROR: Q COEFFICIENT CARD NOT FOUND')
         ENDIF
        ENDIF
       ELSE
        CALL RUNEND('ERROR: NO MODEL FOR KINEMATIC HARDENING')
       ENDIF
      ELSE
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI2)=0.0D0
       GO TO 3
      ENDIF
C
C**** LOOK FOR 'RECRYSTALLIZATION' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIVT',NPRIN,ITAPE)
    3 IF(WORDS(1).EQ.'RECRY') THEN
       KAUXX=1
       IF(KAUXX.GT.KPLA8) KPLA8=1
       IF(WORDS(2).EQ.'MODEL') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        IRECR=INT(PARAM(1))
        IF(IRECR.LE.0.OR.IRECR.GT.1)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF MODEL FOR RECRYSTALLIZAT.')
C
C**** LOOK FOR 'INITIAL_VALUES' CARD
C
        IF(IRECR.EQ.1) THEN
         IF(IVIFL.NE.5)
     .    CALL RUNEND('ERROR: IRECR=1 ONLY WORKS WITH IVIFL=5')
         IF(MBUSS.LT.20)                                 ! see auxl_om.f
     .    CALL RUNEND('ERROR: MBUSS MUST GE 20')
         NPRIN=0
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'INITI') THEN
          NPRIN=1
          CALL LISTEN('PROPIVT',NPRIN,ITAPE)
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI2)=PARAM(1)          ! deformation resistance So
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(2)          ! grain size Lo
          WRITE(LURES,857) PROPS(NPOI1)
          WRITE(LURES,858) PROPS(NPOI1+1)
          INITP=1       ! to be revised  ! ic index for int. variables
         ELSE
          CALL RUNEND('PROPIVT: INITIAL_VALUES CARD NOT FOUND')
         ENDIF
C
C**** LOOK FOR 'DEFORMATION_RESISTANCE' CARD
C
         NPRIN=0
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'DEFOR') THEN
          NPRIN=1
          CALL LISTEN('PROPIVT',NPRIN,ITAPE)
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI2)=PARAM(1)          ! hs
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(2)          ! as
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(3)          ! x
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(4)          ! n2
          WRITE(LURES,859) PROPS(NPOI1)
          WRITE(LURES,860) PROPS(NPOI1+1)
          WRITE(LURES,861) PROPS(NPOI1+2)
          WRITE(LURES,862) PROPS(NPOI1+3)
         ELSE
          CALL RUNEND('PROPIVT: DEFORMATION_RESISTANCE CARD NOT FOUND')
         ENDIF
C
C**** LOOK FOR 'GRAIN_GROWTH' CARD
C
         NPRIN=0
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'GRAIN') THEN
          NPRIN=1
          CALL LISTEN('PROPIVT',NPRIN,ITAPE)
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI2)=PARAM(1)          ! fR
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(2)          ! fC
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(3)          ! fX
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(4)          ! CR
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(5)          ! CC
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(6)          ! p
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(7)          ! q
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(8)          ! Lo_dot
          WRITE(LURES,863) PROPS(NPOI1)
          WRITE(LURES,864) PROPS(NPOI1+1)
          WRITE(LURES,865) PROPS(NPOI1+2)
          WRITE(LURES,866) PROPS(NPOI1+3)
          WRITE(LURES,867) PROPS(NPOI1+4)
          WRITE(LURES,868) PROPS(NPOI1+5)
          WRITE(LURES,869) PROPS(NPOI1+6)
          WRITE(LURES,870) PROPS(NPOI1+7)
         ELSE
          CALL RUNEND('PROPIVT: GRAIN_GROWTH_PARAMETERS CARD NOT FOUND')
         ENDIF
        ENDIF                ! irecr.eq.1
       ELSE
        CALL RUNEND('ERROR: NO MODEL FOR RECRYSTALLIZATION')
       ENDIF
      ELSE
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI2)=0.0D0    ! IRECR=0
       GO TO 4
      ENDIF
C
C**** LOOK FOR 'YIELD_FUNCTION' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIVT',NPRIN,ITAPE)
    4 IF(WORDS(1).EQ.'YIELD') THEN
       IAUXX=0
       IF(WORDS(2).EQ.'VON_M') THEN
        IAUXX=1
        PROPS(36)=32.0D0
        NCRIT=32
        IF(WORDS(3).EQ.'VERSI') THEN
         PROPS(35)=PARAM(1)
         IIAUX=INT(PROPS(35))
         IF(IIAUX.LE.0.OR.IIAUX.GT.4)
     .    CALL RUNEND('ERROR IN VERSION OF VON MISES YIELD CRIT.')
         IF(IIAUX.EQ.4.AND.NPRE3.EQ.0) INITP=1  ! ic index for int. var.
        ELSE
         PROPS(35)=1.0D0       ! default value for version of Von Mises
        ENDIF
       ENDIF
       IF(WORDS(2).EQ.'GURSO') THEN
        KAUXX=1
        IF(KAUXX.GT.KPLA4) KPLA4=1
        IAUXX=1
        PROPS(36)=41.0D0
        NCRIT=41
        IF(WORDS(3).EQ.'VERSI') THEN
         PROPS(35)=PARAM(1)
         IIAUX=INT(PROPS(35))
         IF(IIAUX.LE.0.OR.IIAUX.GT.4)
     .    CALL RUNEND('ERROR IN VERSION OF GURSON YIELD CRIT.')
        ELSE
         PROPS(35)=1.0D0       ! default value for version of Gurson
        ENDIF
       ENDIF
       IF(WORDS(2).EQ.'GREEN') THEN          ! Green sand model
        PROPS(36)=39.0D0
        NCRIT=39
        IF(WORDS(3).EQ.'VERSI') THEN
         PROPS(35)=PARAM(1)
         IIAUX=INT(PROPS(35))
         IF(IIAUX.LE.0.OR.IIAUX.GT.1)
     .    CALL RUNEND('ERROR IN VERSION OF GREEN SAND YIELD CRIT.')
        ELSE
         PROPS(35)=1.0D0       ! default value for version of Green sand
        ENDIF
C
        IF(IIAUX.EQ.1) THEN
C
C**** LOOK FOR 'COEFFICIENT_B1' CARD
C
         NPRIN=0
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'COEB1') THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          NLINE=INT(PARAM(1))
          IF(NLINE.GT.20)
     .     CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN COEF. B1')
          NPRIN=1
          DO LINEA=1,NLINE
           CALL LISTEN('PROPIVT',NPRIN,ITAPE)
           DO I=1,2
            NPOI2=NPOI2+1
            PROPS(NPOI2)=PARAM(I)
           ENDDO
          ENDDO
          WRITE(LURES,831)
          DO IPROP=NPOI1+1,NPOI2,2
           WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
          ENDDO
          WRITE(LURES,899)
         ELSE
          CALL RUNEND('PROPIVT: COEFFICIENT B1 CARD NOT FOUND')
         ENDIF
C
C**** LOOK FOR 'COEFFICIENT_B2' CARD
C
         NPRIN=0
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'COEB2') THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          NLINE=INT(PARAM(1))
          IF(NLINE.GT.20)
     .     CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN COEF. B2')
          NPRIN=1
          DO LINEA=1,NLINE
           CALL LISTEN('PROPIVT',NPRIN,ITAPE)
           DO I=1,2
            NPOI2=NPOI2+1
            PROPS(NPOI2)=PARAM(I)
           ENDDO
          ENDDO
          WRITE(LURES,832)
          DO IPROP=NPOI1+1,NPOI2,2
           WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
          ENDDO
          WRITE(LURES,899)
         ELSE
          CALL RUNEND('PROPIVT: COEFFICIENT B2 CARD NOT FOUND')
         ENDIF
C
C**** LOOK FOR 'COEFFICIENT_B3' CARD
C
         NPRIN=0
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'COEB3') THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          NLINE=INT(PARAM(1))
          IF(NLINE.GT.20)
     .     CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN COEF. B3')
          NPRIN=1
          DO LINEA=1,NLINE
           CALL LISTEN('PROPIVT',NPRIN,ITAPE)
           DO I=1,2
            NPOI2=NPOI2+1
            PROPS(NPOI2)=PARAM(I)
           ENDDO
          ENDDO
          WRITE(LURES,833)
          DO IPROP=NPOI1+1,NPOI2,2
           WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
          ENDDO
          WRITE(LURES,899)
         ELSE
          CALL RUNEND('PROPIVT: COEFFICIENT B3 CARD NOT FOUND')
         ENDIF
C
C**** LOOK FOR 'COEFFICIENT_B4' CARD
C
         NPRIN=0
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'COEB4') THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          NLINE=INT(PARAM(1))
          IF(NLINE.GT.20)
     .     CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN COEF. B4')
          NPRIN=1
          DO LINEA=1,NLINE
           CALL LISTEN('PROPIVT',NPRIN,ITAPE)
           DO I=1,2
            NPOI2=NPOI2+1
            PROPS(NPOI2)=PARAM(I)
           ENDDO
          ENDDO
          WRITE(LURES,834)
          DO IPROP=NPOI1+1,NPOI2,2
           WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
          ENDDO
          WRITE(LURES,899)
         ELSE
          CALL RUNEND('PROPIVT: COEFFICIENT B4 CARD NOT FOUND')
         ENDIF
C
C**** LOOK FOR 'COEFFICIENT_B5' CARD
C
         NPRIN=0
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'COEB5') THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          NLINE=INT(PARAM(1))
          IF(NLINE.GT.20)
     .     CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN COEF. B5')
          NPRIN=1
          DO LINEA=1,NLINE
           CALL LISTEN('PROPIVT',NPRIN,ITAPE)
           DO I=1,2
            NPOI2=NPOI2+1
            PROPS(NPOI2)=PARAM(I)
           ENDDO
          ENDDO
          WRITE(LURES,835)
          DO IPROP=NPOI1+1,NPOI2,2
           WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
          ENDDO
          WRITE(LURES,899)
         ELSE
          CALL RUNEND('PROPIVT: COEFFICIENT B5 CARD NOT FOUND')
         ENDIF
        ENDIF                   ! iiaux.eq.1
       ENDIF                    ! words(2).eq.'green'
       IF(IAUXX.EQ.0)
     .  CALL RUNEND('PROPIVT: WRONG YIELD FUNCTION         ')
      ELSE
       CALL RUNEND('PROPIVT: YIELD FUNCTION CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'FLOW_POTENTIAL' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIVT',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'FLOW_') THEN
       IAUXX=0
       IF(WORDS(2).EQ.'VON_M') THEN
        IAUXX=1
        PROPS(52)=32.0D0
        NCRIP=32
        PROPS(56)=1.0D0         ! assumption
        KFLUJ=1
        IF(WORDS(3).EQ.'VERSI') THEN
         PROPS(51)=PARAM(1)
         IIAUX=INT(PROPS(51))
         IF(IIAUX.LE.0.OR.IIAUX.GT.4)
     .    CALL RUNEND('ERROR IN VERSION OF VON MISES FLOW POT.')
        ELSE
         PROPS(51)=1.0D0        ! default value for version of Von Mises
        ENDIF
       ENDIF
       IF(WORDS(2).EQ.'GURSO') THEN
        KAUXX=1
        IF(KAUXX.GT.KPLA4) KPLA4=1
        IAUXX=1
        PROPS(52)=41.0D0
        NCRIP=41
        IF(WORDS(3).EQ.'VERSI') THEN
         PROPS(51)=PARAM(1)
         IIAUX=INT(PROPS(51))
         IF(IIAUX.LE.0.OR.IIAUX.GT.4)
     .    CALL RUNEND('ERROR IN VERSION OF GURSON FLOW POT.')
        ELSE
         PROPS(51)=1.0D0        ! default value for version of Gurson
        ENDIF
       ENDIF
       IF(WORDS(2).EQ.'GREEN') THEN          ! Green sand model
        PROPS(52)=40.0D0
        NCRIP=40
        IF(WORDS(3).EQ.'VERSI') THEN
         PROPS(51)=PARAM(1)
         IIAUX=INT(PROPS(51))
         IF(IIAUX.LE.0.OR.IIAUX.GT.1)
     .    CALL RUNEND('ERROR IN VERSION OF DRUCKER-PRAGER FLOW POT.')
        ELSE
         PROPS(51)=1.0D0        ! default value for version of Green sand
        ENDIF
C
        IF(IIAUX.EQ.1) THEN
C
C**** LOOK FOR 'COEFFICIENT_B6' CARD
C
         NPRIN=0
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'COEB6') THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          NLINE=INT(PARAM(1))
          IF(NLINE.GT.20)
     .     CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN COEF. B6')
          NPRIN=1
          DO LINEA=1,NLINE
           CALL LISTEN('PROPIVT',NPRIN,ITAPE)
           DO I=1,2
            NPOI2=NPOI2+1
            PROPS(NPOI2)=PARAM(I)
           ENDDO
          ENDDO
          WRITE(LURES,836)
          DO IPROP=NPOI1+1,NPOI2,2
           WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
          ENDDO
          WRITE(LURES,899)
         ELSE
          CALL RUNEND('PROPIVT: COEFFICIENT B6 CARD NOT FOUND')
         ENDIF
C
C**** LOOK FOR 'COEFFICIENT_B7' CARD
C
         NPRIN=0
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'COEB7') THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          NLINE=INT(PARAM(1))
          IF(NLINE.GT.20)
     .     CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN COEF. B7')
          NPRIN=1
          DO LINEA=1,NLINE
           CALL LISTEN('PROPIVT',NPRIN,ITAPE)
           DO I=1,2
            NPOI2=NPOI2+1
            PROPS(NPOI2)=PARAM(I)
           ENDDO
          ENDDO
          WRITE(LURES,837)
          DO IPROP=NPOI1+1,NPOI2,2
           WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
          ENDDO
          WRITE(LURES,899)
         ELSE
          CALL RUNEND('PROPIVT: COEFFICIENT B7 CARD NOT FOUND')
         ENDIF
C
C**** LOOK FOR 'COEFFICIENT_B8' CARD
C
         NPRIN=0
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'COEB8') THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          NLINE=INT(PARAM(1))
          IF(NLINE.GT.20)
     .     CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN COEF. B8')
          NPRIN=1
          DO LINEA=1,NLINE
           CALL LISTEN('PROPIVT',NPRIN,ITAPE)
           DO I=1,2
            NPOI2=NPOI2+1
            PROPS(NPOI2)=PARAM(I)
           ENDDO
          ENDDO
          WRITE(LURES,838)
          DO IPROP=NPOI1+1,NPOI2,2
           WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
          ENDDO
          WRITE(LURES,899)
         ELSE
          CALL RUNEND('PROPIVT: COEFFICIENT B8 CARD NOT FOUND')
         ENDIF
        ENDIF                   ! iiaux.eq.1
       ENDIF                    ! words(2).eq.'green'
       IF(IAUXX.EQ.0)
     .  CALL RUNEND('PROPIVT: WRONG FLOW POTENTIAL FUNCTION         ')
      ELSE
       CALL RUNEND('PROPIVT: FLOW_POTENTIAL CARD NOT FOUND')
      ENDIF
C
      IF(NCRIT.EQ.NCRIP) THEN
       IF(KSYMM.EQ.0)
     .  CALL RUNMEN('WARNING: ASSOCIATE PLAST. & NON SYMMETRIC SOLVER')
      ELSE
       IF(KSYMM.EQ.1)
     .  CALL RUNMEN('WARNING: NON ASSOCIATE PLAST. & SYMMETRIC SOLVER')
      ENDIF
C
      IF(IVIFL.EQ.4) THEN
       IF(NCRIT.NE.32.OR.NCRIP.NE.32)
     .  CALL RUNEND('ERROR: IVIFL=4 & NCRIT ne 32 OR NCRIP ne 32')
      ENDIF
C
      IF(NCRIT.EQ.41.OR.NCRIP.EQ.41) THEN         ! Gurson
C
C**** LOOK FOR 'POROSITY' CARD
C
       NPRIN=0
       CALL LISTEN('PROPIVT',NPRIN,ITAPE)
       IF(WORDS(1).EQ.'POROS')THEN
        IF(WORDS(2).EQ.'MODEL')THEN
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=PARAM(1)
         IIAUX=INT(PARAM(1))
         IF(IIAUX.LE.0.OR.IIAUX.GT.4)
     .    CALL RUNEND('ERROR: WRONG NUMBER OF MODEL FOR POROSITY')
         IF(ISOTT.EQ.3)
     .    CALL RUNMEN('WARNING: ISOT. MODEL 3 WITH POROSITY')
C
C**** LOOK FOR 'POROSITY_PARAMETER' CARD
C
         IF(IIAUX.EQ.1.OR.IIAUX.EQ.2.OR.IIAUX.EQ.3.OR.IIAUX.EQ.4) THEN
          IF(IIAUX.EQ.1.OR.IIAUX.EQ.3) THEN
           IF(KSYMM.EQ.1)
     .      CALL RUNMEN('WARNING: NON ASSOCIATE PLAST. & SYM. SOLVER')
          ENDIF
          NPRIN=0
          CALL LISTEN('PROPIPT',NPRIN,ITAPE)
          IF(WORDS(1).EQ.'POROS') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20)
     .      CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN POROS. PAR.')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPIPT',NPRIN,ITAPE)
            DO I=1,2
             NPOI2=NPOI2+1
             PROPS(NPOI2)=PARAM(I)
            ENDDO
           ENDDO
           WRITE(LURES,820)
           DO IPROP=NPOI1+1,NPOI2,2
            WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
           ENDDO
           WRITE(LURES,899)
          ELSE
           CALL RUNEND('PROPIPT: POROSITY_PARAMETER CARD NOT FOUND')
          ENDIF
         ENDIF
C
C**** LOOK FOR 'Q1_PARAMETER' CARD
C
         IF(IIAUX.EQ.4) THEN
          NPRIN=0
          CALL LISTEN('PROPIPT',NPRIN,ITAPE)
          IF(WORDS(1).EQ.'Q1_PA') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20)
     .      CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN Q1_PAR.')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPIPT',NPRIN,ITAPE)
            DO I=1,2
             NPOI2=NPOI2+1
             PROPS(NPOI2)=PARAM(I)
            ENDDO
           ENDDO
           WRITE(LURES,821)
           DO IPROP=NPOI1+1,NPOI2,2
            WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
           ENDDO
           WRITE(LURES,899)
          ELSE
           CALL RUNEND('PROPIPT: Q1_PARAMETER CARD NOT FOUND')
          ENDIF
         ENDIF
C
C**** LOOK FOR 'Q2_PARAMETER' CARD
C
         IF(IIAUX.EQ.4) THEN
          NPRIN=0
          CALL LISTEN('PROPIPT',NPRIN,ITAPE)
          IF(WORDS(1).EQ.'Q2_PA') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20)
     .      CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN Q1_PAR.')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPIPT',NPRIN,ITAPE)
            DO I=1,2
             NPOI2=NPOI2+1
             PROPS(NPOI2)=PARAM(I)
            ENDDO
           ENDDO
           WRITE(LURES,822)
           DO IPROP=NPOI1+1,NPOI2,2
            WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
           ENDDO
           WRITE(LURES,899)
          ELSE
           CALL RUNEND('PROPIPT: Q2_PARAMETER CARD NOT FOUND')
          ENDIF
         ENDIF
        ELSE
         CALL RUNEND('ERROR: NO MODEL FOR POROSITY')
        ENDIF
       ELSE
        CALL RUNEND('ERROR: POROSITY CARD NOT FOUND')
       ENDIF
      ENDIF                            ! ncrit.eq.41.or.ncrip.eq.41
C
C**** LOOK FOR 'PHASE_CHANGE' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIVT',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'PHASE') THEN
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI1)=PARAM(1)           ! nplatm
       NLINE=INT(PROPS(NPOI1))
       IF(NLINE.GT.5)
     .  CALL RUNEND('ERROR: WRONG NUMBER OF PHASE_CHANGES')
C
       IF(ITERME.GE.0) THEN
        NLINEA=NFPCH/2
        IF(NLINE.GT.NLINEA)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF PHASE_CHANGES')
        CALL RUNMEN('WARNING: CHECK THE CORRESPONDANCE BETWEEN THERMAL &
     . MECHANICAL PHASE-CHANGES    ')
       ENDIF                           ! iterme.gt.0
C
       NPRIN=1
       ILSX=0
       DO LINEA=1,NLINE
        WRITE(LURES,813) LINEA
        CALL LISTEN('PROPIVT',NPRIN,ITAPE)
C
        IF(WORDS(1).EQ.'SS_PH') THEN
C
         NPOI2=NPOI2+1
         PROPS(NPOI2)=0.0D0             ! ITEMRO
         NPOI2=NPOI2+1
         PROPS(NPOI2)=0.0D0
         NPOI2=NPOI2+1
         PROPS(NPOI2)=0.0D0
         NPOI2=NPOI2+1
         PROPS(NPOI2)=0.0D0
C
         NPOI2=NPOI2+1                  ! takes into account sense of pc
         PROPS(NPOI2)=0.0D0             ! ISEPC
         IF(WORDS(2).EQ.'HEATI') PROPS(NPOI2)=1.0D0
C
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'    ') THEN    ! constant expansion
          NPOI2=NPOI2+1
          PROPS(NPOI2)=0.0D0            ! IPCFOM
          NPOI2=NPOI2+1
          PROPS(NPOI2)=0.0D0            ! IPCMOM
          NPOI2=NPOI2+1
          PROPS(NPOI2)=0.0D0            ! ILSPC
          NPOI2=NPOI2+1
          PROPS(NPOI2)=LINEA            ! ISSPC
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(1)         ! expansion
          WRITE(LURES,816)
          WRITE(LURES,815) PROPS(NPOI2)
          WRITE(LURES,899)
         ENDIF
        ENDIF
C
        IF(WORDS(1).EQ.'LS_PH') THEN
         ILSX=ILSX+1
         IF(ILSX.GT.1)
     .    CALL RUNEND('PROPIVT: TWO OR MORE LS PCS ARE NOT IMPLEMENTED')
C
C**** OPTION TO USE A LIQUID-SOLID PHASE-CHANGE FUNCTION DEFINED IN
C     MECHANICAL INPUT FILE
C
         IF(WORDS(2).EQ.'T_SOL'.AND.WORDS(3).EQ.'T_LIQ') THEN
          NPOI2=NPOI2+1
          PROPS(NPOI2)=1.0D0            ! ITEMRO
          ITEMRO=1
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(1)         ! TEMPS
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(2)         ! TEMPL
         ELSE
          NPOI2=NPOI2+1
          PROPS(NPOI2)=0.0D0            ! ITEMRO
          ITEMRO=0
          NPOI2=NPOI2+1
          PROPS(NPOI2)=0.0D0
          NPOI2=NPOI2+1
          PROPS(NPOI2)=0.0D0
         ENDIF
C
C**** SHRINKAGE EFFECT
C
         IF(WORDS(2).EQ.'SHRIN'.OR.WORDS(4).EQ.'SHRIN') THEN
          NPOI2=NPOI2+1
          PROPS(NPOI2)=1.0D0            ! ISHRI
          ISHRI=1
         ELSE
          NPOI2=NPOI2+1
          PROPS(NPOI2)=0.0D0
          ISHRI=0
         ENDIF
C
         NPOI2=NPOI2+1                  ! takes into account sense of pc
         PROPS(NPOI2)=0.0D0             ! ISEPC
         IF(WORDS(2).EQ.'HEATI'.OR.WORDS(3).EQ.'HEATI'.OR.
     .      WORDS(5).EQ.'HEATI') PROPS(NPOI2)=1.0D0
C
         CALL LISTEN('PROPIVT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'    ') THEN    ! constant expansion
          NPOI2=NPOI2+1
          PROPS(NPOI2)=0.0D0            ! IPCFOM
          NPOI2=NPOI2+1
          PROPS(NPOI2)=0.0D0            ! IPCMOM
          NPOI2=NPOI2+1
          PROPS(NPOI2)=LINEA            ! ILSPC
          NPOI2=NPOI2+1
          PROPS(NPOI2)=0.0D0            ! ISSPC
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(1)         ! expansion
          WRITE(LURES,814)
          IF(ITEMRO.EQ.1) WRITE(LURES,819) PROPS(NPOI2-6),PROPS(NPOI2-5)
          WRITE(LURES,815) PROPS(NPOI2)
          WRITE(LURES,899)
         ENDIF
C
         IF(ISHRI.EQ.1) THEN
          KPLA5=1
          WRITE(LURES,826)
          WRITE(LURES,899)
         ENDIF
C
        ENDIF              ! words(1).eq.'ls_ph'
C
       ENDDO               ! linea=1,nline
C
      ELSE
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI1)=0.0D0  ! nplatm
       GO TO 5
      ENDIF
C
C**** LOOK FOR 'DAMAGE_PARAMETERS' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIVT',NPRIN,ITAPE)
    5 IF(WORDS(1).EQ.'DAMAG')THEN
       KAUXX=1
       IF(KAUXX.GT.KPLA3) KPLA3=1
       IF(WORDS(2).EQ.'MODEL')THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)      ! idamg
        IDAMG=INT(PROPS(NPOI1))
        IF(IDAMG.NE.1.AND.IDAMG.NE.2.AND.IDAMG.NE.3.AND.IDAMG.NE.4.AND.
     .     IDAMG.NE.5.AND.IDAMG.NE.21)
     .   CALL RUNEND('PROPIVT: WRONG DAMAGE MODEL')
       ELSE                        ! default: damage model=1
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=1.0D0
        IDAMG=1                    ! idamg
       ENDIF
C
       IF(IDAMG.EQ.1.OR.IDAMG.EQ.2) THEN
C
C**** LOOK FOR 'DAMAGE_STRESS' CARD
C
        NPRIN=0
        CALL LISTEN('PROPIVT',NPRIN,ITAPE)
        IF(WORDS(1).EQ.'DAMAG') THEN
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=PARAM(1)
         NLINE=INT(PARAM(1))
         IF(NLINE.GT.20)
     .    CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN DAMAGE STRESS')
         NPRIN=1
         DO LINEA=1,NLINE
          CALL LISTEN('PROPIVT',NPRIN,ITAPE)
          DO I=1,2
           NPOI2=NPOI2+1
           PROPS(NPOI2)=PARAM(I)
          ENDDO
         ENDDO
         WRITE(LURES,817)
         DO IPROP=NPOI1+1,NPOI2,2
          WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
         ENDDO
         WRITE(LURES,899)
        ELSE
         CALL RUNEND('PROPIVT: DAMAGE STRESS CARD NOT FOUND')
        ENDIF
C
C**** LOOK FOR 'FRACTURE ENERGY' CARD
C
        NPRIN=0
        CALL LISTEN('PROPIVT',NPRIN,ITAPE)
        IF(WORDS(1).EQ.'FRACT') THEN
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=PARAM(1)
         NLINE=INT(PARAM(1))
         IF(NLINE.GT.20)
     .    CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN FRACT. ENERGY')
         NPRIN=1
         DO LINEA=1,NLINE
          CALL LISTEN('PROPIVT',NPRIN,ITAPE)
          DO I=1,2
           NPOI2=NPOI2+1
           PROPS(NPOI2)=PARAM(I)
          ENDDO
         ENDDO
         WRITE(LURES,818)
         DO IPROP=NPOI1+1,NPOI2,2
          WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
         ENDDO
         WRITE(LURES,899)
        ELSE
         CALL RUNEND('PROPIVT: FRACTURE ENERGY CARD NOT FOUND')
        ENDIF
       ENDIF
C
       IF(IDAMG.EQ.3) THEN
C
C**** LOOK FOR 'COEFFICIENT_B9' CARD
C
        NPRIN=0
        CALL LISTEN('PROPIVT',NPRIN,ITAPE)
        IF(WORDS(1).EQ.'COEB9') THEN
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=PARAM(1)
         NLINE=INT(PARAM(1))
         IF(NLINE.GT.20)
     .    CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN COEF. B9')
         NPRIN=1
         DO LINEA=1,NLINE
          CALL LISTEN('PROPIVT',NPRIN,ITAPE)
          DO I=1,2
           NPOI2=NPOI2+1
           PROPS(NPOI2)=PARAM(I)
          ENDDO
         ENDDO
         WRITE(LURES,839)
         DO IPROP=NPOI1+1,NPOI2,2
          WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
         ENDDO
         WRITE(LURES,899)
        ELSE
         CALL RUNEND('PROPIVT: COEFFICIENT B9 CARD NOT FOUND')
        ENDIF
       ENDIF                  ! idamg.eq.3
C
       IF(IDAMG.EQ.4.OR.IDAMG.EQ.5) THEN
        call runend('ERROR: idamg=4,5 not implemented yet in propivt.f')
       ENDIF
C
       IF(IDAMG.EQ.21) THEN
        call runend('PROPIVT: IDAMG=21 not implemented')
       ENDIF
      ELSE
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI1)=0.0D0          ! idamg
       GO TO 6
      ENDIF
C
C**** LOOK FOR FREE ENERGY MODEL
C
      CALL LISTEN('PROPIPT',NPRIN,ITAPE)
    6 IF(WORDS(1).EQ.'FREE_') THEN
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI1)=PARAM(1)       ! ifren
       IFREN=INT(PROPS(NPOI1))
       IF(LARGE.EQ.0) THEN
        IF(IFREN.LE.0.OR.IFREN.GT.3)
     .   CALL RUNEND('ERROR: WRONG FREE ENERGY MODEL')
       ELSE
        IF(IFREN.LE.0.OR.IFREN.GT.10)
     .   CALL RUNEND('ERROR: WRONG FREE ENERGY MODEL')
        IF(NITERC.EQ.4) THEN
         IF(IFREN.EQ.2.OR.IFREN.EQ.3.OR.
     .      IFREN.EQ.5.OR.IFREN.EQ.6.OR.
     .      IFREN.EQ.7.OR.IFREN.EQ.8)
     .    CALL RUNEND('ERROR: IFREN=2,3,5,6,7,8 & NITERC=4 ARE INCOMP.')
        ENDIF
       ENDIF                       ! large.eq.0
       WRITE(LURES,827) IFREN
       WRITE(LURES,899)
      ELSE                         ! defaults
       IF(LARGE.EQ.0) THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=3.0D0         ! ifren
        IFREN=3
       ELSE
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=2.0D0         ! ifren
        IFREN=2
        IF(NITERC.EQ.4)
     .   CALL RUNEND('ERROR: IFREN=2 & NITERC=4 ARE INCOMPATIBLE')
       ENDIF
       WRITE(LURES,827) IFREN
       WRITE(LURES,899)
       GO TO 7
      ENDIF
C
C**** FATIGUE MODEL
C
      CALL LISTEN('PROPIPN',NPRIN,ITAPE)
    7 IF(WORDS(1).EQ.'FATIG') THEN
       IFATI=1                             ! fatigue material
       KAUXX=1
       IF(KAUXX.GT.KPLA6) KPLA6=1
C
       IF(WORDS(2).EQ.'MODEL') THEN
        IFATM=INT(PARAM(1))                ! fatigue model
C
        IF(IFATM.EQ.1) THEN
         CALL PROPFA1(ITAPE,PROPS)
        ENDIF
C
        IF(IFATM.GE.2) THEN
         CALL RUNEND('PROPIPN: IFATM GT 2')
        ENDIF
       ELSE
        CALL RUNEND('INPFAT: FATIGUE_MODEL CARD NOT FOUND')
       ENDIF
      ELSE
       GO TO 8
      ENDIF
C
C**** LOOK FOR 'END_MATERIAL' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIVT',NPRIN,ITAPE)
    8 IF(WORDS(1).NE.'END_M')
     . CALL RUNEND('PROPIVT: END_MATERIAL CARD NOT FOUND   ')
      WRITE(LURES,899)
C
C**** CONTROLS DIMENSION OF PROPS
C
      IF(NPOI2.GT.NPROP) THEN
       WRITE(LURES,900) NPROP,NPOI2
       CALL RUNEND('PROPIVT: TOO MANY PROPERTIES TO READ')
      ENDIF
C
      RETURN
  801 FORMAT(/,'  CONSTANT DENSITY VALUE=',E15.6,/)
  803 FORMAT(/,'  YOUNG MODULUS     Vs.     TEMPERATURE',/)
  804 FORMAT(E15.6,10X,E15.6)
  805 FORMAT(/,'  POISSON RATIO     Vs.     TEMPERATURE',/)
  806 FORMAT(/,'  THERMAL DILAT.    Vs.     TEMPERATURE',/)
  807 FORMAT(/,'  THERMAL HARDENING Vs.     TEMPERATURE',/)
  808 FORMAT(/,'  VISCOSITY FUNC.   Vs.     TEMPERATURE',/)
  809 FORMAT(/,'  EXPONENT FUNC.    Vs.     TEMPERATURE',/)
  810 FORMAT(/,'  HARDENING COEF.   Vs.     TEMPERATURE',/)
  811 FORMAT(/,'  B COEFFICIENT     Vs.     TEMPERATURE',/)
  812 FORMAT(/,'  Q COEFFICIENT     Vs.     TEMPERATURE',/)
  813 FORMAT(/,3X,'PHASE-CHANGE NUMBER  =',I5,/)
  814 FORMAT(/,3X,'LIQUID-SOLID PHASE-CHANGE',/)
  815 FORMAT(/,3X,'CONSTANT EXPANSION (PERCENTAGE) = ',E15.6,/)
  816 FORMAT(/,3X,'SOLID-SOLID PHASE-CHANGE',/)
  817 FORMAT(/,'  DAMAGE STRESS     Vs.     TEMPERATURE',/)
  818 FORMAT(/,'  FRACTURE ENERGY   Vs.     TEMPERATURE',/)
  819 FORMAT(/,3X,'SOLIDUS TEMPERATURE  =',E15.6,/,
     .         3X,'LIQUIDUS TEMPERATURE =',E15.6,/)
  820 FORMAT(/,'  POROSITY PARAM.   Vs.     TEMPERATURE',/)
  821 FORMAT(/,'  Q1 PARAMATER      Vs.     TEMPERATURE',/)
  822 FORMAT(/,'  Q2 PARAMATER      Vs.     TEMPERATURE',/)
  823 FORMAT(/,'  COEFFICIENT A     Vs.     TEMPERATURE',/)
  824 FORMAT(/,'  COEFFICIENT B     Vs.     TEMPERATURE',/)
  825 FORMAT(/,'  COEFFICIENT C     Vs.     TEMPERATURE',/)
  826 FORMAT(/,3X,'SHRINKAGE EFFECTS ARE CONSIDERED',/)
  827 FORMAT(/,3X,'FREE ENERGY MODEL=',I5,/)
  831 FORMAT(/,'  COEFFICIENT B1    Vs.     TEMPERATURE',/)
  832 FORMAT(/,'  COEFFICIENT B2    Vs.     TEMPERATURE',/)
  833 FORMAT(/,'  COEFFICIENT B3    Vs.     TEMPERATURE',/)
  834 FORMAT(/,'  COEFFICIENT B4    Vs.     TEMPERATURE',/)
  835 FORMAT(/,'  COEFFICIENT B5    Vs.     TEMPERATURE',/)
  836 FORMAT(/,'  COEFFICIENT B6    Vs.     TEMPERATURE',/)
  837 FORMAT(/,'  COEFFICIENT B7    Vs.     TEMPERATURE',/)
  838 FORMAT(/,'  COEFFICIENT B8    Vs.     TEMPERATURE',/)
  839 FORMAT(/,'  COEFFICIENT B9    Vs.     TEMPERATURE',/)
  851 FORMAT(/,'  MATERIAL CONSTANT Eo_dot      =',E15.6,/)
  852 FORMAT(/,'  MATERIAL CONSTANT w           =',E15.6,/)
  853 FORMAT(/,'  ACTIVATION ENERGY             =',E15.6,/)
  854 FORMAT(/,'  UNIVERSAL GAS CONSTANT        =',E15.6,/)
  855 FORMAT(/,'  MATERIAL CONSTANT chi         =',E15.6,/)
  856 FORMAT(/,'  MATERIAL CONSTANT n1          =',E15.6,/)
  857 FORMAT(/,'  INITIAL DEFORMATION RESISTANCE=',E15.6,/)
  858 FORMAT(/,'  INITIAL GRAIN SIZE            =',E15.6,/)
  859 FORMAT(/,'  MATERIAL CONSTANT hs          =',E15.6,/)
  860 FORMAT(/,'  MATERIAL CONSTANT as          =',E15.6,/)
  861 FORMAT(/,'  MATERIAL CONSTANT x           =',E15.6,/)
  862 FORMAT(/,'  MATERIAL CONSTANT n2          =',E15.6,/)
  863 FORMAT(/,'  MATERIAL CONSTANT fR          =',E15.6,/)
  864 FORMAT(/,'  MATERIAL CONSTANT fC          =',E15.6,/)
  865 FORMAT(/,'  MATERIAL CONSTANT fX          =',E15.6,/)
  866 FORMAT(/,'  MATERIAL CONSTANT CR          =',E15.6,/)
  867 FORMAT(/,'  MATERIAL CONSTANT CC          =',E15.6,/)
  868 FORMAT(/,'  MATERIAL CONSTANT p           =',E15.6,/)
  869 FORMAT(/,'  MATERIAL CONSTANT q           =',E15.6,/)
  870 FORMAT(/,'  MATERIAL CONSTANT Lo_dot      =',E15.6,/)
  899 FORMAT(/)
  900 FORMAT(//,'TOO MANY MATERIAL PROPERTIES TO READ:',/,
     .      20X,'NPROP =',I5,5X,'NUMBE =',I5,/)
      END
