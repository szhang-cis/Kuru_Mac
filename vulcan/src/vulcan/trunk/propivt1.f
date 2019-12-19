      SUBROUTINE PROPIVT1(ITAPE,PROPS)
C***********************************************************************
C
C**** THIS ROUTINE READS THE MATERIAL PROPERTIES FOR THE ISOTROPIC, 
C     VISCOPLASTIC & TEMPERATURE-DEPENDENT SG CAST IRON SOLID MODEL
C     (ELEMENT NO. 30)
C
C**** OUTPUT (see pointe.f)
C
C     KPLA1: isotropic hardening
C     KPLA2: kinematic hardening
C     KPLA5: shrinkage
C     KPLA6: fatigue
C     KPLA7: strain-rate-dependent hardening models
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** THERMAL VARIABLES
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
      PROPS(1)=0.0             ! isotropic
      PROPS(2)=30.0            ! viscoplastic
      PROPS(3)=2.0             ! temperature-dependent
      PROPS(4)=1.0             ! sg cast iron
      PROPS(59)=3.0            ! solidification
C
      PROPS(37)=0.0            ! non prescribed strains (see incstn.f)
      PROPS(19)=0.0
C
      PROPS(21)=0.0            ! no external damping
      PROPS(22)=0.0
C
C**** READ & WRITE MATERIAL PROPERTIES
C
      NPRIN=0
C
C**** LOOK FOR (INITIAL) 'DENSITY' CARD
C
      CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'DENSI') THEN
       PROPS(12)=PARAM(1)
       WRITE(LURES,801) PROPS(12)
       WRITE(LURES,899)
      ELSE
       IF(KDYNA.EQ.1) THEN
        CALL RUNEND('PROPIVT1: DENSITY CARD NOT FOUND')
       ELSE
        PROPS(12)=0.0
        CALL RUNMEN('WARNING: ZERO DENSITY VALUE IS ASSUMED')
        GO TO 1
       ENDIF
      ENDIF
C
C**** LOOK FOR 'YOUNG_MODULUS' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
    1 IF(WORDS(1).EQ.'YOUNG') THEN
       IF(WORDS(2).EQ.'FERRI') THEN
        NPOI1=61
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN YOUNG MODULUS')
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
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
        CALL RUNEND('PROPIVT1: FERRITE CARD NOT FOUND')
       ENDIF
      ELSE
       CALL RUNEND('PROPIVT1: YOUNG_MODULUS CARD NOT FOUND')
      ENDIF
C
      NPRIN=0
      CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'YOUNG') THEN
       IF(WORDS(2).EQ.'AUSTE') THEN
        NPOI1=NPOI2+1
	NPOI2=NPOI1
	PROPS(NPOI1)=PARAM(1)
	NLINE=INT(PARAM(1))
	IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN YOUNG MODULUS')
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
         DO I=1,2
          NPOI2=NPOI2+1
	  PROPS(NPOI2)=PARAM(I)
         ENDDO
	ENDDO
	WRITE(LURES,813)
	DO IPROP=NPOI1+1,NPOI2,2
	 WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
	ENDDO
	WRITE(LURES,899)
       ELSE
        CALL RUNEND('PROPIVT1: AUSTENITE CARD NOT FOUND')
       ENDIF
      ELSE
       CALL RUNEND('PROPIVT1: YOUNG_MODULUS CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'POISSON_RATIO' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'POISS')THEN
       IF(WORDS(2).EQ.'FERRI') THEN
	NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN POISSON_RATIO')
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
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
        CALL RUNEND('PROPIVT1: FERRITE CARD NOT FOUND')
       ENDIF
      ELSE
       CALL RUNEND('PROPIVT1: POISSON_RATIO CARD NOT FOUND')
      ENDIF
C
      NPRIN=0
      CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'POISS')THEN
       IF(WORDS(2).EQ.'AUSTE') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN POISSON_RATIO')
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
         DO I=1,2
	  NPOI2=NPOI2+1
	  PROPS(NPOI2)=PARAM(I)
	 ENDDO
	ENDDO
        WRITE(LURES,815)
        DO IPROP=NPOI1+1,NPOI2,2
         WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
        ENDDO
        WRITE(LURES,899)
       ELSE
        CALL RUNEND('PROPIVT1: AUSTENITE CARD NOT FOUND')
       ENDIF
      ELSE
       CALL RUNEND('PROPIVT1: POISSON_RATIO CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'REFERENCE_TEMPERATURE' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'REFER') THEN
       PROPS(11)=PARAM(1)
      ELSE
       CALL RUNEND('PROPIVT1: REFERENCE TEMPERATURE CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'THERMAL_DILATATION' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'THERM')THEN
       IF(WORDS(2).EQ.'FERRI') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN THERMAL DIL.')
C
        NPOI2=NPOI2+1
        PROPS(NPOI2)=0.0                            ! nalpuf
        IF(WORDS(2).EQ.'TANGE') PROPS(NPOI2)=1.0    ! nalpuf
C
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
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
        CALL RUNEND('PROPIVT1: FERRITE CARD NOT FOUND')
       ENDIF
      ELSE
       CALL RUNEND('PROPIVT1: THERMAL_DILAT. CARD NOT FOUND')
      ENDIF
C
      NPRIN=0
      CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'THERM')THEN
       IF(WORDS(2).EQ.'AUSTE') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN THERMAL DIL.')
C
        NPOI2=NPOI2+1
        PROPS(NPOI2)=0.0                            ! nalpua
        IF(WORDS(2).EQ.'TANGE') PROPS(NPOI2)=1.0    ! nalpua
C
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
         DO I=1,2
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(I)
         ENDDO
        ENDDO
        WRITE(LURES,816)
        DO IPROP=NPOI1+2,NPOI2,2
         WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
        ENDDO
	WRITE(LURES,899)
       ELSE
        CALL RUNEND('PROPIVT1: AUSTENITE CARD NOT FOUND')
       ENDIF
      ELSE
       CALL RUNEND('PROPIVT1: THERMAL_DILAT. CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'TH_HD_FUNCTION' CARD (THERMAL HARDENING)
C
      NPRIN=0
      CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
      IF((WORDS(1).EQ.'HD_FU').OR.(WORDS(1).EQ.'TH_HD'))THEN
       IF(WORDS(2).EQ.'FERRI') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN TH. HD. FUNC.')
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
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
        CALL RUNEND('PROPIVT1: FERRITE CARD NOT FOUND')
       ENDIF
      ELSE
       CALL RUNEND('PROPIVT1: THERMAL HARDEN. CARD NOT FOUND')
      ENDIF
C
      NPRIN=0
      CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
      IF((WORDS(1).EQ.'HD_FU').OR.(WORDS(1).EQ.'TH_HD'))THEN
       IF(WORDS(2).EQ.'AUSTE') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN TH. HD. FUNC.')
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
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
        CALL RUNEND('PROPIVT1: AUSTENITE CARD NOT FOUND')
       ENDIF
      ELSE
       CALL RUNEND('PROPIVT1: THERMAL HARDEN. CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'FLOW_PARAMETERS' CARD (PARAMETERS OF VISCOUS LAW)
C
      NPRIN=0
      CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'FLOW_')THEN
       IF(WORDS(2).EQ.'MODEL')THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)      ! ivifl
        IVIFL=INT(PARAM(1))
       ELSE                        ! default: flow model=1
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=1.0
        IVIFL=1
       ENDIF
      ELSE
       CALL RUNEND('PROPIVT1: FLOW_PARAMETERS CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'VISCOSITY' CARD
C
      IF(IVIFL.EQ.1.OR.IVIFL.EQ.2.OR.IVIFL.EQ.3) THEN
       NPRIN=0
       CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
       IF(WORDS(1).EQ.'VISCO') THEN
        IF(WORDS(2).EQ.'FERRI') THEN
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=PARAM(1)
         NLINE=INT(PARAM(1))
         IF(NLINE.GT.20)
     .    CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN VISCOSITY')
         NPRIN=1
         DO LINEA=1,NLINE
          CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
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
         CALL RUNEND('PROPIVT1: FERRITE CARD NOT FOUND')
        ENDIF
       ELSE
        CALL RUNEND('PROPIVT1: VISCOSITY CARD NOT FOUND')
       ENDIF
C
       NPRIN=0
       CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
       IF(WORDS(1).EQ.'VISCO') THEN
        IF(WORDS(2).EQ.'AUSTE') THEN
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=PARAM(1)
         NLINE=INT(PARAM(1))
         IF(NLINE.GT.20)
     .    CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN VISCOSITY')
         NPRIN=1
         DO LINEA=1,NLINE
          CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
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
         CALL RUNEND('PROPIVT1: AUSTENITE CARD NOT FOUND')
        ENDIF
       ELSE
        CALL RUNEND('PROPIVT1: VISCOSITY CARD NOT FOUND')
       ENDIF
C
      ENDIF               ! ivifl.eq.1.or.ivifl.eq.2.or.ivifl.eq.3
C
C**** LOOK FOR ZABARAS' COEFFICIENTS (A, B & C)
C
      IF(IVIFL.EQ.4)
     . CALL RUNEND('ERROR: IVIFL=4 NOT IMPLEMENTED - see propivt.f')
C
C**** LOOK FOR 'EXPONENT' CARD (EXPONENT OF VISCOUS LAW)
C
      NPRIN=0
      CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'EXPON')THEN
       IF(WORDS(2).EQ.'FERRI') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN EXPONENT')
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
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
        CALL RUNEND('PROPIVT1: FERRITE CARD NOT FOUND')
       ENDIF
      ELSE
       CALL RUNEND('PROPIVT1: EXPONENT CARD NOT FOUND')
      ENDIF
C
      NPRIN=0
      CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'EXPON')THEN
       IF(WORDS(2).EQ.'AUSTE') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN EXPONENT')
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
         DO I=1,2
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(I)
         ENDDO
        ENDDO
        WRITE(LURES,819)
        DO IPROP=NPOI1+1,NPOI2,2
         WRITE(LURES,804) PROPS(IPROP),PROPS(IPROP+1)
        ENDDO
        WRITE(LURES,899)
       ELSE
        CALL RUNEND('PROPIVT1: AUSTENITE CARD NOT FOUND')
       ENDIF
      ELSE
       CALL RUNEND('PROPIVT1: EXPONENT CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'ISOTROPIC_HARDENING' CARD
C
      NPRIN=0
      ISOTT=0
      CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'ISOTR')THEN
       KAUXX=1
       IF(KAUXX.GT.KPLA1) KPLA1=1
       IF(WORDS(2).EQ.'MODEL')THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        IIAUX=INT(PARAM(1))
        ISOTT=IIAUX
        IF((IIAUX.LE.0.OR.IIAUX.GT.8).AND.(IIAUX.LE.9.OR.IIAUX.GT.11))
     .   CALL RUNEND('ERROR: WRONG NUMBER OF MODEL FOR ISOT. HARDENING')
        IF(IIAUX.EQ.10) THEN
         IF(KDYNA.EQ.0)
     .   CALL RUNMEN('WARNING: ISOT. HARD. MODEL=10 & STATIC ANALYSIS')
         KAUXX=1
         IF(KAUXX.GT.KPLA7) KPLA7=1
        ENDIF
C
C**** LOOK FOR 'HD_COEFFICIENT' CARD (HARDENING COEFFICIENT)
C
        IF(IIAUX.EQ.1.OR.IIAUX.EQ.2.OR.IIAUX.EQ.4.OR.IIAUX.EQ.5.OR.
     .     IIAUX.EQ.8.OR.IIAUX.EQ.10.OR.IIAUX.EQ.11) THEN
         NPRIN=0
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
         IF((WORDS(1).EQ.'HD_MO').OR.(WORDS(1).EQ.'HD_CO'))THEN
          IF(WORDS(2).EQ.'FERRI') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20.OR.NLINE.LT.0)
     .      CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN HARD. COEF.')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
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
           CALL RUNEND('ERROR: FERRITE CARD NOT FOUND')
          ENDIF
         ELSE
          CALL RUNEND('ERROR: HARDEN. COEFF. CARD NOT FOUND')
         ENDIF
C
         NPRIN=0
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
         IF((WORDS(1).EQ.'HD_MO').OR.(WORDS(1).EQ.'HD_CO'))THEN
          IF(WORDS(2).EQ.'AUSTE') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20.OR.NLINE.LT.0)
     .      CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN HARD. COEF.')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
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
           CALL RUNEND('ERROR: AUSTENITE CARD NOT FOUND')
          ENDIF
         ELSE
          CALL RUNEND('ERROR: HARDEN. COEFF. CARD NOT FOUND')
         ENDIF
        ENDIF            ! iiaux.eq.1.or.iiaux.eq.2.or.iiaux.eq.4....
C
C**** LOOK FOR 'B_COEFFICIENT' CARD (B COEFFICIENT)
C
        IF(IIAUX.EQ.3.OR.IIAUX.EQ.6.OR.IIAUX.EQ.7.OR.IIAUX.EQ.8.OR.
     .     IIAUX.EQ.10.OR.IIAUX.EQ.11) THEN
         NPRIN=0
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'B_COE')THEN
          IF(WORDS(2).EQ.'FERRI') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20.OR.NLINE.LT.0)
     .      CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN B COEF.')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
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
           CALL RUNEND('ERROR: FERRITE CARD NOT FOUND')
          ENDIF
         ELSE
          CALL RUNEND('ERROR: B COEFFICIENT. CARD NOT FOUND')
         ENDIF
C
         NPRIN=0
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'B_COE')THEN
          IF(WORDS(2).EQ.'AUSTE') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20.OR.NLINE.LT.0)
     .      CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN B COEF.')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
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
           CALL RUNEND('ERROR: AUSTENITE CARD NOT FOUND')
          ENDIF
         ELSE
          CALL RUNEND('ERROR: B COEFFICIENT. CARD NOT FOUND')
         ENDIF
C
C**** LOOK FOR 'Q_COEFFICIENT' CARD (Q COEFFICIENT)
C
         NPRIN=0
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'Q_COE')THEN
          IF(WORDS(2).EQ.'FERRI') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20.OR.NLINE.LT.0)
     .      CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN Q COEF.')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
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
           CALL RUNEND('ERROR: FERRITE CARD NOT FOUND')
          ENDIF
         ELSE
          CALL RUNEND('ERROR: Q COEFFICIENT CARD NOT FOUND')
         ENDIF
C
         NPRIN=0
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'Q_COE')THEN
          IF(WORDS(2).EQ.'AUSTE') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20.OR.NLINE.LT.0)
     .      CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN Q COEF.')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
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
           CALL RUNEND('ERROR: AUSTENITE CARD NOT FOUND')
          ENDIF
         ELSE
          CALL RUNEND('ERROR: Q COEFFICIENT CARD NOT FOUND')
         ENDIF
        ENDIF       ! iiaux.eq.3.or.iiaux.eq.6....
       ELSE
        CALL RUNEND('ERROR: NO MODEL FOR ISOTROPIC HARDENING')
       ENDIF
      ELSE
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI2)=0.0
       GO TO 2
      ENDIF
C
C**** LOOK FOR 'KINEMATIC_HARDENING' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
    2 IF(WORDS(1).EQ.'KINEM')THEN
       KAUXX=1
       IF(KAUXX.GT.KPLA2) KPLA2=1
       IF(WORDS(2).EQ.'MODEL')THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        IIAUX=INT(PARAM(1))
        IF(IIAUX.LE.0.OR.IIAUX.GT.4)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF MODEL FOR KINE. HARDENING')
C
C**** LOOK FOR 'HD_COEFFICIENT' CARD (HARDENING COEFFICIENT)
C
        IF(IIAUX.EQ.1.OR.IIAUX.EQ.2) THEN
         NPRIN=0
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
         IF((WORDS(1).EQ.'HD_MO').OR.(WORDS(1).EQ.'HD_CO'))THEN
          IF(WORDS(2).EQ.'FERRI') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20.OR.NLINE.LT.0)
     .      CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN HARD. COEF.')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
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
           CALL RUNEND('ERROR: FERRITE CARD NOT FOUND')
          ENDIF
         ELSE
          CALL RUNEND('ERROR: HARDEN. COEFF. CARD NOT FOUND')
         ENDIF
C
         NPRIN=0
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
         IF((WORDS(1).EQ.'HD_MO').OR.(WORDS(1).EQ.'HD_CO'))THEN
          IF(WORDS(2).EQ.'AUSTE') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20.OR.NLINE.LT.0)
     .      CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN HARD. COEF.')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
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
           CALL RUNEND('ERROR: AUSTENITE CARD NOT FOUND')
          ENDIF
         ELSE
          CALL RUNEND('ERROR: HARDEN. COEFF. CARD NOT FOUND')
         ENDIF
        ENDIF            ! iiaux.eq.1.or.iiaux.eq.2
C
C**** LOOK FOR 'B_COEFFICIENT' CARD (B COEFFICIENT)
C
        IF(IIAUX.EQ.3.OR.IIAUX.EQ.4) THEN
         NPRIN=0
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'B_COE')THEN
          IF(WORDS(2).EQ.'FERRI') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20.OR.NLINE.LT.0)
     .      CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN B COEF.')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
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
           CALL RUNEND('ERROR: FERRITE CARD NOT FOUND')
          ENDIF
         ELSE
          CALL RUNEND('ERROR: B COEFFICIENT. CARD NOT FOUND')
         ENDIF
C
         NPRIN=0
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'B_COE')THEN
          IF(WORDS(2).EQ.'AUSTE') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20.OR.NLINE.LT.0)
     .      CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN B COEF.')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
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
           CALL RUNEND('ERROR: AUSTENITE CARD NOT FOUND')
          ENDIF
         ELSE
          CALL RUNEND('ERROR: B COEFFICIENT. CARD NOT FOUND')
         ENDIF
C
C**** LOOK FOR 'Q_COEFFICIENT' CARD (Q COEFFICIENT)
C
         NPRIN=0
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'Q_COE')THEN
          IF(WORDS(2).EQ.'FERRI') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20.OR.NLINE.LT.0)
     .      CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN Q COEF.')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
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
           CALL RUNEND('ERROR: FERRITE CARD NOT FOUND')
          ENDIF
         ELSE
          CALL RUNEND('ERROR: Q COEFFICIENT CARD NOT FOUND')
         ENDIF
C
         NPRIN=0
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'Q_COE')THEN
          IF(WORDS(2).EQ.'AUSTE') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20.OR.NLINE.LT.0)
     .      CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN Q COEF.')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
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
           CALL RUNEND('ERROR: AUSTENITE CARD NOT FOUND')
          ENDIF
         ELSE
          CALL RUNEND('ERROR: Q COEFFICIENT CARD NOT FOUND')
         ENDIF
        ENDIF       ! iiaux.eq.3
       ELSE
        CALL RUNEND('ERROR: NO MODEL FOR ISOTROPIC HARDENING')
       ENDIF
      ELSE
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI2)=0.0
       GO TO 3
      ENDIF
C
C**** LOOK FOR 'YIELD_FUNCTION' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
    3 IF(WORDS(1).EQ.'YIELD') THEN
       IF(WORDS(2).EQ.'VON_M') THEN        ! Von-Mises' type
        PROPS(36)=38.0                     ! SG cast iron yield function
        NCRIT=38
        IF(WORDS(3).EQ.'VERSI') THEN
         PROPS(35)=PARAM(1)
         IIAUX=INT(PROPS(35))
         IF(IIAUX.LE.0.OR.IIAUX.GT.3)
     .    CALL RUNEND('ERROR IN VERSION OF VON MISES YIELD CRIT.')
        ELSE
         PROPS(35)=1.0          ! default value for version of Von Mises
        ENDIF
       ELSE
        CALL RUNEND('PROPIVT1: ONLY V.MISES YIELD FUNC. IS AVAILABLE')
       ENDIF
      ELSE
       CALL RUNEND('PROPIVT1: YIELD FUNCTION CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'FLOW_POTENTIAL' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'FLOW_') THEN
       IF(WORDS(2).EQ.'VON_M') THEN        ! Von-Mises' type
        PROPS(52)=38.0                     ! SG cast iron flow potential
        NCRIP=38
        PROPS(56)=1.0D0         ! assumption
        KFLUJ=1
        IF(WORDS(3).EQ.'VERSI') THEN
         PROPS(51)=PARAM(1)
         IIAUX=INT(PROPS(51))
         IF(IIAUX.LE.0.OR.IIAUX.GT.3)
     .    CALL RUNEND('ERROR IN VERSION OF VON MISES FLOW POT.')
        ELSE
         PROPS(51)=1.0          ! default value for version of Von Mises
        ENDIF
       ELSE
        CALL RUNEND('PROPIVT1: ONLY V.MISES FLOW POTEN. IS AVAILABLE')
       ENDIF
      ELSE
       CALL RUNEND('PROPIVT1: FLOW_POTENTIAL CARD NOT FOUND')
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
C**** LOOK FOR 'PHASE_CHANGE' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
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
        WRITE(LURES,833) LINEA
        CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
C
        IF(WORDS(1).EQ.'SS_PH') THEN
C
         NPOI2=NPOI2+1
         PROPS(NPOI2)=0.0               ! ITEMRO
         NPOI2=NPOI2+1
         PROPS(NPOI2)=0.0
         NPOI2=NPOI2+1
         PROPS(NPOI2)=0.0
         NPOI2=NPOI2+1
         PROPS(NPOI2)=0.0
C
         NPOI2=NPOI2+1                  ! takes into account sense of pc
         PROPS(NPOI2)=0.0D0             ! ISEPC
         IF(WORDS(2).EQ.'HEATI') PROPS(NPOI2)=1.0D0
C
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'    ') THEN    ! constant expansion
          NPOI2=NPOI2+1
          PROPS(NPOI2)=0.0              ! IPCFOM
          NPOI2=NPOI2+1
          PROPS(NPOI2)=0.0              ! IPCMOM
          NPOI2=NPOI2+1
          PROPS(NPOI2)=0.0              ! ILSPC
          NPOI2=NPOI2+1
          PROPS(NPOI2)=LINEA            ! ISSPC
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(1)         ! expansion
          WRITE(LURES,836)
          WRITE(LURES,837) PROPS(NPOI2)
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
          PROPS(NPOI2)=1.0              ! ITEMRO
          ITEMRO=1
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(1)         ! TEMPS
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(2)         ! TEMPL
         ELSE
          NPOI2=NPOI2+1
          PROPS(NPOI2)=0.0              ! ITEMRO
          ITEMRO=0
          NPOI2=NPOI2+1
          PROPS(NPOI2)=0.0
          NPOI2=NPOI2+1
          PROPS(NPOI2)=0.0
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
         CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'    ') THEN    ! constant expansion
          NPOI2=NPOI2+1
          PROPS(NPOI2)=0.0              ! IPCFOM
          NPOI2=NPOI2+1
          PROPS(NPOI2)=0.0              ! IPCMOM
          NPOI2=NPOI2+1
          PROPS(NPOI2)=LINEA            ! ILSPC
          NPOI2=NPOI2+1
          PROPS(NPOI2)=0.0              ! ISSPC
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(1)         ! expansion
          WRITE(LURES,834)
          IF(ITEMRO.EQ.1) WRITE(LURES,838) PROPS(NPOI2-6),PROPS(NPOI2-5)
          WRITE(LURES,835) PROPS(NPOI2)
          WRITE(LURES,899)
         ENDIF
C
         IF(ISHRI.EQ.1) THEN
          KPLA5=1
          WRITE(LURES,824)
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
       PROPS(NPOI1)=0.0    ! nplatm
       GO TO 4
      ENDIF
C
C**** LOOK FOR FREE ENERGY MODEL
C
      CALL LISTEN('PROPIPT',NPRIN,ITAPE)
    4 IF(WORDS(1).EQ.'FREE_') THEN
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
       WRITE(LURES,825) IFREN
       WRITE(LURES,899)
      ELSE                         ! defaults
       IF(LARGE.EQ.0) THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=3.0           ! ifren
        IFREN=3
       ELSE
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=2.0           ! ifren
        IFREN=2
        IF(NITERC.EQ.4)
     .   CALL RUNEND('ERROR: IFREN=2 & NITERC=4 ARE INCOMPATIBLE')
       ENDIF
       WRITE(LURES,825) IFREN
       WRITE(LURES,899)
       GO TO 5
      ENDIF
C
C**** FATIGUE MODEL
C
      CALL LISTEN('PROPIPN',NPRIN,ITAPE)
    5 IF(WORDS(1).EQ.'FATIG') THEN
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
       GO TO 6
      ENDIF
C
C**** LOOK FOR 'END_MATERIAL' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIVT1',NPRIN,ITAPE)
    6 IF(WORDS(1).NE.'END_M')
     . CALL RUNEND('PROPIVT1: END_MATERIAL CARD NOT FOUND   ')
      WRITE(LURES,899)
C
C**** CONTROLS DIMENSION OF PROPS
C
      IF(NPOI2.GT.NPROP) THEN
       WRITE(LURES,900) NPROP,NPOI2
       CALL RUNEND('PROPIVT1: TOO MANY PROPERTIES TO READ')
      ENDIF
C
      RETURN
  801 FORMAT(/,'  INITIAL DENSITY VALUE=',E15.6,/)
  803 FORMAT(/,'  YOUNG MODULUS OF FERRITE     Vs.     TEMPERATURE',/)
  813 FORMAT(/,'  YOUNG MODULUS OF AUSTENITE   Vs.     TEMPERATURE',/)
  804 FORMAT(E15.6,21X,E15.6)
  805 FORMAT(/,'  POISSON RATIO OF FERRITE     Vs.     TEMPERATURE',/)
  815 FORMAT(/,'  POISSON RATIO OF AUSTENITE   Vs.     TEMPERATURE',/)
  806 FORMAT(/,'  THERMAL DILAT. OF FERRITE    Vs.     TEMPERATURE',/)
  816 FORMAT(/,'  THERMAL DILAT. OF AUSTENITE  Vs.     TEMPERATURE',/)
  807 FORMAT(/,'  THERMAL HARDENING OF FERRITE Vs.     TEMPERATURE',/)
  817 FORMAT(/,'  THERMAL HARDEN. OF AUSTENITE Vs.     TEMPERATURE',/)
  808 FORMAT(/,'  VISCOSITY FUNC. OF FERRITE   Vs.     TEMPERATURE',/)
  818 FORMAT(/,'  VISCOSITY FUNC. OF AUSTENITE Vs.     TEMPERATURE',/)
  809 FORMAT(/,'  EXPONENT FUNC. OF FERRITE    Vs.     TEMPERATURE',/)
  819 FORMAT(/,'  EXPONENT FUNC. OF AUSTENITE  Vs.     TEMPERATURE',/)
  810 FORMAT(/,'  HARDENING COEF. OF FERRITE   Vs.     TEMPERATURE',/)
  820 FORMAT(/,'  HARDENING COEF. OF ASUTENITE Vs.     TEMPERATURE',/)
  811 FORMAT(/,'  B COEFFICIENT OF FERRITE     Vs.     TEMPERATURE',/)
  821 FORMAT(/,'  B COEFFICIENT OF AUSTENITE   Vs.     TEMPERATURE',/)
  812 FORMAT(/,'  Q COEFFICIENT OF FERRITE     Vs.     TEMPERATURE',/)
  822 FORMAT(/,'  Q COEFFICIENT OF AUSTENITE   Vs.     TEMPERATURE',/)
  824 FORMAT(/,3X,'SHRINKAGE EFFECTS ARE CONSIDERED',/)
  825 FORMAT(/,3X,'FREE ENERGY MODEL=',I5,/)
  833 FORMAT(/,3X,'PHASE-CHANGE NUMBER  =',I5,/)
  834 FORMAT(/,3X,'LIQUID-SOLID PHASE-CHANGE',/)
  835 FORMAT(/,3X,'LIQUID-SOLID METALLURGICAL EXPANSION (%) = ',E15.6,/)
  836 FORMAT(/,3X,'SOLID-SOLID (EUTECTOID) PHASE-CHANGE',/)
  837 FORMAT(/,3X,'EUTECTOID METALLURGICAL EXPANSION (%) = ',E15.6,/)
  838 FORMAT(/,3X,'SOLIDUS TEMPERATURE  =',E15.6,/,
     .         3X,'LIQUIDUS TEMPERATURE =',E15.6,/)
  899 FORMAT(/)
  900 FORMAT(//,'TOO MANY MATERIAL PROPERTIES TO READ:',/,
     .      20X,'NPROP =',I5,5X,'NUMBE =',I5,/)
      END
