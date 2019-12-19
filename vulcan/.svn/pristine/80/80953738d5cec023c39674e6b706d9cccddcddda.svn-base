      SUBROUTINE PROPOPT(ITAPE,PROPS)
C***********************************************************************
C
C**** THIS ROUTINE READS THE MATERIAL PROPERTIES FOR THE ORTHOTROPIC,
C     PLASTIC & TEMPERATURE-DEPENDENT SOLID MODEL WITH DAMAGE
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
      PROPS(1)=1.0             ! orthotropic
      PROPS(2)=20.0            ! plastic
      PROPS(3)=2.0             ! temperature-dependent
      PROPS(59)=3.0            ! solidification
C
      PROPS(37)=0.0            ! non prescribed strains (see incstn.f)
      PROPS(19)=0.0
C
      PROPS(21)=0.0            ! no external damping
      PROPS(22)=0.0
C
      IF(NANIS.EQ.0)
     . CALL RUNEND('ERROR: ANISO. CARD MUST BE INPUT IN PROBLEM DATA')
      IANIS=1                  ! at least one material is orthotropic
C
      IF(NDIME.EQ.1)
     . CALL RUNEND('ERROR: 1D ANISOTROPY MODEL NOT IMPLEMENTED')
C
C**** READ & WRITE MATERIAL PROPERTIES
C
      NPRIN=0
C
C**** LOOK FOR 'DENSITY' CARD
C
      CALL LISTEN('PROPOPT',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'DENSI') THEN
       NLINE=INT(PARAM(1))
       IF(NLINE.GT.20)
     .  CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN DENSITY')
       NPRIN=1
       DO LINEA=1,NLINE
        CALL LISTEN('PROPOPT',NPRIN,ITAPE)
        PROPS(12)=PARAM(1)
       ENDDO
       WRITE(LURES,801) PROPS(12)
       WRITE(LURES,899)
      ELSE
       IF(KDYNA.EQ.1) THEN
        CALL RUNEND('PROPOPT: DENSITY CARD NOT FOUND')
       ELSE
        PROPS(12)=0.0D0
        CALL RUNMEN('WARNING: ZERO DENSITY VALUE IS ASSUMED')
        GO TO 1
       ENDIF
      ENDIF
C
C**** LOOK FOR 'YOUNG_MODULI' CARD
C
      NPRIN=0
      CALL LISTEN('PROPOPT',NPRIN,ITAPE)
    1 IF(WORDS(1).EQ.'YOUNG') THEN
       CALL LISTEN('PROPOPT',NPRIN,ITAPE)
       IF(WORDS(1).EQ.'COMPO'.AND.WORDS(2).EQ.'XX') THEN
        NPOI1=61
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN YOUNG MOD. XX')
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPOPT',NPRIN,ITAPE)
         DO I=1,2
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(I)
         ENDDO
        ENDDO
        WRITE(LURES,803)
        DO IPROP=NPOI1+1,NPOI2,2
         WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
        ENDDO
        WRITE(LURES,899)
       ELSE
        CALL RUNEND('PROPOPT: COMPONENT XX CARD NOT FOUND')
       ENDIF
       NPRIN=0
       CALL LISTEN('PROPOPT',NPRIN,ITAPE)
       IF(WORDS(1).EQ.'COMPO'.AND.WORDS(2).EQ.'YY') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN YOUNG MOD. YY')
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPOPT',NPRIN,ITAPE)
         DO I=1,2
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(I)
         ENDDO
        ENDDO
        WRITE(LURES,804)
        DO IPROP=NPOI1+1,NPOI2,2
         WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
        ENDDO
        WRITE(LURES,899)
       ELSE
        CALL RUNEND('PROPOPT: COMPONENT YY CARD NOT FOUND')
       ENDIF
       NPRIN=0
       CALL LISTEN('PROPOPT',NPRIN,ITAPE)
       IF(WORDS(1).EQ.'COMPO'.AND.WORDS(2).EQ.'ZZ') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN YOUNG MOD. ZZ')
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPOPT',NPRIN,ITAPE)
         DO I=1,2
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(I)
         ENDDO
        ENDDO
        WRITE(LURES,805)
        DO IPROP=NPOI1+1,NPOI2,2
         WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
        ENDDO
        WRITE(LURES,899)
       ELSE
        CALL RUNEND('PROPOPT: COMPONENT ZZ CARD NOT FOUND')
       ENDIF
      ELSE
       CALL RUNEND('PROPOPT: YOUNG_MODULI CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'POISSON_RATII' CARD
C
      NPRIN=0
      CALL LISTEN('PROPOPT',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'POISS') THEN
       CALL LISTEN('PROPOPT',NPRIN,ITAPE)
       IF(WORDS(1).EQ.'COMPO'.AND.WORDS(2).EQ.'XY') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN POISSON R. XY')
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPOPT',NPRIN,ITAPE)
         DO I=1,2
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(I)
         ENDDO
        ENDDO
        WRITE(LURES,806)
        DO IPROP=NPOI1+1,NPOI2,2
         WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
        ENDDO
        WRITE(LURES,899)
       ELSE
        CALL RUNEND('PROPOPT: COMPONENT XY CARD NOT FOUND')
       ENDIF
       NPRIN=0
       CALL LISTEN('PROPOPT',NPRIN,ITAPE)
       IF(WORDS(1).EQ.'COMPO'.AND.WORDS(2).EQ.'XZ') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN POISSON R. XZ')
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPOPT',NPRIN,ITAPE)
         DO I=1,2
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(I)
         ENDDO
        ENDDO
        WRITE(LURES,807)
        DO IPROP=NPOI1+1,NPOI2,2
         WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
        ENDDO
        WRITE(LURES,899)
       ELSE
        CALL RUNEND('PROPOPT: COMPONENT XZ CARD NOT FOUND')
       ENDIF
       NPRIN=0
       CALL LISTEN('PROPOPT',NPRIN,ITAPE)
       IF(WORDS(1).EQ.'COMPO'.AND.WORDS(2).EQ.'YZ') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN POISSON R. YZ')
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPOPT',NPRIN,ITAPE)
         DO I=1,2
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(I)
         ENDDO
        ENDDO
        WRITE(LURES,808)
        DO IPROP=NPOI1+1,NPOI2,2
         WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
        ENDDO
        WRITE(LURES,899)
       ELSE
        CALL RUNEND('PROPOPT: COMPONENT YZ CARD NOT FOUND')
       ENDIF
      ELSE
       CALL RUNEND('PROPOPT: POISSON_RATII CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'SHEAR_MODULI' CARD
C
      NPRIN=0
      CALL LISTEN('PROPOPT',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'SHEAR') THEN
       CALL LISTEN('PROPOPT',NPRIN,ITAPE)
       IF(WORDS(1).EQ.'COMPO'.AND.WORDS(2).EQ.'XY') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN SHEAR MOD. XY')
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPOPT',NPRIN,ITAPE)
         DO I=1,2
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(I)
         ENDDO
        ENDDO
        WRITE(LURES,809)
        DO IPROP=NPOI1+1,NPOI2,2
         WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
        ENDDO
        WRITE(LURES,899)
       ELSE
        CALL RUNEND('PROPOPT: COMPONENT XY CARD NOT FOUND')
       ENDIF
       IF(NDIME.EQ.3) THEN
        NPRIN=0
        CALL LISTEN('PROPOPT',NPRIN,ITAPE)
        IF(WORDS(1).EQ.'COMPO'.AND.WORDS(2).EQ.'XZ') THEN
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=PARAM(1)
         NLINE=INT(PARAM(1))
         IF(NLINE.GT.20)
     .    CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN SHEAR MOD. XZ')
         NPRIN=1
         DO LINEA=1,NLINE
          CALL LISTEN('PROPOPT',NPRIN,ITAPE)
          DO I=1,2
           NPOI2=NPOI2+1
           PROPS(NPOI2)=PARAM(I)
          ENDDO
         ENDDO
         WRITE(LURES,810)
         DO IPROP=NPOI1+1,NPOI2,2
          WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
         ENDDO
         WRITE(LURES,899)
        ELSE
         CALL RUNEND('PROPOPT: COMPONENT XZ CARD NOT FOUND')
        ENDIF
        NPRIN=0
        CALL LISTEN('PROPOPT',NPRIN,ITAPE)
        IF(WORDS(1).EQ.'COMPO'.AND.WORDS(2).EQ.'YZ') THEN
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=PARAM(1)
         NLINE=INT(PARAM(1))
         IF(NLINE.GT.20)
     .    CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN SHEAR MOD. YZ')
         NPRIN=1
         DO LINEA=1,NLINE
          CALL LISTEN('PROPOPT',NPRIN,ITAPE)
          DO I=1,2
           NPOI2=NPOI2+1
           PROPS(NPOI2)=PARAM(I)
          ENDDO
         ENDDO
         WRITE(LURES,811)
         DO IPROP=NPOI1+1,NPOI2,2
          WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
         ENDDO
         WRITE(LURES,899)
        ELSE
         CALL RUNEND('PROPOPT: COMPONENT YZ CARD NOT FOUND')
        ENDIF
       ENDIF
      ELSE
       CALL RUNEND('PROPOPT: SHEAR_MODULI CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'REFERENCE_TEMPERATURE' CARD
C
      NPRIN=0
      CALL LISTEN('PROPOPT',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'REFER') THEN
       PROPS(11)=PARAM(1)
      ELSE
       CALL RUNEND('PROPOPT: REFERENCE TEMPERATURE CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'THERMAL_DILATATION' CARD
C
      NPRIN=0
      CALL LISTEN('PROPOPT',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'THERM') THEN
C
       NPOI2=NPOI2+1
       PROPS(NPOI2)=0.0D0                          ! nalpu
       IF(WORDS(2).EQ.'TANGE') PROPS(NPOI2)=1.0D0  ! nalpu
C
       CALL LISTEN('PROPOPT',NPRIN,ITAPE)
       IF(WORDS(1).EQ.'COMPO'.AND.WORDS(2).EQ.'XX') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN THERMAL DIL. XX')
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPOPT',NPRIN,ITAPE)
         DO I=1,2
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(I)
         ENDDO
        ENDDO
        WRITE(LURES,812)
        DO IPROP=NPOI1+1,NPOI2,2
         WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
        ENDDO
        WRITE(LURES,899)
       ELSE
        CALL RUNEND('PROPOPT: COMPONENT XX CARD NOT FOUND')
       ENDIF
       NPRIN=0
       CALL LISTEN('PROPOPT',NPRIN,ITAPE)
       IF(WORDS(1).EQ.'COMPO'.AND.WORDS(2).EQ.'YY') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN THERMAL DIL. YY')
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPOPT',NPRIN,ITAPE)
         DO I=1,2
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(I)
         ENDDO
        ENDDO
        WRITE(LURES,813)
        DO IPROP=NPOI1+1,NPOI2,2
         WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
        ENDDO
        WRITE(LURES,899)
       ELSE
        CALL RUNEND('PROPOPT: COMPONENT YY CARD NOT FOUND')
       ENDIF
       NPRIN=0
       CALL LISTEN('PROPOPT',NPRIN,ITAPE)
       IF(WORDS(1).EQ.'COMPO'.AND.WORDS(2).EQ.'ZZ') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        NLINE=INT(PARAM(1))
        IF(NLINE.GT.20)
     .   CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN THERMAL DIL. ZZ')
        NPRIN=1
        DO LINEA=1,NLINE
         CALL LISTEN('PROPOPT',NPRIN,ITAPE)
         DO I=1,2
          NPOI2=NPOI2+1
          PROPS(NPOI2)=PARAM(I)
         ENDDO
        ENDDO
        WRITE(LURES,814)
        DO IPROP=NPOI1+1,NPOI2,2
         WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
        ENDDO
        WRITE(LURES,899)
       ELSE
        CALL RUNEND('PROPOPT: COMPONENT ZZ CARD NOT FOUND')
       ENDIF
      ELSE
       CALL RUNEND('PROPOPT: THERMAL_DILAT. CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'TH_HD_FUNCTION' CARD (THERMAL HARDENING)
C
C     Note: Elastic limit of the 0 direction
C
      NPRIN=0
      CALL LISTEN('PROPOPT',NPRIN,ITAPE)
      IF((WORDS(1).EQ.'HD_FU').OR.(WORDS(1).EQ.'TH_HD')) THEN
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI1)=PARAM(1)
       NLINE=INT(PARAM(1))
       IF(NLINE.GT.20)
     .  CALL RUNEND('ERROR: WRONG NUMB. OF POINTS IN TH. HD. FUNC.')
       NPRIN=1
       DO LINEA=1,NLINE
        CALL LISTEN('PROPOPT',NPRIN,ITAPE)
        DO I=1,2
         NPOI2=NPOI2+1
         PROPS(NPOI2)=PARAM(I)
        ENDDO
       ENDDO
       WRITE(LURES,815)
       DO IPROP=NPOI1+1,NPOI2,2
        WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
       ENDDO
       WRITE(LURES,899)
      ELSE
       CALL RUNEND('PROPOPT: THERMAL HARDEN. CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'ISOTROPIC_HARDENING' CARD
C
C     Notes: Hardening of the 0 direction
C            Only IFVER=1,4 (& consequently IPVER=1,4) are implemented;
C            see solidi.f
C
      ISOTT=0
      PROPS(35)=1.0D0       ! IFVER
      NPRIN=0
      CALL LISTEN('PROPOPT',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'ISOTR') THEN
       KAUXX=1
       IF(KAUXX.GT.KPLA1) KPLA1=1
       IF(WORDS(2).EQ.'MODEL') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        IIAUX=INT(PARAM(1))
        ISOTT=IIAUX
        IF(IIAUX.LE.0.OR.IIAUX.GT.13)
     .   CALL RUNEND('ERROR: WRONG NUMBER MODEL FOR ISOT. HARD.')
        IF(IIAUX.EQ.7.OR.IIAUX.EQ.9.OR.IIAUX.EQ.10.OR.
     .     IIAUX.EQ.11.OR.IIAUX.EQ.12) PROPS(35)=4.0D0
        IF(IIAUX.EQ.10.OR.IIAUX.EQ.11) THEN
         IF(KDYNA.EQ.0)
     .    CALL RUNMEN('WARNING: ISOT.HARD. MODEL=10 & STATIC ANAL.')
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
         CALL LISTEN('PROPOPT',NPRIN,ITAPE)
         IF((WORDS(1).EQ.'HD_MO').OR.(WORDS(1).EQ.'HD_CO')) THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          NLINE=INT(PARAM(1))
          IF(NLINE.GT.20.OR.NLINE.LT.0)
     .     CALL RUNEND('ERROR: WRONG NUMB.OF POINTS IN HARD. COEF.')
          NPRIN=1
          DO LINEA=1,NLINE
           CALL LISTEN('PROPOPT',NPRIN,ITAPE)
           DO I=1,2
            NPOI2=NPOI2+1
            PROPS(NPOI2)=PARAM(I)
           ENDDO
          ENDDO
          WRITE(LURES,816)
          DO IPROP=NPOI1+1,NPOI2,2
           WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
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
         CALL LISTEN('PROPOPT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'B_COE') THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          NLINE=INT(PARAM(1))
          IF(NLINE.GT.20.OR.NLINE.LT.0)
     .     CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN B COEF.')
          NPRIN=1
          DO LINEA=1,NLINE
           CALL LISTEN('PROPOPT',NPRIN,ITAPE)
           DO I=1,2
            NPOI2=NPOI2+1
            PROPS(NPOI2)=PARAM(I)
           ENDDO
          ENDDO
          WRITE(LURES,817)
          DO IPROP=NPOI1+1,NPOI2,2
           WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
          ENDDO
          WRITE(LURES,899)
         ELSE
          CALL RUNEND('ERROR: B COEFFICIENT. CARD NOT FOUND')
         ENDIF
C
C**** LOOK FOR 'Q_COEFFICIENT' CARD (Q COEFFICIENT)
C
         NPRIN=0
         CALL LISTEN('PROPOPT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'Q_COE') THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          NLINE=INT(PARAM(1))
          IF(NLINE.GT.20.OR.NLINE.LT.0)
     .     CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN Q COEF.')
          NPRIN=1
          DO LINEA=1,NLINE
           CALL LISTEN('PROPOPT',NPRIN,ITAPE)
           DO I=1,2
            NPOI2=NPOI2+1
            PROPS(NPOI2)=PARAM(I)
           ENDDO
          ENDDO
          WRITE(LURES,818)
          DO IPROP=NPOI1+1,NPOI2,2
           WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
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
      CALL LISTEN('PROPOPT',NPRIN,ITAPE)
    2 IF(WORDS(1).EQ.'KINEM') THEN

       call runend('error: kinematic hardening not implemented')

       KAUXX=1
       IF(KAUXX.GT.KPLA2) KPLA2=1
       IF(WORDS(2).EQ.'MODEL') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        IIAUX=INT(PARAM(1))
        IF(IIAUX.LE.0.OR.IIAUX.GT.4)
     .   CALL RUNEND('ERROR: WRONG NUMBER MODEL FOR KIN. HARDENING')
C
C**** LOOK FOR 'HD_COEFFICIENT' CARD (HARDENING COEFFICIENT)
C
        IF(IIAUX.EQ.1.OR.IIAUX.EQ.2) THEN
         NPRIN=0
         CALL LISTEN('PROPOPT',NPRIN,ITAPE)
         IF((WORDS(1).EQ.'HD_MO').OR.(WORDS(1).EQ.'HD_CO')) THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          NLINE=INT(PARAM(1))
          IF(NLINE.GT.20.OR.NLINE.LT.0)
     .     CALL RUNEND('ERROR: WRONG NUMB. OF POINTS IN HARD.COEF.')
          NPRIN=1
          DO LINEA=1,NLINE
           CALL LISTEN('PROPOPT',NPRIN,ITAPE)
           DO I=1,2
            NPOI2=NPOI2+1
            PROPS(NPOI2)=PARAM(I)
           ENDDO
          ENDDO
          WRITE(LURES,819)
          DO IPROP=NPOI1+1,NPOI2,2
           WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
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
         CALL LISTEN('PROPOPT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'B_COE') THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          NLINE=INT(PARAM(1))
          IF(NLINE.GT.20.OR.NLINE.LT.0)
     .     CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN B COEF.')
          NPRIN=1
          DO LINEA=1,NLINE
           CALL LISTEN('PROPOPT',NPRIN,ITAPE)
           DO I=1,2
            NPOI2=NPOI2+1
            PROPS(NPOI2)=PARAM(I)
           ENDDO
          ENDDO
          WRITE(LURES,820)
          DO IPROP=NPOI1+1,NPOI2,2
           WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
          ENDDO
          WRITE(LURES,899)
         ELSE
          CALL RUNEND('ERROR: B COEFFICIENT. CARD NOT FOUND')
         ENDIF
C
C**** LOOK FOR 'Q_COEFFICIENT' CARD (Q COEFFICIENT)
C
         NPRIN=0
         CALL LISTEN('PROPOPT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'Q_COE') THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          NLINE=INT(PARAM(1))
          IF(NLINE.GT.20.OR.NLINE.LT.0)
     .     CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN Q COEF.')
          NPRIN=1
          DO LINEA=1,NLINE
           CALL LISTEN('PROPOPT',NPRIN,ITAPE)
           DO I=1,2
            NPOI2=NPOI2+1
            PROPS(NPOI2)=PARAM(I)
           ENDDO
          ENDDO
          WRITE(LURES,821)
          DO IPROP=NPOI1+1,NPOI2,2
           WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
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
C**** LOOK FOR 'YIELD_FUNCTION' CARD
C
      NPRIN=0
      CALL LISTEN('PROPOPT',NPRIN,ITAPE)
    3 IF(WORDS(1).EQ.'YIELD') THEN
       IAUXX=0
       IF(WORDS(2).EQ.'HILL4') THEN             ! Hill 48
        IAUXX=1
        PROPS(36)=42.0D0                        ! NCRIT
        CALL LISTEN('PROPOPT',NPRIN,ITAPE)
        IAUXY=0
C
        IF(WORDS(1).EQ.'GENER') THEN            ! general case (2D & 3D)
         IF(NTYPE.EQ.3)
     .    CALL RUNEND('PROPOPT: NTYPE=3 & GENERAL CASE NOT IMPL. YET')
         IF(NTYPE.EQ.5)
     .    CALL RUNEND('PROPOPT: NTYPE=5 NOT POSSIBLE WITH HILL48 F-N')
         PROPS(45)=1.0D0                        ! NC42T
         IAUXY=1
C
C**** LOOK FOR F-N HILL PARAMETERS (F, G, H, N, M & L)
C
         NPRIN=0
         CALL LISTEN('PROPOPT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'HILL_') THEN
          CALL LISTEN('PROPOPT',NPRIN,ITAPE)
          IF(WORDS(1).EQ.'F_PAR') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20)
     .      CALL RUNEND('ERROR: WRONG NUMB. OF POINTS IN F_PARAMETER')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPOPT',NPRIN,ITAPE)
            DO I=1,2
             NPOI2=NPOI2+1
             PROPS(NPOI2)=PARAM(I)
            ENDDO
           ENDDO
           WRITE(LURES,1831)
           DO IPROP=NPOI1+1,NPOI2,2
            WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
           ENDDO
           WRITE(LURES,899)
          ELSE
           CALL RUNEND('PROPOPT: F_PARAMETER CARD NOT FOUND')
          ENDIF
          CALL LISTEN('PROPOPT',NPRIN,ITAPE)
          IF(WORDS(1).EQ.'G_PAR') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20)
     .      CALL RUNEND('ERROR: WRONG NUMB. OF POINTS IN G_PARAMETER')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPOPT',NPRIN,ITAPE)
            DO I=1,2
             NPOI2=NPOI2+1
             PROPS(NPOI2)=PARAM(I)
            ENDDO
           ENDDO
           WRITE(LURES,1832)
           DO IPROP=NPOI1+1,NPOI2,2
            WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
           ENDDO
           WRITE(LURES,899)
          ELSE
           CALL RUNEND('PROPOPT: G_PARAMETER CARD NOT FOUND')
          ENDIF
          CALL LISTEN('PROPOPT',NPRIN,ITAPE)
          IF(WORDS(1).EQ.'H_PAR') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20)
     .      CALL RUNEND('ERROR: WRONG NUMB. OF POINTS IN H_PARAMETER')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPOPT',NPRIN,ITAPE)
            DO I=1,2
             NPOI2=NPOI2+1
             PROPS(NPOI2)=PARAM(I)
            ENDDO
           ENDDO
           WRITE(LURES,1833)
           DO IPROP=NPOI1+1,NPOI2,2
            WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
           ENDDO
           WRITE(LURES,899)
          ELSE
           CALL RUNEND('PROPOPT: H_PARAMETER CARD NOT FOUND')
          ENDIF
          CALL LISTEN('PROPOPT',NPRIN,ITAPE)
          IF(WORDS(1).EQ.'N_PAR') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20)
     .      CALL RUNEND('ERROR: WRONG NUMB. OF POINTS IN N_PARAMETER')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPOPT',NPRIN,ITAPE)
            DO I=1,2
             NPOI2=NPOI2+1
             PROPS(NPOI2)=PARAM(I)
            ENDDO
           ENDDO
           WRITE(LURES,1834)
           DO IPROP=NPOI1+1,NPOI2,2
            WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
           ENDDO
           WRITE(LURES,899)
          ELSE
           CALL RUNEND('PROPOPT: N_PARAMETER CARD NOT FOUND')
          ENDIF
          IF(NTYPE.EQ.4) THEN                   ! 3D
           CALL LISTEN('PROPOPT',NPRIN,ITAPE)
           IF(WORDS(1).EQ.'M_PAR') THEN
            NPOI1=NPOI2+1
            NPOI2=NPOI1
            PROPS(NPOI1)=PARAM(1)
            NLINE=INT(PARAM(1))
            IF(NLINE.GT.20)
     .       CALL RUNEND('ERROR: WRONG NUMB. OF POINTS IN M_PARAMETER')
            NPRIN=1
            DO LINEA=1,NLINE
             CALL LISTEN('PROPOPT',NPRIN,ITAPE)
             DO I=1,2
              NPOI2=NPOI2+1
              PROPS(NPOI2)=PARAM(I)
             ENDDO
            ENDDO
            WRITE(LURES,1835)
            DO IPROP=NPOI1+1,NPOI2,2
             WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
            ENDDO
            WRITE(LURES,899)
           ELSE
            CALL RUNEND('PROPOPT: M_PARAMETER CARD NOT FOUND')
           ENDIF
           CALL LISTEN('PROPOPT',NPRIN,ITAPE)
           IF(WORDS(1).EQ.'L_PAR') THEN
            NPOI1=NPOI2+1
            NPOI2=NPOI1
            PROPS(NPOI1)=PARAM(1)
            NLINE=INT(PARAM(1))
            IF(NLINE.GT.20)
     .       CALL RUNEND('ERROR: WRONG NUMB. OF POINTS IN L_PARAMETER')
            NPRIN=1
            DO LINEA=1,NLINE
             CALL LISTEN('PROPOPT',NPRIN,ITAPE)
             DO I=1,2
              NPOI2=NPOI2+1
              PROPS(NPOI2)=PARAM(I)
             ENDDO
            ENDDO
            WRITE(LURES,1836)
            DO IPROP=NPOI1+1,NPOI2,2
             WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
            ENDDO
            WRITE(LURES,899)
           ELSE
            CALL RUNEND('PROPOPT: L_PARAMETER CARD NOT FOUND')
           ENDIF
          ENDIF
         ELSE
          CALL RUNEND('PROPOPT: HILL PARAMETERS CARD NOT FOUND')
         ENDIF
        ENDIF                                   ! general case (F-N)
C
        IF(WORDS(1).EQ.'PLANE') THEN            ! plane stress case
         IF(NTYPE.EQ.2.OR.NTYPE.EQ.5)
     .    CALL RUNEND('PROPOPT: NTYPE=2,5 NOT POSSIBLE WITH HILL48 PS')
C
         IF(WORDS(2).EQ.'STREN') THEN           ! strength based
          IF(NTYPE.EQ.3)
     .     CALL RUNEND('PROPOPT: NTYPE=3 & PL_ST-S NOT IMPL. YET')
          PROPS(45)=2.0D0                       ! NC42T
          IAUXY=1
C
C**** LOOK FOR 'STRENGTH_COEFFICIENTS' CARD
C
C     Note: 3 forms of considering the biaxial strength:
C           1) experimentally measured value
C           2) estimation (as an average; not consistent):
C              1/sigma_b^2=1/2(1/sigma_0^2+1/sigma_90^2)
C           3) in terms of Lankford coefficients:
C              sigma_b=(1/(1+R0)*R0/R90+1/(1+R0))^-1/2*sigma_0
C
          NPRIN=0
          CALL LISTEN('PROPOPT',NPRIN,ITAPE)
          IF(WORDS(1).EQ.'STREN') THEN
           CALL LISTEN('PROPOPT',NPRIN,ITAPE)
           IF(WORDS(1).EQ.'S_45_') THEN
            NPOI1=NPOI2+1
            NPOI2=NPOI1
            PROPS(NPOI1)=PARAM(1)
            NLINE=INT(PARAM(1))
            IF(NLINE.GT.20)
     .       CALL RUNEND('ERROR: WRONG NUMB. OF POINTS IN S_45 COEF.')
            NPRIN=1
            DO LINEA=1,NLINE
             CALL LISTEN('PROPOPT',NPRIN,ITAPE)
             DO I=1,2
              NPOI2=NPOI2+1
              PROPS(NPOI2)=PARAM(I)
             ENDDO
            ENDDO
            WRITE(LURES,1837)
            DO IPROP=NPOI1+1,NPOI2,2
             WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
            ENDDO
            WRITE(LURES,899)
           ELSE
            CALL RUNEND('PROPOPT: STRENGTH 45 COEFF. CARD NOT FOUND')
           ENDIF
           IF(NTYPE.NE.3) THEN        ! only one S is needed for axisym.
            CALL LISTEN('PROPOPT',NPRIN,ITAPE)
            IF(WORDS(1).EQ.'S_90_') THEN
             NPOI1=NPOI2+1
             NPOI2=NPOI1
             PROPS(NPOI1)=PARAM(1)
             NLINE=INT(PARAM(1))
             IF(NLINE.GT.20)
     .        CALL RUNEND('ERROR: WRONG NUMB. OF POINTS IN S_90 COEF.')
             NPRIN=1
             DO LINEA=1,NLINE
              CALL LISTEN('PROPOPT',NPRIN,ITAPE)
              DO I=1,2
               NPOI2=NPOI2+1
               PROPS(NPOI2)=PARAM(I)
              ENDDO
             ENDDO
             WRITE(LURES,1838)
             DO IPROP=NPOI1+1,NPOI2,2
              WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
             ENDDO
             WRITE(LURES,899)
            ELSE
             CALL RUNEND('PROPOPT: STRENGTH 90 COEFF. CARD NOT FOUND')
            ENDIF
            CALL LISTEN('PROPOPT',NPRIN,ITAPE)
            IF(WORDS(1).EQ.'S_BIA') THEN
             NPOI1=NPOI2+1
             NPOI2=NPOI1
             PROPS(NPOI1)=PARAM(1)
             NLINE=INT(PARAM(1))
             IF(NLINE.GT.20)
     .        CALL RUNEND('ERROR: WRONG NUMB. OF POINTS IN S_BIA COEF.')
             NPRIN=1
             DO LINEA=1,NLINE
              CALL LISTEN('PROPOPT',NPRIN,ITAPE)
              DO I=1,2
               NPOI2=NPOI2+1
               PROPS(NPOI2)=PARAM(I)
              ENDDO
             ENDDO
             WRITE(LURES,1839)
             DO IPROP=NPOI1+1,NPOI2,2
              WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
             ENDDO
             WRITE(LURES,899)
            ELSE
             CALL RUNEND('PROPOPT: STRENGTH BIA COEFF. CARD NOT FOUND')
            ENDIF
           ENDIF
          ELSE
           CALL RUNEND('PROPOPT: STRENGTH COEFF. CARD NOT FOUND')
          ENDIF
         ENDIF                                  ! strength based
C
         IF(WORDS(2).EQ.'LANKF') THEN           ! R based
          PROPS(45)=3.0D0                       ! NC42T
          IAUXY=1
C
C**** LOOK FOR 'LANKFORD_COEFFICIENTS' CARD
C
          NPRIN=0
          CALL LISTEN('PROPOPT',NPRIN,ITAPE)
          IF(WORDS(1).EQ.'LANKF') THEN
           CALL LISTEN('PROPOPT',NPRIN,ITAPE)
           IF(WORDS(1).EQ.'R_0_C') THEN
            NPOI1=NPOI2+1
            NPOI2=NPOI1
            PROPS(NPOI1)=PARAM(1)
            NLINE=INT(PARAM(1))
            IF(NLINE.GT.20)
     .       CALL RUNEND('ERROR: WRONG NUMB. OF POINTS IN R0 COEF.')
            NPRIN=1
            DO LINEA=1,NLINE
             CALL LISTEN('PROPOPT',NPRIN,ITAPE)
             DO I=1,2
              NPOI2=NPOI2+1
              PROPS(NPOI2)=PARAM(I)
             ENDDO
            ENDDO
            WRITE(LURES,831)
            DO IPROP=NPOI1+1,NPOI2,2
             WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
            ENDDO
            WRITE(LURES,899)
           ELSE
            CALL RUNEND('PROPOPT: R0 LANKFORD COEFF. CARD NOT FOUND')
           ENDIF
           IF(NTYPE.NE.3) THEN        ! only one R is needed for axisym.
            CALL LISTEN('PROPOPT',NPRIN,ITAPE)
            IF(WORDS(1).EQ.'R_45_') THEN
             NPOI1=NPOI2+1
             NPOI2=NPOI1
             PROPS(NPOI1)=PARAM(1)
             NLINE=INT(PARAM(1))
             IF(NLINE.GT.20)
     .        CALL RUNEND('ERROR: WRONG NUMB. OF POINTS IN R45 COEF.')
             NPRIN=1
             DO LINEA=1,NLINE
              CALL LISTEN('PROPOPT',NPRIN,ITAPE)
              DO I=1,2
               NPOI2=NPOI2+1
               PROPS(NPOI2)=PARAM(I)
              ENDDO
             ENDDO
             WRITE(LURES,832)
             DO IPROP=NPOI1+1,NPOI2,2
              WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
             ENDDO
             WRITE(LURES,899)
            ELSE
             CALL RUNEND('PROPOPT: R45 LANKFORD COEFF. CARD NOT FOUND')
            ENDIF
            CALL LISTEN('PROPOPT',NPRIN,ITAPE)
            IF(WORDS(1).EQ.'R_90_') THEN
             NPOI1=NPOI2+1
             NPOI2=NPOI1
             PROPS(NPOI1)=PARAM(1)
             NLINE=INT(PARAM(1))
             IF(NLINE.GT.20)
     .        CALL RUNEND('ERROR: WRONG NUMB. OF POINTS IN R90 COEF.')
             NPRIN=1
             DO LINEA=1,NLINE
              CALL LISTEN('PROPOPT',NPRIN,ITAPE)
              DO I=1,2
               NPOI2=NPOI2+1
               PROPS(NPOI2)=PARAM(I)
              ENDDO
             ENDDO
             WRITE(LURES,833)
             DO IPROP=NPOI1+1,NPOI2,2
              WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
             ENDDO
             WRITE(LURES,899)
            ELSE
             CALL RUNEND('PROPOPT: R90 LANKFORD COEFF. CARD NOT FOUND')
            ENDIF
           ENDIF
          ELSE
           CALL RUNEND('PROPOPT: LANKFORD COEFF. CARD NOT FOUND')
          ENDIF
         ENDIF                                  ! R based
        ENDIF                                   ! plane stress case
C
        IF(IAUXY.EQ.0)
     .   CALL RUNEND('PROPOPT: WRONG VERSION OF HILL 48 YIELD FUNCTION')
       ENDIF
       IF(IAUXX.EQ.0)
     .  CALL RUNEND('PROPOPT: WRONG YIELD FUNCTION         ')
       NCRIT=INT(PROPS(36))
      ELSE
       CALL RUNEND('PROPOPT: YIELD FUNCTION CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'FLOW_POTENTIAL' CARD
C
      PROPS(51)=PROPS(35)       ! IPVER (assumption!)
      NPRIN=0
      CALL LISTEN('PROPOPT',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'FLOW_') THEN
       IAUXX=0
       PROPS(47)=1.0D0          ! NASNA: assumed non-associate plastic.
       IF(WORDS(2).EQ.'ASSOC') THEN
        IAUXX=1
        PROPS(47)=0.0D0         ! NASNA: associate plasticity
        PROPS(52)=PROPS(36)     ! NCRIP
        PROPS(46)=PROPS(45)     ! NC42P
        PROPS(56)=1.0D0         ! assumed KFLUJ
        IF(WORDS(3).EQ.'INCRE') ! plastic flow computed incrementally
     .   PROPS(56)=0.0D0
       ENDIF
       IF(WORDS(2).EQ.'HILL4') THEN             ! Hill 48
        IAUXX=1
        PROPS(52)=42.0D0        ! NCRIP
        PROPS(56)=1.0D0         ! assumed KFLUJ
        IF(WORDS(3).EQ.'INCRE') ! plastic flow computed incrementally
     .   PROPS(56)=0.0D0
        CALL LISTEN('PROPOPT',NPRIN,ITAPE)
        IAUXY=0
C
        IF(WORDS(1).EQ.'GENER') THEN            ! general case (2D & 3D)
         IF(NTYPE.EQ.1.OR.NTYPE.EQ.5)
     .    CALL RUNEND('PROPOPT: NTYPE=1,5 NOT POSSIBLE WITH HILL48 3D')
         PROPS(46)=1.0D0                        ! NC42P
         IAUXY=1
C
C**** LOOK FOR F-N HILL PARAMETERS (F, G, H, N, M & L)
C
         NPRIN=0
         CALL LISTEN('PROPOPT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'HILL_') THEN
          CALL LISTEN('PROPOPT',NPRIN,ITAPE)
          IF(WORDS(1).EQ.'F_PAR') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20)
     .      CALL RUNEND('ERROR: WRONG NUMB. OF POINTS IN F_PARAMETER')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPOPT',NPRIN,ITAPE)
            DO I=1,2
             NPOI2=NPOI2+1
             PROPS(NPOI2)=PARAM(I)
            ENDDO
           ENDDO
           WRITE(LURES,2831)
           DO IPROP=NPOI1+1,NPOI2,2
            WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
           ENDDO
           WRITE(LURES,899)
          ELSE
           CALL RUNEND('PROPOPT: F_PARAMETER CARD NOT FOUND')
          ENDIF
          CALL LISTEN('PROPOPT',NPRIN,ITAPE)
          IF(WORDS(1).EQ.'G_PAR') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20)
     .      CALL RUNEND('ERROR: WRONG NUMB. OF POINTS IN G_PARAMETER')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPOPT',NPRIN,ITAPE)
            DO I=1,2
             NPOI2=NPOI2+1
             PROPS(NPOI2)=PARAM(I)
            ENDDO
           ENDDO
           WRITE(LURES,2832)
           DO IPROP=NPOI1+1,NPOI2,2
            WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
           ENDDO
           WRITE(LURES,899)
          ELSE
           CALL RUNEND('PROPOPT: G_PARAMETER CARD NOT FOUND')
          ENDIF
          CALL LISTEN('PROPOPT',NPRIN,ITAPE)
          IF(WORDS(1).EQ.'H_PAR') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20)
     .      CALL RUNEND('ERROR: WRONG NUMB. OF POINTS IN H_PARAMETER')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPOPT',NPRIN,ITAPE)
            DO I=1,2
             NPOI2=NPOI2+1
             PROPS(NPOI2)=PARAM(I)
            ENDDO
           ENDDO
           WRITE(LURES,2833)
           DO IPROP=NPOI1+1,NPOI2,2
            WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
           ENDDO
           WRITE(LURES,899)
          ELSE
           CALL RUNEND('PROPOPT: H_PARAMETER CARD NOT FOUND')
          ENDIF
          CALL LISTEN('PROPOPT',NPRIN,ITAPE)
          IF(WORDS(1).EQ.'N_PAR') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20)
     .      CALL RUNEND('ERROR: WRONG NUMB. OF POINTS IN N_PARAMETER')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPOPT',NPRIN,ITAPE)
            DO I=1,2
             NPOI2=NPOI2+1
             PROPS(NPOI2)=PARAM(I)
            ENDDO
           ENDDO
           WRITE(LURES,2834)
           DO IPROP=NPOI1+1,NPOI2,2
            WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
           ENDDO
           WRITE(LURES,899)
          ELSE
           CALL RUNEND('PROPOPT: N_PARAMETER CARD NOT FOUND')
          ENDIF
          IF(NTYPE.EQ.4) THEN                   ! 3D
           CALL LISTEN('PROPOPT',NPRIN,ITAPE)
           IF(WORDS(1).EQ.'M_PAR') THEN
            NPOI1=NPOI2+1
            NPOI2=NPOI1
            PROPS(NPOI1)=PARAM(1)
            NLINE=INT(PARAM(1))
            IF(NLINE.GT.20)
     .       CALL RUNEND('ERROR: WRONG NUMB. OF POINTS IN M_PARAMETER')
            NPRIN=1
            DO LINEA=1,NLINE
             CALL LISTEN('PROPOPT',NPRIN,ITAPE)
             DO I=1,2
              NPOI2=NPOI2+1
              PROPS(NPOI2)=PARAM(I)
             ENDDO
            ENDDO
            WRITE(LURES,2835)
            DO IPROP=NPOI1+1,NPOI2,2
             WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
            ENDDO
            WRITE(LURES,899)
           ELSE
            CALL RUNEND('PROPOPT: M_PARAMETER CARD NOT FOUND')
           ENDIF
           CALL LISTEN('PROPOPT',NPRIN,ITAPE)
           IF(WORDS(1).EQ.'L_PAR') THEN
            NPOI1=NPOI2+1
            NPOI2=NPOI1
            PROPS(NPOI1)=PARAM(1)
            NLINE=INT(PARAM(1))
            IF(NLINE.GT.20)
     .       CALL RUNEND('ERROR: WRONG NUMB. OF POINTS IN L_PARAMETER')
            NPRIN=1
            DO LINEA=1,NLINE
             CALL LISTEN('PROPOPT',NPRIN,ITAPE)
             DO I=1,2
              NPOI2=NPOI2+1
              PROPS(NPOI2)=PARAM(I)
             ENDDO
            ENDDO
            WRITE(LURES,2836)
            DO IPROP=NPOI1+1,NPOI2,2
             WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
            ENDDO
            WRITE(LURES,899)
           ELSE
            CALL RUNEND('PROPOPT: L_PARAMETER CARD NOT FOUND')
           ENDIF
          ENDIF
         ELSE
          CALL RUNEND('PROPOPT: HILL PARAMETERS CARD NOT FOUND')
         ENDIF
        ENDIF                                   ! general case (F-N)
C
        IF(WORDS(1).EQ.'PLANE') THEN            ! plane stress case
         IF(NTYPE.EQ.2.OR.NTYPE.EQ.5)
     .    CALL RUNEND('PROPOPT: NTYPE=2,5 NOT POSSIBLE WITH HILL48 PS')
         IF(WORDS(2).EQ.'STREN')
     .    CALL RUNEND('PROPOPT: STRENGTH BASED FLOW POT. NOT POSSIBLE')
C
         IF(WORDS(2).EQ.'LANKF') THEN           ! R based
          PROPS(46)=3.0D0                       ! NC42P
          IAUXY=1
C
C**** LOOK FOR 'LANKFORD_COEFFICIENTS' CARD
C
          NPRIN=0
          CALL LISTEN('PROPOPT',NPRIN,ITAPE)
          IF(WORDS(1).EQ.'LANKF') THEN
           CALL LISTEN('PROPOPT',NPRIN,ITAPE)
           IF(WORDS(1).EQ.'R_0_C') THEN
            NPOI1=NPOI2+1
            NPOI2=NPOI1
            PROPS(NPOI1)=PARAM(1)
            NLINE=INT(PARAM(1))
            IF(NLINE.GT.20)
     .       CALL RUNEND('ERROR: WRONG NUMB. OF POINTS IN R0 COEF.')
            NPRIN=1
            DO LINEA=1,NLINE
             CALL LISTEN('PROPOPT',NPRIN,ITAPE)
             DO I=1,2
              NPOI2=NPOI2+1
              PROPS(NPOI2)=PARAM(I)
             ENDDO
            ENDDO
            WRITE(LURES,834)
            DO IPROP=NPOI1+1,NPOI2,2
             WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
            ENDDO
            WRITE(LURES,899)
           ELSE
            CALL RUNEND('PROPOPT: R0 LANKFORD COEFF. CARD NOT FOUND')
           ENDIF
           IF(NTYPE.NE.3) THEN        ! only one R is needed for axisym.
            CALL LISTEN('PROPOPT',NPRIN,ITAPE)
            IF(WORDS(1).EQ.'R_45_') THEN
             NPOI1=NPOI2+1
             NPOI2=NPOI1
             PROPS(NPOI1)=PARAM(1)
             NLINE=INT(PARAM(1))
             IF(NLINE.GT.20)
     .        CALL RUNEND('ERROR: WRONG NUMB. OF POINTS IN R45 COEF.')
             NPRIN=1
             DO LINEA=1,NLINE
              CALL LISTEN('PROPOPT',NPRIN,ITAPE)
              DO I=1,2
               NPOI2=NPOI2+1
               PROPS(NPOI2)=PARAM(I)
              ENDDO
             ENDDO
             WRITE(LURES,835)
             DO IPROP=NPOI1+1,NPOI2,2
              WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
             ENDDO
             WRITE(LURES,899)
            ELSE
             CALL RUNEND('PROPOPT: R45 LANKFORD COEFF. CARD NOT FOUND')
            ENDIF
            CALL LISTEN('PROPOPT',NPRIN,ITAPE)
            IF(WORDS(1).EQ.'R_90_') THEN
             NPOI1=NPOI2+1
             NPOI2=NPOI1
             PROPS(NPOI1)=PARAM(1)
             NLINE=INT(PARAM(1))
             IF(NLINE.GT.20)
     .        CALL RUNEND('ERROR: WRONG NUMB. OF POINTS IN R90 COEF.')
             NPRIN=1
             DO LINEA=1,NLINE
              CALL LISTEN('PROPOPT',NPRIN,ITAPE)
              DO I=1,2
               NPOI2=NPOI2+1
               PROPS(NPOI2)=PARAM(I)
              ENDDO
             ENDDO
             WRITE(LURES,836)
             DO IPROP=NPOI1+1,NPOI2,2
              WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
             ENDDO
             WRITE(LURES,899)
            ELSE
             CALL RUNEND('PROPOPT: R90 LANKFORD COEFF. CARD NOT FOUND')
            ENDIF
           ENDIF
          ELSE
           CALL RUNEND('PROPOPT: LANKFORD COEFF. CARD NOT FOUND')
          ENDIF
         ENDIF                     ! R based
        ENDIF                      ! plane stress case
C
        IF(IAUXY.EQ.0)
     .   CALL RUNEND('PROPOPT: WRONG VERSION OF HILL 48 POT. FUNCTION')
       ENDIF
       IF(IAUXX.EQ.0)
     .  CALL RUNEND('PROPOPT: WRONG FLOW POTENTIAL FUNCTION')
       NCRIP=INT(PROPS(52))
       NASNA=INT(PROPS(47))
      ELSE
       CALL RUNEND('PROPOPT: FLOW POTENTIAL CARD NOT FOUND')
      ENDIF
C
C**** SOME CONTROLS
C
      IF(NASNA.EQ.0) THEN
       IF(KSYMM.EQ.0)
     .  CALL RUNMEN('WARNING: ASSOCIATE PLAST. & NON SYMMETRIC SOLVER')
      ELSE
       IF(KSYMM.EQ.1)
     .  CALL RUNMEN('WARNING: NON ASSOCIATE PLAST. & SYMMETRIC SOLVER')
      ENDIF
C
      IF(NCRIT.EQ.41.OR.NCRIP.EQ.41) THEN         ! Gurson
C
C**** LOOK FOR 'POROSITY' CARD
C
       NPRIN=0
       CALL LISTEN('PROPOPT',NPRIN,ITAPE)
       IF(WORDS(1).EQ.'POROS') THEN
        IF(WORDS(2).EQ.'MODEL') THEN
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
          CALL LISTEN('PROPOPT',NPRIN,ITAPE)
          IF(WORDS(1).EQ.'POROS') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20)
     .      CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN POROS. PAR.')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPOPT',NPRIN,ITAPE)
            DO I=1,2
             NPOI2=NPOI2+1
             PROPS(NPOI2)=PARAM(I)
            ENDDO
           ENDDO
           WRITE(LURES,841)
           DO IPROP=NPOI1+1,NPOI2,2
            WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
           ENDDO
           WRITE(LURES,899)
          ELSE
           CALL RUNEND('PROPOPT: POROSITY_PARAMETER CARD NOT FOUND')
          ENDIF
         ENDIF
C
C**** LOOK FOR 'Q1_PARAMETER' CARD
C
         IF(IIAUX.EQ.4) THEN
          NPRIN=0
          CALL LISTEN('PROPOPT',NPRIN,ITAPE)
          IF(WORDS(1).EQ.'Q1_PA') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20)
     .      CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN Q1_PAR.')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPOPT',NPRIN,ITAPE)
            DO I=1,2
             NPOI2=NPOI2+1
             PROPS(NPOI2)=PARAM(I)
            ENDDO
           ENDDO
           WRITE(LURES,842)
           DO IPROP=NPOI1+1,NPOI2,2
            WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
           ENDDO
           WRITE(LURES,899)
          ELSE
           CALL RUNEND('PROPOPT: Q1_PARAMETER CARD NOT FOUND')
          ENDIF
         ENDIF
C
C**** LOOK FOR 'Q2_PARAMETER' CARD
C
         IF(IIAUX.EQ.4) THEN
          NPRIN=0
          CALL LISTEN('PROPOPT',NPRIN,ITAPE)
          IF(WORDS(1).EQ.'Q2_PA') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           NLINE=INT(PARAM(1))
           IF(NLINE.GT.20)
     .      CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN Q1_PAR.')
           NPRIN=1
           DO LINEA=1,NLINE
            CALL LISTEN('PROPOPT',NPRIN,ITAPE)
            DO I=1,2
             NPOI2=NPOI2+1
             PROPS(NPOI2)=PARAM(I)
            ENDDO
           ENDDO
           WRITE(LURES,843)
           DO IPROP=NPOI1+1,NPOI2,2
            WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
           ENDDO
           WRITE(LURES,899)
          ELSE
           CALL RUNEND('PROPOPT: Q2_PARAMETER CARD NOT FOUND')
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
      CALL LISTEN('PROPOPT',NPRIN,ITAPE)
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
        WRITE(LURES,861) LINEA
        CALL LISTEN('PROPOPT',NPRIN,ITAPE)
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
         CALL LISTEN('PROPOPT',NPRIN,ITAPE)
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
          WRITE(LURES,864)
          WRITE(LURES,863) PROPS(NPOI2)
          WRITE(LURES,899)
         ENDIF
        ENDIF
C
        IF(WORDS(1).EQ.'LS_PH') THEN
         ILSX=ILSX+1
         IF(ILSX.GT.1)
     .    CALL RUNEND('PROPOPT: TWO OR MORE LS PCS ARE NOT IMPLEMENTED')
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
         CALL LISTEN('PROPOPT',NPRIN,ITAPE)
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
          WRITE(LURES,862)
          IF(ITEMRO.EQ.1) WRITE(LURES,865) PROPS(NPOI2-6),PROPS(NPOI2-5)
          WRITE(LURES,815) PROPS(NPOI2)
          WRITE(LURES,899)
         ENDIF
C
         IF(ISHRI.EQ.1) THEN
          KPLA5=1
          WRITE(LURES,866)
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
       GO TO 4
      ENDIF
C
C**** LOOK FOR 'DAMAGE_PARAMETERS' CARD
C
      NPRIN=0
      CALL LISTEN('PROPOPT',NPRIN,ITAPE)
    4 IF(WORDS(1).EQ.'DAMAG') THEN
       KAUXX=1
       IF(KAUXX.GT.KPLA3) KPLA3=1
       IF(WORDS(2).EQ.'MODEL') THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)      ! idamg
        IDAMG=INT(PROPS(NPOI1))
        IF(IDAMG.NE.1.AND.IDAMG.NE.2.AND.IDAMG.NE.3.AND.IDAMG.NE.4.AND.
     .     IDAMG.NE.5.AND.IDAMG.NE.21)
     .   CALL RUNEND('PROPOPT: WRONG DAMAGE MODEL')
       ELSE                        ! default: damage model=1
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=1.0D0         ! idamg
        IDAMG=1
       ENDIF
C
       IF(IDAMG.EQ.1.OR.IDAMG.EQ.2) THEN
        call runmen('WARNING: plastic & viscopl. models are different')
C
C**** LOOK FOR 'DAMAGE_PARAMETER_A' CARD
C
        NPRIN=0
        CALL LISTEN('PROPOPT',NPRIN,ITAPE)
        IF(WORDS(1).EQ.'DAMAG') THEN
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=PARAM(1)
         NLINE=INT(PARAM(1))
         IF(NLINE.GT.20)
     .    CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN A PARAMTER')
         NPRIN=1
         DO LINEA=1,NLINE
          CALL LISTEN('PROPOPT',NPRIN,ITAPE)
          DO I=1,2
           NPOI2=NPOI2+1
           PROPS(NPOI2)=PARAM(I)
          ENDDO
         ENDDO
         WRITE(LURES,871)
         DO IPROP=NPOI1+1,NPOI2,2
          WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
         ENDDO
         WRITE(LURES,899)
        ELSE
         CALL RUNEND('PROPOPT: A PARAMETER NOT FOUND')
        ENDIF
       ENDIF
C
       IF(IDAMG.EQ.3) THEN
C
C**** LOOK FOR 'COEFFICIENT_B9' CARD
C
        NPRIN=0
        CALL LISTEN('PROPOPT',NPRIN,ITAPE)
        IF(WORDS(1).EQ.'COEB9') THEN
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=PARAM(1)
         NLINE=INT(PARAM(1))
         IF(NLINE.GT.20)
     .    CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN COEF. B9')
         NPRIN=1
         DO LINEA=1,NLINE
          CALL LISTEN('PROPOPT',NPRIN,ITAPE)
          DO I=1,2
           NPOI2=NPOI2+1
           PROPS(NPOI2)=PARAM(I)
          ENDDO
         ENDDO
         WRITE(LURES,839)
         DO IPROP=NPOI1+1,NPOI2,2
          WRITE(LURES,802) PROPS(IPROP),PROPS(IPROP+1)
         ENDDO
         WRITE(LURES,899)
        ELSE
         CALL RUNEND('PROPOPT: COEFFICIENT B9 CARD NOT FOUND')
        ENDIF
       ENDIF                  ! idamg.eq.3
C
       IF(IDAMG.EQ.4.OR.IDAMG.EQ.5) THEN
        call runend('ERROR: idamg=4,5 not implemented yet in propipt.f')
       ENDIF
C
       IF(IDAMG.EQ.21) THEN
        NPRIN=0
        CALL LISTEN('PROPOPT',NPRIN,ITAPE)
        IF(WORDS(1).EQ.'GAMMA') THEN
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI2)=PARAM(1)     ! gamma-
         NPOI2=NPOI2+1
         PROPS(NPOI2)=PARAM(2)     ! gamma+
         WRITE(LURES,872) IDAMG
         WRITE(LURES,873) PROPS(NPOI2-1)
         WRITE(LURES,874) PROPS(NPOI2)
         WRITE(LURES,899)
        ELSE
         CALL RUNEND('PROPOPT: GAMMA PARAMETERS NOT FOUND')
        ENDIF
       ENDIF
      ELSE
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI1)=0.0D0          ! idamg
       GO TO 5
      ENDIF
C
C**** LOOK FOR FREE ENERGY MODEL
C
      CALL LISTEN('PROPOPT',NPRIN,ITAPE)
    5 IF(WORDS(1).EQ.'FREE_') THEN
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
       WRITE(LURES,881) IFREN
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
       WRITE(LURES,881) IFREN
       WRITE(LURES,899)
       GO TO 6
      ENDIF
C
C**** FATIGUE MODEL
C
      CALL LISTEN('PROPOPT',NPRIN,ITAPE)
    6 IF(WORDS(1).EQ.'FATIG') THEN
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
         CALL RUNEND('PROPOPT: IFATM GT 2')
        ENDIF
       ELSE
        CALL RUNEND('ERROR: FATIGUE_MODEL CARD NOT FOUND')
       ENDIF
      ELSE
       GO TO 7
      ENDIF
C
C**** LOOK FOR 'END_MATERIAL' CARD
C
      NPRIN=0
      CALL LISTEN('PROPOPT',NPRIN,ITAPE)
    7 IF(WORDS(1).NE.'END_M')
     . CALL RUNEND('PROPOPT: END_MATERIAL CARD NOT FOUND   ')
      WRITE(LURES,899)
C
C**** CONTROLS DIMENSION OF PROPS
C
      IF(NPOI2.GT.NPROP) THEN
       WRITE(LURES,900) NPROP,NPOI2
       CALL RUNEND('PROPOPT: TOO MANY PROPERTIES TO READ')
      ENDIF
C
      RETURN
  801 FORMAT(/,'  CONSTANT DENSITY VALUE=',E15.6,/)
  802 FORMAT(E15.6,10X,E15.6)
  803 FORMAT(/,'  YOUNG MODULUS XX     Vs.   TEMPERATURE',/)
  804 FORMAT(/,'  YOUNG MODULUS YY     Vs.   TEMPERATURE',/)
  805 FORMAT(/,'  YOUNG MODULUS ZZ     Vs.   TEMPERATURE',/)
  806 FORMAT(/,'  POISSON RATIO XY     Vs.   TEMPERATURE',/)
  807 FORMAT(/,'  POISSON RATIO XZ     Vs.   TEMPERATURE',/)
  808 FORMAT(/,'  POISSON RATIO YZ     Vs.   TEMPERATURE',/)
  809 FORMAT(/,'  SHEAR MODULUS XY     Vs.   TEMPERATURE',/)
  810 FORMAT(/,'  SHEAR MODULUS XZ     Vs.   TEMPERATURE',/)
  811 FORMAT(/,'  SHEAR MODULUS YZ     Vs.   TEMPERATURE',/)
  812 FORMAT(/,'  THERMAL DILAT. XX    Vs.   TEMPERATURE',/)
  813 FORMAT(/,'  THERMAL DILAT. YY    Vs.   TEMPERATURE',/)
  814 FORMAT(/,'  THERMAL DILAT. ZZ    Vs.   TEMPERATURE',/)
  815 FORMAT(/,'  THERMAL HARDENING    Vs.   TEMPERATURE',/)
  816 FORMAT(/,'  HARDENING COEF.      Vs.   TEMPERATURE',/)
  817 FORMAT(/,'  B COEFFICIENT        Vs.   TEMPERATURE',/)
  818 FORMAT(/,'  Q COEFFICIENT        Vs.   TEMPERATURE',/)
  819 FORMAT(/,'  K HARDENING COEF.    Vs.   TEMPERATURE',/)
  820 FORMAT(/,'  K B COEFFICIENT      Vs.   TEMPERATURE',/)
  821 FORMAT(/,'  K Q COEFFICIENT      Vs.   TEMPERATURE',/)
  831 FORMAT(/,'  R0 LANKFORD F-COEF.  Vs.   TEMPERATURE',/)
  832 FORMAT(/,'  R45 LANKFORD F-COEF. Vs.   TEMPERATURE',/)
  833 FORMAT(/,'  R90 LANKFORD F-COEF. Vs.   TEMPERATURE',/)
  834 FORMAT(/,'  R0 LANKFORD G-COEF.  Vs.   TEMPERATURE',/)
  835 FORMAT(/,'  R45 LANKFORD G-COEF. Vs.   TEMPERATURE',/)
  836 FORMAT(/,'  R90 LANKFORD G-COEF. Vs.   TEMPERATURE',/)
 1831 FORMAT(/,'  F F-PARAMETER        Vs.   TEMPERATURE',/)
 1832 FORMAT(/,'  G F-PARAMETER        Vs.   TEMPERATURE',/)
 1833 FORMAT(/,'  H F-PARAMETER        Vs.   TEMPERATURE',/)
 1834 FORMAT(/,'  N F-PARAMETER        Vs.   TEMPERATURE',/)
 1835 FORMAT(/,'  M F-PARAMETER        Vs.   TEMPERATURE',/)
 1836 FORMAT(/,'  L F-PARAMETER        Vs.   TEMPERATURE',/)
 1837 FORMAT(/,'  S45 STRENGTH F-COEF. Vs.   TEMPERATURE',/)
 1838 FORMAT(/,'  S90 STRENGTH F-COEF. Vs.   TEMPERATURE',/)
 1839 FORMAT(/,'  SB STRENGTH F-COEF.  Vs.   TEMPERATURE',/)
 2831 FORMAT(/,'  F G-PARAMETER        Vs.   TEMPERATURE',/)
 2832 FORMAT(/,'  G G-PARAMETER        Vs.   TEMPERATURE',/)
 2833 FORMAT(/,'  H G-PARAMETER        Vs.   TEMPERATURE',/)
 2834 FORMAT(/,'  N G-PARAMETER        Vs.   TEMPERATURE',/)
 2835 FORMAT(/,'  M G-PARAMETER        Vs.   TEMPERATURE',/)
 2836 FORMAT(/,'  L G-PARAMETER        Vs.   TEMPERATURE',/)
  839 FORMAT(/,'  COEFFICIENT B9       Vs.   TEMPERATURE',/)
  841 FORMAT(/,'  POROSITY PARAM.      Vs.   TEMPERATURE',/)
  842 FORMAT(/,'  Q1 PARAMATER         Vs.   TEMPERATURE',/)
  843 FORMAT(/,'  Q2 PARAMATER         Vs.   TEMPERATURE',/)
  851 FORMAT(/,'  AV. LANKFORD COE.    Vs.   TEMPERATURE',/)
  861 FORMAT(/,3X,'PHASE-CHANGE NUMBER  =',I5,/)
  862 FORMAT(/,3X,'LIQUID-SOLID PHASE-CHANGE',/)
  863 FORMAT(/,3X,'CONSTANT EXPANSION (PERCENTAGE) = ',E15.6,/)
  864 FORMAT(/,3X,'SOLID-SOLID PHASE-CHANGE',/)
  865 FORMAT(/,3X,'SOLIDUS TEMPERATURE  =',E15.6,/,
     .         3X,'LIQUIDUS TEMPERATURE =',E15.6,/)
  866 FORMAT(/,3X,'SHRINKAGE EFFECTS ARE CONSIDERED',/)
  871 FORMAT(/,3X,'A PARAMATER      Vs.     TEMPERATURE',/)
  872 FORMAT(/,3X,'PARAMETERS OF THE DAMAGE MODEL=',I5,/)
  873 FORMAT(/,3X,'COMPRESSION GAMMA PARAMETER    =',E15.6,/)
  874 FORMAT(/,3X,'TENSION GAMMA PARAMETER        =',E15.6,/)
  881 FORMAT(/,3X,'FREE ENERGY MODEL=',I5,/)
  899 FORMAT(/)
  900 FORMAT(//,'TOO MANY MATERIAL PROPERTIES TO READ:',/,
     .      20X,'NPROP =',I5,5X,'NUMBE =',I5,/)
      END
