      SUBROUTINE PROPIPN(ITAPE,PROPS)
C***********************************************************************
C
C**** THIS ROUTINE READS THE MATERIAL PROPERTIES FOR THE ISOTROPIC, 
C     PLASTIC & NON TEMPERATURE-DEPENDENT SOLID MODEL WITH DAMAGE
C     (ELEMENT NO. 30)
C
C**** OUTPUT (see pointe.f)
C
C     KPLA1: isotropic hardening
C     KPLA2: kinematic hardening
C     KPLA3: damage
C     KPLA4: porosity (Gurson model)
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
      PROPS(2)=20.0            ! plastic
      PROPS(3)=1.0             ! no temperature-dependent
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
C**** LOOK FOR 'DENSITY' CARD
C
      CALL LISTEN('PROPIPN',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'DENSI') THEN
       PROPS(12)=PARAM(1)
       WRITE(LURES,801) PROPS(12)
       WRITE(LURES,899)
      ELSE
       IF(KDYNA.EQ.1) THEN
        CALL RUNEND('PROPIPN: DENSITY CARD NOT FOUND')
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
      CALL LISTEN('PROPIPN',NPRIN,ITAPE)
    1 IF(WORDS(1).EQ.'YOUNG') THEN
       NPOI1=61
       NPOI2=NPOI1
       PROPS(NPOI1)=PARAM(1)
       WRITE(LURES,803) PROPS(NPOI1)
       WRITE(LURES,899)
      ELSE
       CALL RUNEND('PROPIPN: YOUNG_MODULUS CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'POISSON_RATIO' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIPN',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'POISS')THEN
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI1)=PARAM(1)
       WRITE(LURES,805) PROPS(NPOI1)
       WRITE(LURES,899)
      ELSE
       CALL RUNEND('PROPIPN: POISSON_RATIO CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'REFERENCE_TEMPERATURE' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIPN',NPRIN,ITAPE)
      NTEMP=0
      IF(WORDS(1).EQ.'REFER') THEN
       PROPS(11)=PARAM(1)
       WRITE(LURES,804) PROPS(11)
      ELSE
       NTEMP=1
       GO TO 2
      ENDIF
C
C**** LOOK FOR 'THERMAL_DILATATION' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIPN',NPRIN,ITAPE)
    2 IF(WORDS(1).EQ.'THERM')THEN
       IF(NTEMP.EQ.1)
     .  CALL RUNEND('PROPIPN: REFERENCE TEMPERATURE CARD NOT FOUND')
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI1)=PARAM(1)
       WRITE(LURES,806) PROPS(NPOI1)
       WRITE(LURES,899)
      ELSE
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI2)=0.0
       IF(NTEMP.EQ.0)
     .  CALL RUNMEN('WARNING: THERMAL_DILAT. CARD NOT FOUND')
       GO TO 3
      ENDIF
C
C**** LOOK FOR 'HD_FUNCTION' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIPN',NPRIN,ITAPE)
    3 IF((WORDS(1).EQ.'HD_FU').OR.(WORDS(1).EQ.'TH_HD'))THEN
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI1)=PARAM(1)
       WRITE(LURES,807) PROPS(NPOI1)
       WRITE(LURES,899)
      ELSE
       CALL RUNEND('PROPIPN: HARDEN. FUNCTION CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'ISOTROPIC_HARDENING' CARD
C
      NPRIN=0
      ISOTT=0
      CALL LISTEN('PROPIPN',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'ISOTR')THEN
       KAUXX=1
       IF(KAUXX.GT.KPLA1) KPLA1=1
       IF(WORDS(2).EQ.'MODEL')THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)
        IIAUX=INT(PARAM(1))
        ISOTT=IIAUX
        IF(IIAUX.LE.0.OR.IIAUX.GT.13)
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
         CALL LISTEN('PROPIPN',NPRIN,ITAPE)
         IF((WORDS(1).EQ.'HD_MO').OR.(WORDS(1).EQ.'HD_CO'))THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          WRITE(LURES,810) PROPS(NPOI1)
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
         CALL LISTEN('PROPIPN',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'B_COE')THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          WRITE(LURES,811) PROPS(NPOI1)
          WRITE(LURES,899)
         ELSE
          CALL RUNEND('ERROR: B COEFFICIENT. CARD NOT FOUND')
         ENDIF
C
C**** LOOK FOR 'Q_COEFFICIENT' CARD (Q COEFFICIENT)
C
         NPRIN=0
         CALL LISTEN('PROPIPN',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'Q_COE')THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          WRITE(LURES,812) PROPS(NPOI1)
          WRITE(LURES,899)
         ELSE
          CALL RUNEND('ERROR: Q COEFFICIENT CARD NOT FOUND')
         ENDIF
        ENDIF
C
C**** LOOK FOR 'HD_COEFFICIENT' CARD (Ep & T DEP. HARDENING COEFFICIENT)
C
        IF(IIAUX.EQ.9) THEN
         NPRIN=0
         CALL LISTEN('PROPIPN',NPRIN,ITAPE)
         IF((WORDS(1).EQ.'HD_MO').OR.(WORDS(1).EQ.'HD_CO'))THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          NLINE=INT(PARAM(1))
          IF(NLINE.GT.20)
     .     CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN HA. CO. FUNC.')
          NPRIN=1
          DO LINEA=1,NLINE
           CALL LISTEN('PROPIPT',NPRIN,ITAPE)
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
          CALL RUNEND('ERROR: HARDEN. COEFF. CARD NOT FOUND')
         ENDIF
        ENDIF
       ELSE
        CALL RUNEND('ERROR: NO MODEL FOR ISOTROPIC HARDENING')
       ENDIF
      ELSE
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI2)=0.0
       GO TO 4
      ENDIF
C
C**** LOOK FOR 'KINEMATIC_HARDENING' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIPT',NPRIN,ITAPE)
    4 IF(WORDS(1).EQ.'KINEM')THEN
       KAUXX=1
       IF(KAUXX.GT.KPLA2) KPLA2=1
       IF(WORDS(2).EQ.'MODEL')THEN
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
         CALL LISTEN('PROPIPT',NPRIN,ITAPE)
         IF((WORDS(1).EQ.'HD_MO').OR.(WORDS(1).EQ.'HD_CO'))THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          WRITE(LURES,810) PROPS(NPOI1)
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
         CALL LISTEN('PROPIPT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'B_COE')THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          WRITE(LURES,811) PROPS(NPOI1)
          WRITE(LURES,899)
         ELSE
          CALL RUNEND('ERROR: B COEFFICIENT. CARD NOT FOUND')
         ENDIF
C
C**** LOOK FOR 'Q_COEFFICIENT' CARD (Q COEFFICIENT)
C
         NPRIN=0
         CALL LISTEN('PROPIPT',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'Q_COE')THEN
          NPOI1=NPOI2+1
          NPOI2=NPOI1
          PROPS(NPOI1)=PARAM(1)
          WRITE(LURES,812) PROPS(NPOI1)
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
       PROPS(NPOI2)=0.0
       GO TO 5
      ENDIF
C
C**** LOOK FOR 'YIELD_FUNCTION' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIPN',NPRIN,ITAPE)
    5 IF(WORDS(1).EQ.'YIELD') THEN
       IAUXX=0
       IF(WORDS(2).EQ.'VON_M') THEN
        IAUXX=1
        PROPS(36)=32.0
        NCRIT=32
        IF(WORDS(3).EQ.'VERSI') THEN
         PROPS(35)=PARAM(1)
         IIAUX=INT(PROPS(35))
         IF(IIAUX.LE.0.OR.IIAUX.GT.4)
     .    CALL RUNEND('ERROR IN VERSION OF VON MISES YIELD CRIT.')
         IF(IIAUX.EQ.4.AND.NPRE3.EQ.0) INITP=1  ! ic index for int. var.
        ELSE
         PROPS(35)=1.0          ! default value for version of Von Mises
        ENDIF
       ENDIF
       IF(WORDS(2).EQ.'GURSO') THEN
        KAUXX=1
        IF(KAUXX.GT.KPLA4) KPLA4=1
        IAUXX=1
        PROPS(36)=41.0
        NCRIT=41
        IF(WORDS(3).EQ.'VERSI') THEN
         PROPS(35)=PARAM(1)
         IIAUX=INT(PROPS(35))
         IF(IIAUX.LE.0.OR.IIAUX.GT.4)
     .    CALL RUNEND('ERROR IN VERSION OF GURSON YIELD CRIT.')
        ELSE
         PROPS(35)=1.0          ! default value for version of Gurson
        ENDIF
       ENDIF
       IF(IAUXX.EQ.0)
     .  CALL RUNEND('PROPIPN: WRONG YIELD FUNCTION         ')
      ELSE
       CALL RUNEND('PROPIPN: YIELD FUNCTION CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'FLOW_POTENTIAL' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIPN',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'FLOW_') THEN
       IAUXX=0
       PROPS(47)=0.0D0          ! NASNA: assumed associate plastic.
       IF(WORDS(2).EQ.'VON_M') THEN
        IAUXX=1
        PROPS(52)=32.0
        NCRIP=32
        PROPS(56)=1.0D0         ! assumption
        KFLUJ=1
        IF(WORDS(3).EQ.'VERSI') THEN
         PROPS(51)=PARAM(1)
         IIAUX=INT(PROPS(51))
         IF(IIAUX.LE.0.OR.IIAUX.GT.4)
     .    CALL RUNEND('ERROR IN VERSION OF VON MISES FLOW POT.')
        ELSE
         PROPS(51)=1.0          ! default value for version of Von Mises
        ENDIF
       ENDIF
       IF(WORDS(2).EQ.'GURSO') THEN
        KAUXX=1
        IF(KAUXX.GT.KPLA4) KPLA4=1
        IAUXX=1
        PROPS(52)=41.0
        NCRIP=41
        IF(WORDS(3).EQ.'VERSI') THEN
         PROPS(51)=PARAM(1)
         IIAUX=INT(PROPS(51))
         IF(IIAUX.LE.0.OR.IIAUX.GT.4)
     .    CALL RUNEND('ERROR IN VERSION OF GURSON FLOW POT.')
        ELSE
         PROPS(51)=1.0          ! default value for version of Gurson
        ENDIF
       ENDIF
       IF(IAUXX.EQ.0)
     .  CALL RUNEND('PROPIPN: WRONG FLOW POTENTIAL FUNCTION         ')
      ELSE
       CALL RUNEND('PROPIPN: FLOW_POTENTIAL CARD NOT FOUND')
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
      IF(NCRIT.EQ.41.OR.NCRIP.EQ.41) THEN         ! Gurson
C
C**** LOOK FOR 'POROSITY' CARD
C
       NPRIN=0
       CALL LISTEN('PROPINT',NPRIN,ITAPE)
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
          CALL LISTEN('PROPINT',NPRIN,ITAPE)
          IF(WORDS(1).EQ.'POROS') THEN
           NPOI1=NPOI2+1
           NPOI2=NPOI1
           PROPS(NPOI1)=PARAM(1)
           WRITE(LURES,814) PROPS(NPOI1)
           WRITE(LURES,899)
          ELSE
           CALL RUNEND('PROPINT: POROSITY_PARAMETER CARD NOT FOUND')
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
           WRITE(LURES,815) PROPS(NPOI1)
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
           WRITE(LURES,816) PROPS(NPOI1)
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
C**** LOOK FOR 'LS_PHASE_CHANGE' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIPN',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'LS_PH') THEN
       CALL LISTEN('PROPIPN',NPRIN,ITAPE)
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI2)=PARAM(1)            ! solidus
       NPOI2=NPOI2+1
       PROPS(NPOI2)=PARAM(2)            ! liquidus
      ELSE
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI2)=10.0E+10
       NPOI2=NPOI2+1
       PROPS(NPOI2)=10.0E+10
       GO TO 6
      ENDIF
C
C**** LOOK FOR 'DAMAGE_PARAMETERS' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIPT',NPRIN,ITAPE)
    6 IF(WORDS(1).EQ.'DAMAG')THEN
       KAUXX=1
       IF(KAUXX.GT.KPLA3) KPLA3=1
       IF(WORDS(2).EQ.'MODEL')THEN
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=PARAM(1)      ! idamg
        IDAMG=INT(PROPS(NPOI1))
        IF(IDAMG.NE.1.AND.IDAMG.NE.2.AND.IDAMG.NE.3.AND.IDAMG.NE.4.AND.
     .     IDAMG.NE.5.AND.IDAMG.NE.6.AND.IDAMG.NE.21)
     .   CALL RUNEND('PROPIPT: WRONG DAMAGE MODEL')
       ELSE                        ! default: damage model=1
        NPOI1=NPOI2+1
        NPOI2=NPOI1
        PROPS(NPOI1)=1.0           ! idamg
        IDAMG=1
       ENDIF
C
       IF(IDAMG.EQ.1.OR.IDAMG.EQ.2) THEN
        call runmen('WARNING: plastic & viscopl.  models are different')
C
C**** LOOK FOR 'DAMAGE_PARAMETER_A' CARD
C
        NPRIN=0
        CALL LISTEN('PROPIPN',NPRIN,ITAPE)
        IF(WORDS(1).EQ.'DAMAG') THEN
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=PARAM(1)
         NPRIN=1
         WRITE(LURES,820) IDAMG
         WRITE(LURES,821) PROPS(NPOI1)
         WRITE(LURES,899)
        ELSE
         CALL RUNEND('PROPIPN: COEFFICIENT DAMAGE_P_A CARD NOT FOUND')
        ENDIF
       ENDIF                  ! idamg.eq.1.or.idamg.eq.2
C
       IF(IDAMG.EQ.3) THEN

        call runend('ERROR: idamg=3 not implemented yet in propipn.f')

C
C**** LOOK FOR 'COEFFICIENT_B9' CARD
C
        NPRIN=0
        CALL LISTEN('PROPIPT',NPRIN,ITAPE)
        IF(WORDS(1).EQ.'COEB9') THEN
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=PARAM(1)
         NLINE=INT(PARAM(1))
         IF(NLINE.GT.20)
     .    CALL RUNEND('ERROR: WRONG NUMBER OF POINTS IN COEF. B9')
         NPRIN=1
         DO LINEA=1,NLINE
          CALL LISTEN('PROPIPT',NPRIN,ITAPE)
          DO I=1,2
           NPOI2=NPOI2+1
           PROPS(NPOI2)=PARAM(I)
          ENDDO
         ENDDO
c        WRITE(LURES,839)
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
C
C**** LOOK FOR 'EPSILON_D_PARAMETER' CARD
C
        NPRIN=0
        CALL LISTEN('PROPIPN',NPRIN,ITAPE)
        IF(WORDS(1).EQ.'EPSIL') THEN
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=PARAM(1)
         NPRIN=1
         WRITE(LURES,820) IDAMG
         WRITE(LURES,822) PROPS(NPOI1)
         WRITE(LURES,899)
        ELSE
         CALL RUNEND('PROPIPN: COEFFICIENT EPSILON_D CARD NOT FOUND')
        ENDIF
C
C**** LOOK FOR 'EPSILON_R_PARAMETER' CARD
C
        NPRIN=0
        CALL LISTEN('PROPIPN',NPRIN,ITAPE)
        IF(WORDS(1).EQ.'EPSIL') THEN
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=PARAM(1)
         NPRIN=1
         WRITE(LURES,823) PROPS(NPOI1)
         WRITE(LURES,899)
        ELSE
         CALL RUNEND('PROPIPN: COEFFICIENT EPSILON_R CARD NOT FOUND')
        ENDIF
C
C**** LOOK FOR 'CRITICAL_DAMAGE_PARAMETER' CARD
C
        NPRIN=0
        CALL LISTEN('PROPIPN',NPRIN,ITAPE)
        IF(WORDS(1).EQ.'CRITI') THEN
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=PARAM(1)
         NPRIN=1
         WRITE(LURES,824) PROPS(NPOI1)
         WRITE(LURES,899)
        ELSE
         CALL RUNEND('PROPIPN: CRITICAL_DAMAGE CARD NOT FOUND')
        ENDIF
C
C**** LOOK FOR 'FINAL_DAMAGE_PARAMETER' CARD
C
        NPRIN=0
        CALL LISTEN('PROPIPN',NPRIN,ITAPE)
        IF(WORDS(1).EQ.'FINAL') THEN
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=PARAM(1)
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=PARAM(2)
         NPRIN=1
         WRITE(LURES,825) PROPS(NPOI1-1),PROPS(NPOI1)
         WRITE(LURES,899)
        ELSE
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=1.0D0   ! final damage (theoretical value)
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=0.0D0   ! smoothing exponent ("theoretical value")
         CALL RUNMEN('WARNING: FINAL_DAMAGE_VALUE SET TO 1.0')
         GO TO 7
        ENDIF
       ENDIF                  ! idamg.eq.4.or.idamg.eq.5
C
       IF(IDAMG.EQ.6) THEN
C
C**** LOOK FOR 'S_PARAMETER' CARD
C
        NPRIN=0
        CALL LISTEN('PROPIPN',NPRIN,ITAPE)
        IF(WORDS(1).EQ.'S_PAR') THEN
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=PARAM(1)
         NPRIN=1
         WRITE(LURES,820) IDAMG
         WRITE(LURES,822) PROPS(NPOI1)
         WRITE(LURES,899)
        ELSE
         CALL RUNEND('PROPIPN: COEFFICIENT S_PARAMETER CARD NOT FOUND')
        ENDIF
C
C**** LOOK FOR 'S_EXPONENT' CARD
C
        NPRIN=0
        CALL LISTEN('PROPIPN',NPRIN,ITAPE)
        IF(WORDS(1).EQ.'S_EXP') THEN
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=PARAM(1)
         NPRIN=1
         WRITE(LURES,829) PROPS(NPOI1)
         WRITE(LURES,899)
        ELSE
         CALL RUNEND('PROPIPN: COEFFICIENT S_EXPONENT CARD NOT FOUND')
        ENDIF
C
C**** LOOK FOR 'B_EXPONENT' CARD
C
        NPRIN=0
        CALL LISTEN('PROPIPN',NPRIN,ITAPE)
        IF(WORDS(1).EQ.'B_EXP') THEN
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=PARAM(1)
         NPRIN=1
         WRITE(LURES,830) PROPS(NPOI1)
         WRITE(LURES,899)
        ELSE
         CALL RUNEND('PROPIPN: CRITICAL_DAMAGE CARD NOT FOUND')
        ENDIF
C
C**** LOOK FOR 'A_PARAMETERS' CARD
C
        NPRIN=0
        CALL LISTEN('PROPIPN',NPRIN,ITAPE)
        IF(WORDS(1).EQ.'A_PAR') THEN
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=PARAM(1)
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=PARAM(2)
         NPRIN=1
         WRITE(LURES,831) PROPS(NPOI1-1),PROPS(NPOI1)
         WRITE(LURES,899)
        ELSE
         CALL RUNEND('PROPIPN: A_PARAMETERS CARD NOT FOUND')
        ENDIF
C
C**** LOOK FOR 'FINAL_DAMAGE_PARAMETER' CARD
C
        NPRIN=0
        CALL LISTEN('PROPIPN',NPRIN,ITAPE)
        IF(WORDS(1).EQ.'FINAL') THEN
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=PARAM(1)
         NPRIN=1
         WRITE(LURES,832) PROPS(NPOI1)
         WRITE(LURES,899)
        ELSE
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI1)=1.0D0   ! final damage (theoretical value)
         CALL RUNMEN('WARNING: FINAL_DAMAGE_VALUE SET TO 1.0')
         GO TO 7
        ENDIF
       ENDIF                  ! idamg.eq.6
C
       IF(IDAMG.EQ.21) THEN
        NPRIN=0
        CALL LISTEN('PROPIPT',NPRIN,ITAPE)
        IF(WORDS(1).EQ.'GAMMA')THEN
         NPOI1=NPOI2+1
         NPOI2=NPOI1
         PROPS(NPOI2)=PARAM(1)     ! gamma-
         NPOI2=NPOI2+1
         PROPS(NPOI2)=PARAM(2)     ! gamma+
         WRITE(LURES,820) IDAMG
         WRITE(LURES,826) PROPS(NPOI2-1)
         WRITE(LURES,827) PROPS(NPOI2)
         WRITE(LURES,899)
        ELSE
         CALL RUNEND('PROPIPT: GAMMA PARAMETERS NOT FOUND')
        ENDIF
       ENDIF
      ELSE
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI1)=0.0            ! idamg
       GO TO 7
      ENDIF
C
C**** LOOK FOR FREE ENERGY MODEL
C
      CALL LISTEN('PROPIPT',NPRIN,ITAPE)
    7 IF(WORDS(1).EQ.'FREE_') THEN
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
       ENDIF
       WRITE(LURES,817) IFREN
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
       ENDIF
       WRITE(LURES,817) IFREN
       WRITE(LURES,899)
       GO TO 8
      ENDIF
C
C**** FATIGUE MODEL
C
      CALL LISTEN('PROPIPN',NPRIN,ITAPE)
    8 IF(WORDS(1).EQ.'FATIG') THEN
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
       GO TO 9
      ENDIF
C
C**** LOOK FOR 'END_MATERIAL' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIPN',NPRIN,ITAPE)
    9 IF(WORDS(1).NE.'END_M')
     . CALL RUNEND('PROPIPN: END_MATERIAL CARD NOT FOUND   ')
      WRITE(LURES,899)
C
C**** CONTROLS DIMENSION OF PROPS
C
      IF(NPOI2.GT.NPROP) THEN
       WRITE(LURES,900) NPROP,NPOI2
       CALL RUNEND('PROPIPN: TOO MANY PROPERTIES TO READ')
      ENDIF
C
      RETURN
  801 FORMAT(/,'  CONSTANT DENSITY VALUE=',E15.6,/)
  802 FORMAT(E15.6,10X,E15.6)
  803 FORMAT(/,'  YOUNG MODULUS         =',E15.6,/)
  804 FORMAT(/,'  REFERENCE TEMPERATURE =',E15.6,/)
  805 FORMAT(/,'  POISSON RATIO         =',E15.6,/)
  806 FORMAT(/,'  THERMAL DILATATION    =',E15.6,/)
  807 FORMAT(/,'  HARDENING FUNCTION    =',E15.6,/)
  810 FORMAT(/,'  HARDENING COEFFICIENT =',E15.6,/)
  811 FORMAT(/,'  B COEFFICIENT         =',E15.6,/)
  812 FORMAT(/,'  Q COEFFICIENT         =',E15.6,/)
  813 FORMAT(/,'  HARDENING COEF.   Vs.     EFF. PL. DEF.',/)
  814 FORMAT(/,'  POROSITY PARAMETER    =',E15.6,/)
  815 FORMAT(/,'  Q1 PARAMATER          =',E15.6,/)
  816 FORMAT(/,'  Q2 PARAMATER          =',E15.6,/)
  817 FORMAT(/,3X,'FREE ENERGY MODEL=',I5,/)
  820 FORMAT(/,3X,'DAMAGE PARAMETERS FOR MODEL=',I5,/)
  821 FORMAT(/,'  DAMAGE PARAMETER A    =',E15.6,/)
  822 FORMAT(/,'  EPSILON_D PARAMETER   =',E15.6,/)
  823 FORMAT(/,'  EPSILON_R PARAMETER   =',E15.6,/)
  824 FORMAT(/,'  CRITICAL DAMAGE       =',E15.6,/)
  825 FORMAT(/,'  FINAL DAMAGE          =',E15.6,/,
     .       /,'  SMOOTHING EXPONENT    =',E15.6,/)
  826 FORMAT(/,'  COMPRESSION GAMMA PAR.=',E15.6,/)
  827 FORMAT(/,'  TENSION GAMMA PARAM.  =',E15.6,/)
  828 FORMAT(/,'  S PARAMETER           =',E15.6,/)
  829 FORMAT(/,'  S EXPONENT            =',E15.6,/)
  830 FORMAT(/,'  B EXPONENT            =',E15.6,/)
  831 FORMAT(/,'  A PARAMETERS          =',E15.6,2X,E15.6,/)
  832 FORMAT(/,'  FINAL DAMAGE          =',E15.6,/)
  899 FORMAT(/)
  900 FORMAT(//,'TOO MANY MATERIAL PROPERTIES TO READ:',/,
     .      20X,'NPROP =',I5,5X,'NUMBE =',I5,/)
      END
