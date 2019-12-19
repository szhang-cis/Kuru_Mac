      SUBROUTINE PROP04(ITAPE,PROPS)
C***********************************************************************
C
C**** THIS ROUTINE READS THE MATERIAL PROPERTIES FOR CONTACT MODEL
C     (ELEMENT NO. 4)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'auxl_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'inpo_om.f'
C
      DIMENSION PROPS(*)
C
C**** READ & WRITE MATERIAL PROPERTIES
C
      NPRIN=0
      PROPS(1)=1.0
      ICOMO=0    ! regularization model flag (only for penalty approach)
C
C**** LOOK FOR 'NORMAL_STIFFNESS' CARD
C
      CALL LISTEN('PROP04',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'NORMA') THEN
       IF(WORDS(2).EQ.'MODEL') THEN
        PROPS(7)=PARAM(1)
        ICOMO=INT(PROPS(7))
        IF(NPOIC.GT.0)
     .   CALL RUNEND('ERROR: CONTACT MODEL NOT POSSIBLE WITH NPOIC>0')
        CALL LISTEN('PROP04',NPRIN,ITAPE)
        IF(ICOMO.EQ.1) THEN
         PROPS(2)=PARAM(1)                  ! En
        ENDIF
        IF(ICOMO.EQ.2) THEN
         PROPS(2)=PARAM(1)                  ! En
         PROPS(8)=PARAM(2)                  ! \bar g
        ENDIF
        IF(ICOMO.EQ.3) THEN
         PROPS(2)=PARAM(1)                  ! En
         PROPS(9)=PARAM(2)                  ! En0
        ENDIF
        IF(ICOMO.EQ.4) THEN
         PROPS(2)=PARAM(1)                  ! En
         PROPS(8)=PARAM(2)                  ! \bar g
         PROPS(9)=PARAM(3)                  ! En0
        ENDIF
        IF(ICOMO.EQ.5) THEN                 ! mode I debounding model
         IF(WORDS(1).EQ.'COMPR') THEN
          PROPS(2)=PARAM(1)                 ! E_n_compression
         ELSE
          CALL RUNEND('ERROR: NO COMPRESSIVE STIFFNESS CARD FOUND')
         ENDIF
         CALL LISTEN('PROP04',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'MAXIM') THEN
          PROPS(17)=PARAM(1)                ! sigma_max
         ELSE
          CALL RUNEND('ERROR: NO MAX. NORM. CONTACT STRESS CARD FOUND')
         ENDIF
         CALL LISTEN('PROP04',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'NORMA') THEN
          PROPS(8)=PARAM(1)                 ! Gnc
         ELSE
          CALL RUNEND('ERROR: NO NORM.CRIT. FRACTURE ENERGY CARD FOUND')
         ENDIF
         CALL LISTEN('PROP04',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'PERCE') THEN
          PROPS(9)=PARAM(1)                 ! u_n-bar/u^c_n
         ELSE
          CALL RUNEND('ERROR: NO PERCENT. GAP AT DEBOUNDING CARD FOUND')
         ENDIF
         CALL LISTEN('PROP04',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'CRITI') THEN
          PROPS(18)=PARAM(1)                ! critical dn
         ELSE
          CALL RUNEND('ERROR: NO CRITICAL DEBOUNDING CARD FOUND')
         ENDIF
         CALL LISTEN('PROP04',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'TANGE') THEN
          PROPS(54)=PARAM(1)                ! E_t
         ELSE
          CALL RUNEND('ERROR: NO TANGENTIAL STIFFNESS CARD FOUND')
         ENDIF
        ENDIF
       ELSE
        PROPS(2)=PARAM(1)                   ! En
        IF(NPOIC.GT.0) THEN
         IF(IAUGM.EQ.0.OR.IAUGM.EQ.1) THEN
          PROPS(3)=PARAM(2)
          IMODEL=INT(PROPS(3))              ! contact model
          PROPS(4)=PARAM(3)
          IMODEX=INT(PROPS(4))              ! regularization model
          IF(IMODEL.EQ.0) THEN
           IF(KSOLV.EQ.0.OR.KSOLV.EQ.1)
     .      CALL RUNMEN('WARNING: NEGATIVE PIVOTS FOR MODEL 0 OF ALM')
          ELSE
           IF(KSOLV.EQ.0.OR.KSOLV.EQ.1) THEN
            IF(KSYMM.EQ.1)
     .       CALL RUNMEN('WARNING: SYMMETRIC SOLVER FOR MODEL 1 OF ALM')
           ENDIF
           IF(KSOLV.EQ.2)
     .      CALL RUNMEN('WARNING: SOLVER NOT ADEQUATE FOR MODEL 1 ALM')
          ENDIF
         ENDIF
         IF(IAUGM.EQ.4.OR.IAUGM.EQ.5) THEN
          PROPS(3)=1.0                      ! contact model
          PROPS(4)=1.0                      ! regularization model
          IF(KSOLV.EQ.0.OR.KSOLV.EQ.1) THEN
           IF(KSYMM.EQ.1)
     .      CALL RUNMEN('ERROR: SYMMETRIC SOLVER FOR MODEL 4,5 OF ALM')
          ENDIF
          IF(KSOLV.EQ.2)
     .     CALL RUNEND('ERROR: SOLVER NOT ADEQUATE FOR MODEL 1 OF ALM')
         ENDIF
        ENDIF
       ENDIF
      ELSE
       CALL RUNEND('PROP04: NORMAL_STIFFNESS CARD NOT FOUND')
      ENDIF
C
C**** READ CONTACT TEMPERATURE 'TEMPA' OR 'TEMPB' CARD
C
      CALL LISTEN('PROP04',NPRIN,ITAPE)
      ITEAB=0
      IF(WORDS(1).EQ.'TEMPA') THEN
       ITEAB=1
       PROPS(5)=PARAM(1)
       PROPS(6)=1.0      ! inote
      ELSE
       IF(WORDS(1).EQ.'TEMPB') THEN
        ITEAB=1
        PROPS(5)=PARAM(1)
        PROPS(6)=2.0     ! inote
       ELSE
        PROPS(5)=10.0E+10
        PROPS(6)=1.0     ! inote
        GO TO 1
       ENDIF
      ENDIF
C
C**** LOOK FOR 'JACOBIAN_REGULARIZATION_FACTOR' CARD
C
      CALL LISTEN('PROP04',NPRIN,ITAPE)
    1 IF(WORDS(1).EQ.'JACOB') THEN
       PROPS(10)=1.0
       PROPS(11)=PARAM(1)           ! factor 1
       PROPS(12)=PARAM(2)           ! factor 2
       PROPS(13)=PARAM(3)           ! nnn20
       NNN20=INT(PROPS(13))
       IF(NNN20.LE.0) PROPS(13)=20.0
      ELSE
       GO TO 2
      ENDIF
C
C**** LOOK FOR 'MINIMUM_&_MAXIMUM_ADMISSIBLE_CONTACT_GAPS' CARD
C
      CALL LISTEN('PROP04',NPRIN,ITAPE)
    2 IF(WORDS(1).EQ.'MINIM') THEN
       PROPS(14)=1.0
       PROPS(15)=PARAM(1)      ! minimum admissible contact gap (TOLGA)
       PROPS(16)=PARAM(2)      ! maximum admissible contact gap (TOLGAM)
      ELSE
       GO TO 3
      ENDIF
C
C**** LOOK FOR 'FRICTION' CARD
C
      CALL LISTEN('PROP04',NPRIN,ITAPE)
    3 IF(WORDS(1).EQ.'FRICT') THEN
       PROPS(51)=1.0D0
       IFRIC=INT(PROPS(51))
       IFRIM=0
       IF(WORDS(2).EQ.'COULO') THEN
        PROPS(52)=1.0D0
        IFRIM=INT(PROPS(52))
        PROPS(53)=PARAM(1)
        IFRCO=INT(PROPS(53))
        IF(IFRCO.EQ.1) THEN
         CALL LISTEN('PROP04',NPRIN,ITAPE)
         IF(WORDS(1).EQ.'TANGE') THEN
          PROPS(54)=PARAM(1)
          CALL LISTEN('PROP04',NPRIN,ITAPE)
          IF(WORDS(1).EQ.'FRICT') THEN
           PROPS(55)=PARAM(1)
          ELSE
           CALL RUNEND('PROP04: FRICTIONAL COEFFICIENT CARD NOT FOUND')
          ENDIF
         ELSE
          CALL RUNEND('PROP04: TANGENTIAL STIFFNESS CARD NOT FOUND')
         ENDIF
         IF((NDIME-1+NDIME*NDIME).GT.NHIST)    ! p_t & C_f; see frin32.f
     .    CALL RUNEND('ERROR: NDIME-1+NDIME*NDIME > NHIST (prop04)')
        ELSE
         CALL RUNEND('ERROR: ONLY COULOMB LAW No. 1 IS IMPLEM. NOW')
        ENDIF
       ELSE
        CALL RUNEND('ERROR: ONLY COULOMB LAW IS IMPLEM. NOW')
       ENDIF
      ELSE
       GO TO 4
      ENDIF
C
C**** LOOK FOR 'END_MATERIAL' CARD
C
      NPRIN=0
      CALL LISTEN('PROP04',NPRIN,ITAPE)
    4 IF(WORDS(1).NE.'END_M')
     . CALL RUNEND('PROP04: END_MATERIAL CARD NOT FOUND   ')
      WRITE(LURES,802)
C
C**** WRITES PROPERTIES
C
      IF(ICOMO.EQ.0) THEN
       WRITE(LURES,810)
       WRITE(LURES,803) PROPS(2)
      ENDIF
      IF(ICOMO.EQ.1) THEN
       WRITE(LURES,811)
       WRITE(LURES,803) PROPS(2)
      ENDIF
      IF(ICOMO.EQ.2) THEN
       WRITE(LURES,812)
       WRITE(LURES,803) PROPS(2),PROPS(8)
      ENDIF
      IF(ICOMO.EQ.3) THEN
       WRITE(LURES,813)
       WRITE(LURES,803) PROPS(2),PROPS(9)
      ENDIF
      IF(ICOMO.EQ.4) THEN
       WRITE(LURES,814)
       WRITE(LURES,803) PROPS(2),PROPS(8),PROPS(9)
      ENDIF
      IF(ICOMO.EQ.5) THEN
       WRITE(LURES,815)
       RIGIN=PROPS(17)*PROPS(17)/(2.0D0*PROPS(9)*PROPS(8))   ! K_n
       GAPCN=2.0D0*PROPS(8)/PROPS(17)                        ! u^c_n
       WRITE(LURES,816) PROPS(2),PROPS(17),PROPS(8),PROPS(9),
     .                  RIGIN,GAPCN,PROPS(18),PROPS(54)
       IF(INT(PROPS(51)).EQ.1)
     .  CALL RUNMEN('WARNING: DEBOUNDING & FRICTION ARE INCOMPATIBLE')
       PROPS(51)=1.0                         ! ifric=1
       PROPS(52)=2.0                         ! ifrim=2 (stick condition)
       CALL RUNMEN('WARNING: NO STABILIZATION IS USED')
       PROPS(5)=-10.0E+10                    ! negative infinite temper.
      ENDIF
C
      IF(ITEAB.EQ.1) THEN
       WRITE(LURES,820)
       WRITE(LURES,803) PROPS(5)
      ENDIF
C
      III10=INT(PROPS(10))
      IF(III10.EQ.0) THEN
       WRITE(LURES,830)
      ELSE
       WRITE(LURES,831)
       WRITE(LURES,803) PROPS(11),PROPS(12)
      ENDIF
C
      III14=INT(PROPS(14))
      IF(III14.EQ.0) THEN
       WRITE(LURES,840)              ! defaults: see mass32.f & frin32.f
      ELSE
       WRITE(LURES,841)
       WRITE(LURES,803) PROPS(15),PROPS(16)
      ENDIF
C
      IFRIC=INT(PROPS(51))
      IF(IFRIC.EQ.1) THEN
       IFRIM=INT(PROPS(52))
       IF(IFRIM.EQ.1) THEN
        IFRCO=INT(PROPS(53))
        IF(IFRCO.EQ.1) THEN
         WRITE(LURES,851)
         WRITE(LURES,803) PROPS(54),PROPS(55)
        ENDIF
       ENDIF
      ENDIF
C
      RETURN
C
  802 FORMAT(/)
  803 FORMAT(E15.6,10X,E15.6,10X,E15.6,10X,E15.6,/)
  810 FORMAT(/,'NORMAL STIFFNESS',/)
  811 FORMAT(/,'NORMAL STIFFNESS (MODEL=1)',/)
  812 FORMAT(/,'NORMAL STIFFNESS (MODEL=2)',/)
  813 FORMAT(/,'NORMAL STIFFNESS (MODEL=3)',/)
  814 FORMAT(/,'NORMAL STIFFNESS (MODEL=4)',/)
  815 FORMAT(/,'NORMAL STIFFNESS (MODEL=5)',/)
  816 FORMAT('NORMAL CONTACT COMPRESSIVE STIFFNESS=',E15.6,/,
     .       'MAXIMUM NORMAL CONTACT STRESS=',E15.6,/,
     .       'NORMAL CRITICAL FRACTURE ENERGY=',E15.6,/,
     .       'PERCENTAGE OF CONTACT GAP AT DEBOUNDING=',E15.6,/,
     .       'NORMAL CONTACT TENSILE STIFFNESS (COMPUTED)=',E15.6,/,
     .       'CONTACT GAP AT MAXIMUM STRESS (COMPUTED)=',E15.6,/,
     .       'CRITICAL DEBOUNDING PARAMETER=',E15.6,/,
     .       'TANGENTIAL CONTACT STIFFNESS IN COMPRESSION=',E15.6,/)
  820 FORMAT(/,'CONTACT TEMPERATURE',/)
  830 FORMAT(/,'JACOBIAN REGULARIZATION FACTORS (DEFAULTS)',/)
  831 FORMAT(/,'JACOBIAN REGULARIZATION FACTORS',/)
  840 FORMAT(/,'MINIMUM & MAXIMUM ADMISSIBLE CONTACT GAPS (DEFAULTS)',
     .       /,'-1.0D-08  1.0D+20 (INFINITE)',/)
  841 FORMAT(/,'MINIMUM & MAXIMUM ADMISSIBLE CONTACT GAPS',/)
  851 FORMAT(/,'FRICTION: COULOMB LAW NO. 1',/,
     .         'TANGENTIAL STIFFNESS  FRICTIONAL COEFFICIENT',/)
      END
