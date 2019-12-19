      SUBROUTINE PROPIEN(ITAPE,PROPS)
C***********************************************************************
C
C**** THIS ROUTINE READS THE MATERIAL PROPERTIES FOR THE ISOTROPIC, 
C     ELASTIC & NON TEMPERATURE-DEPENDENT SOLID MODEL
C     (ELEMENT NO. 30)
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
      PROPS(2)=10.0            ! elastic
      PROPS(3)=1.0             ! no temperature-dependent
      PROPS(59)=3.0            ! solidification
C
      PROPS(37)=0.0            ! non prescribed strains (see incstn.f)
      PROPS(19)=0.0
C
      PROPS(21)=0.0            ! no external damping
      PROPS(22)=0.0
C
      KPLA1=0
      KPLA2=0
      KPLA3=0
      KPLA4=0
      KPLA5=0
      KPLA6=0
      KPLA7=0
C
C**** READ & WRITE MATERIAL PROPERTIES
C
      NPRIN=0
C
C**** LOOK FOR 'DENSITY' CARD
C
      CALL LISTEN('PROPIEN',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'DENSI') THEN
       PROPS(12)=PARAM(1)
       WRITE(LURES,801) PROPS(12)
       WRITE(LURES,899)
      ELSE
       IF(KDYNA.EQ.1) THEN
        CALL RUNEND('PROPIEN: DENSITY CARD NOT FOUND')
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
      CALL LISTEN('PROPIEN',NPRIN,ITAPE)
    1 IF(WORDS(1).EQ.'YOUNG') THEN
       NPOI1=61
       NPOI2=NPOI1
       PROPS(NPOI1)=PARAM(1)
       WRITE(LURES,803) PROPS(NPOI1)
       WRITE(LURES,899)
      ELSE
       CALL RUNEND('PROPIEN: YOUNG_MODULUS CARD NOT FOUND')
      ENDIF
C
C**** LOOK FOR 'POISSON_RATIO' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIEN',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'POISS') THEN
       NPOI1=NPOI2+1
       NPOI2=NPOI1
       PROPS(NPOI1)=PARAM(1)
       WRITE(LURES,805) PROPS(NPOI1)
       WRITE(LURES,899)
      ELSE
       CALL RUNEND('PROPIEN: POISSON_RATIO CARD NOT FOUND')
      ENDIF
C
      NPOI1=NPOI2+1
      NPOI2=NPOI1
      PROPS(NPOI1)=1.0D+20                     ! infinity value for C^th
C
      PROPS(36)=32.0D0                         ! default: Von Mises
      PROPS(52)=32.0D0                         ! default: Von Mises
C
C**** LOOK FOR FREE ENERGY MODEL
C
      CALL LISTEN('PROPIEN',NPRIN,ITAPE)
      IF(WORDS(1).EQ.'FREE_') THEN
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
       GO TO 9
      ENDIF
C
C**** LOOK FOR 'END_MATERIAL' CARD
C
      NPRIN=0
      CALL LISTEN('PROPIEN',NPRIN,ITAPE)
    9 IF(WORDS(1).NE.'END_M')
     . CALL RUNEND('PROPIEN: END_MATERIAL CARD NOT FOUND   ')
      WRITE(LURES,899)
C
C**** CONTROLS DIMENSION OF PROPS
C
      IF(NPOI2.GT.NPROP) THEN
       WRITE(LURES,900) NPROP,NPOI2
       CALL RUNEND('PROPIEN: TOO MANY PROPERTIES TO READ')
      ENDIF
C
      RETURN
  801 FORMAT(/,'  CONSTANT DENSITY VALUE=',E15.6,/)
  802 FORMAT(E15.6,10X,E15.6)
  803 FORMAT(/,'  YOUNG MODULUS         =',E15.6,/)
  804 FORMAT(/,'  REFERENCE TEMPERATURE =',E15.6,/)
  805 FORMAT(/,'  POISSON RATIO         =',E15.6,/)
  817 FORMAT(/,3X,'FREE ENERGY MODEL=',I5,/)
  899 FORMAT(/)
  900 FORMAT(//,'TOO MANY MATERIAL PROPERTIES TO READ:',/,
     .      20X,'NPROP =',I5,5X,'NUMBE =',I5,/)
      END
