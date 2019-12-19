      SUBROUTINE INPSET(ITAPE,PROEL)
C***********************************************************************
C
C**** THIS ROUTINE READS THE MATERIAL PROPERTIES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'inpo_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION PROEL(NPREL,*)
      DIMENSION NDATO(2,5)
C
      DATA NDATO/3,3,   ! 2D PLANE STRESS
     .           3,4,   ! 2D PLANE STRAIN
     .           4,4,   ! 2D AXISYMMETRIC
     .           6,6,   ! 3D
     .           1,1/   ! 1D
C
      IERR1=1
      IF(ITAPE.NE.LUDAT) THEN
       ITAPE=LUSET       ! external .set file
       OPEN(UNIT=LUSET,FILE=CB1,STATUS='OLD',ERR=1000)
       IERR1=0
C
 1000  IF(IERR1.NE.0) THEN
        IF(IERR1.EQ.1) WRITE(LURES,901)
        CALL RUNEND('ERROR IN OPENING FILE         ')
       ENDIF
      ENDIF
C
      DO 10 IGRUP=1,NGRUP
C
C**** READ ELEMENTS CHARACTERISTICS
C
       NPRIN=1
       CALL LISTEN('INPSET',NPRIN,ITAPE)
       ISETS=INT(PARAM(1))
       IMATS=INT(PARAM(2))
       ITYPE=INT(PARAM(3))
       NTYPE=INT(PARAM(4))
       IRULE=INT(PARAM(5))
       IGAUS=INT(PARAM(6))
       THICK=    PARAM(7)
C
       IF(ISETS.LE.0.OR.ISETS.GT.NGRUP)
     .  CALL RUNEND('WRONG ISETS VALUE')
       IF(IMATS.LE.0.OR.IMATS.GT.NMATS)
     .  CALL RUNEND('WRONG IMATS VALUE')
       IF(ITYPE.NE.3.AND.ITYPE.NE.4.AND.ITYPE.NE.30.AND.ITYPE.NE.32.AND.
     .  ITYPE.NE.33)
     .  CALL RUNEND('WRONG ITYPE VALUE')
       IF(NTYPE.LE.0.OR.NTYPE.GT.5)
     .  CALL RUNEND('WRONG NTYPE VALUE')
       IF(IRULE.LE.0.OR.IRULE.GT.14)
     .  CALL RUNEND('WRONG IRULE VALUE')
       IF(IRULE.LE.10.AND.IRULE.GT.8)
     .  CALL RUNEND('WRONG IRULE VALUE')
       IF(IGAUS.LE.0.OR.IGAUS.GT.NGAUS)
     .  CALL RUNEND('WRONG IGAUS VALUE')
C
       PROEL(1,ISETS)=FLOAT(IMATS)
       PROEL(3,ISETS)=FLOAT(IRULE)
       PROEL(4,ISETS)=FLOAT(IGAUS)
       IF(ITYPE.EQ.0) ITYPE=1
       PROEL(5,ISETS)=FLOAT(ITYPE)
       PROEL(6,ISETS)=FLOAT(NTYPE)
C
C**** SOLID ELEMENTS
C
       IF(ITYPE.EQ.30) THEN
        IF(NTYPE.EQ.1.OR.NTYPE.EQ.2) THEN      ! 2D
         IF(THICK.EQ.0.0D0)
     .    CALL RUNEND('ERROR: THICK VALUE EQUAL ZERO')
         PROEL(7,ISETS)=THICK
        ENDIF
        IF(NTYPE.EQ.3.OR.NTYPE.EQ.4) THEN      ! axisymm. & 3D
         IF(THICK.NE.0.0D0)
     .    CALL RUNMEN('WARNING: THICK VALUE NOT USED')
         PROEL(7,ISETS)=1.0D0
        ENDIF
        IF(NTYPE.EQ.5) THEN                    ! 1D
         IF(THICK.EQ.0.0D0)
     .    CALL RUNEND('ERROR: CROSS AREA VALUE EQUAL ZERO')
         PROEL(7,ISETS)=THICK
        ENDIF
       ENDIF
C
C**** CONTACT & LINKING ELEMENTS
C
       IF(ITYPE.EQ.3.OR.ITYPE.EQ.4) THEN
        PROEL(7,ISETS)=1.0D0
       ENDIF
C
       IF(ITYPE.EQ.32) THEN
        IF(NTYPE.EQ.1.OR.NTYPE.EQ.2.OR.NTYPE.EQ.3) THEN     ! 2D & axys.
         IF(NOCOI.GT.0) THEN
          ISETC=INT(PARAM(7))
          IMAES=INT(PARAM(8))
          THICK=    PARAM(9)
          PROEL(22,ISETS)=1.0D0                ! non-coinc. contact set
          IF(NOCOI.EQ.2.AND.IMAES.EQ.0)
     .     PROEL(22,ISETS)=0.0D0               ! coinc. contact set
         ENDIF
         IF(THICK.NE.0.0D0)
     .    CALL RUNMEN('WARNING: THICK VALUE NOT USED')
         PROEL(7,ISETS)=1.0D0
        ENDIF
        IF(NTYPE.EQ.4) THEN                    ! impossible
         CALL RUNEND('ERROR: 3D FOR CONTACT ELEMENT IS NOT POSSIBLE')
        ENDIF
        IF(NTYPE.EQ.5) THEN                    ! 1D
         IF(NOCOI.GT.0) THEN
          ISETC=INT(PARAM(7))
          IMAES=INT(PARAM(8))
          THICK=    PARAM(9)
          PROEL(22,ISETS)=1.0D0                ! non-coinc. contact set
          IF(NOCOI.EQ.1.AND.IMAES.EQ.0) THEN   ! old version
           THICK=PARAM(7)
           ISETC=0
          ENDIF
          IF(NOCOI.EQ.2.AND.IMAES.EQ.0) THEN
           THICK=PARAM(7)
           ISETC=0
           PROEL(22,ISETS)=0.0D0               ! coinc. contact set
          ENDIF
         ENDIF
         IF(THICK.EQ.0.0D0)
     .    CALL RUNEND('ERROR: CROSS AREA VALUE EQUAL ZERO')
         PROEL(7,ISETS)=THICK
        ENDIF
        IF(NOCOI.GT.0) THEN
         PROEL(19,ISETS)=FLOAT(ISETC)
         PROEL(20,ISETS)=FLOAT(IMAES)
        ENDIF
       ENDIF
C
C**** DG ELEMENTS
C
       IF(ITYPE.EQ.33) THEN
        IF(NTYPE.EQ.1.OR.NTYPE.EQ.2.OR.NTYPE.EQ.3) THEN     ! 2D & axys.
         IF(THICK.NE.0.0D0)
     .    CALL RUNMEN('WARNING: THICK VALUE NOT USED')
         PROEL(7,ISETS)=1.0D0
        ENDIF
        IF(NTYPE.EQ.4) THEN                    ! impossible
         CALL RUNEND('ERROR: 3D FOR DG ELEMENT IS NOT POSSIBLE')
        ENDIF
        IF(NTYPE.EQ.5) THEN                    ! 1D
         IF(THICK.EQ.0.0D0)
     .    CALL RUNEND('ERROR: CROSS AREA VALUE EQUAL ZERO')
         PROEL(7,ISETS)=THICK
        ENDIF
       ENDIF
C
C**** PRINT ELEMENTS CHARACTERISTICS
C
       IF(ITYPE.EQ.3) THEN                     ! link
        WRITE(LURES,900)
        WRITE(LURES,913) ISETS,IMATS,ITYPE
       ENDIF
C
       IF(ITYPE.EQ.4) THEN                     ! contact
        WRITE(LURES,900)
        WRITE(LURES,913) ISETS,IMATS,ITYPE
       ENDIF
C
       IF(ITYPE.EQ.30) THEN                    ! solid
        WRITE(LURES,900)
        IF(NTYPE.EQ.5) THEN                    ! 1D
         WRITE(LURES,910) ISETS,IMATS,ITYPE,NTYPE,IRULE,IGAUS,THICK
        ENDIF
        IF(NTYPE.EQ.1.OR.NTYPE.EQ.2) THEN      ! 2D
         WRITE(LURES,910) ISETS,IMATS,ITYPE,NTYPE,IRULE,IGAUS,THICK
        ENDIF
        IF(NTYPE.EQ.3.OR.NTYPE.EQ.4) THEN      ! axisymm. & 3D
         WRITE(LURES,912) ISETS,IMATS,ITYPE,NTYPE,IRULE,IGAUS
        ENDIF
       ENDIF
C
       IF(ITYPE.EQ.32) THEN                    ! contact & linking
        IF(NTYPE.EQ.5) THEN                    ! 1D
         IF(NOCOI.EQ.0.OR.(NOCOI.EQ.2.AND.IMAES.EQ.0)) THEN
          WRITE(LURES,900)
          WRITE(LURES,910) ISETS,IMATS,ITYPE,NTYPE,IRULE,IGAUS,THICK
         ELSE
          WRITE(LURES,902)
          WRITE(LURES,920) ISETS,IMATS,ITYPE,NTYPE,IRULE,IGAUS,ISETC,
     .                     IMAES,THICK
         ENDIF
        ENDIF
        IF(NTYPE.EQ.1.OR.NTYPE.EQ.2) THEN      ! 2D
         IF(NOCOI.EQ.0.OR.(NOCOI.EQ.2.AND.IMAES.EQ.0)) THEN
          WRITE(LURES,900)
          WRITE(LURES,912) ISETS,IMATS,ITYPE,NTYPE,IRULE,IGAUS
         ELSE
          WRITE(LURES,902)
          WRITE(LURES,922) ISETS,IMATS,ITYPE,NTYPE,IRULE,IGAUS,ISETC,
     .                     IMAES
         ENDIF
        ENDIF
        IF(NTYPE.EQ.3) THEN                    ! axisymm.
         IF(NOCOI.EQ.0.OR.(NOCOI.EQ.2.AND.IMAES.EQ.0)) THEN
          WRITE(LURES,900)
          WRITE(LURES,912) ISETS,IMATS,ITYPE,NTYPE,IRULE,IGAUS
         ELSE
          WRITE(LURES,902)
          WRITE(LURES,922) ISETS,IMATS,ITYPE,NTYPE,IRULE,IGAUS,ISETC,
     .                     IMAES
         ENDIF
        ENDIF
       ENDIF
C
       IF(ITYPE.EQ.33) THEN                    ! Discontinuous Galerkin
        IF(NTYPE.EQ.5) THEN                    ! 1D
         WRITE(LURES,900)
         WRITE(LURES,910) ISETS,IMATS,ITYPE,NTYPE,IRULE,IGAUS,THICK
        ENDIF
        IF(NTYPE.EQ.1.OR.NTYPE.EQ.2) THEN      ! 2D
         WRITE(LURES,900)
         WRITE(LURES,912) ISETS,IMATS,ITYPE,NTYPE,IRULE,IGAUS
        ENDIF
        IF(NTYPE.EQ.3) THEN                    ! axisymm.
         WRITE(LURES,900)
         WRITE(LURES,912) ISETS,IMATS,ITYPE,NTYPE,IRULE,IGAUS
        ENDIF
       ENDIF
C
C**** ESTABLISH DEFAULT TUNING PARAMETERS
C
       IF(ITYPE.EQ.3.OR.ITYPE.EQ.33) THEN      ! link
        PROEL( 8,ISETS)=0.0D0
        PROEL( 9,ISETS)=0.0D0
        PROEL(10,ISETS)=0.0D0
        PROEL(11,ISETS)=0.0D0
        PROEL(12,ISETS)=0.0D0
        PROEL(13,ISETS)=0.0D0
       ENDIF
C
       IF(ITYPE.EQ.4.OR.ITYPE.EQ.32) THEN      ! contact
        PROEL( 8,ISETS)=0.0D0      ! iiten
        PROEL( 9,ISETS)=10.0D0     ! iitef
        PROEL(10,ISETS)=1.0E-06    ! trupl
        PROEL(11,ISETS)=0.0D0      ! trupm
        PROEL(12,ISETS)=-1.0E-08   ! tolga
        PROEL(13,ISETS)=1.0E+20    ! tolgam
       ENDIF
C
       IF(ITYPE.EQ.30) THEN                    ! solid
        PROEL( 8,ISETS)=0.0D0      ! kpart
        PROEL( 9,ISETS)=1.0D0      ! kparb
        PROEL(10,ISETS)=1.0D0      ! kpari
        PROEL(11,ISETS)=10.0E+07   ! estab
        PROEL(12,ISETS)=0.0D0
        PROEL(13,ISETS)=0.0D0
       ENDIF
C
C**** DETERMINES NKOST (COMPONENTS OF MATERIAL CONSTITUTIVE TENSOR)
C
       NSTRE=NDATO(1,NTYPE)
       NSTRS=NDATO(2,NTYPE)
       NKOSA=(NSTRS-1)*NSTRS/2+NSTRS    ! symm. part of const. tensor
       IF(KSYMM.EQ.0) NKOSA=NSTRS*NSTRS ! unsymm. part of const. tensor
       IF(NKOSA.GT.NKOST) NKOST=NKOSA
C
C**** INITIALISES FATIGUE PARAMETERS
C
c      PROEL(15,ISETS)=0.0D0
c      PROEL(16,ISETS)=0.0D0
c      PROEL(17,ISETS)=0.0D0
c      PROEL(18,ISETS)=0.0D0
C
   10 CONTINUE
C
C**** LOOK FOR THE 'END_SETS' CARD
C
      NPRIN=0
      JTAPE=LUDAT
      CALL LISTEN('INPSET',NPRIN,JTAPE)
      IF(WORDS(1).NE.'END_S')
     . CALL RUNEND('INPSET: END_SETS CARD NOT FOUND')
C
      RETURN
  900 FORMAT(//,6X,'ELEMENT CHARACTERISTICS',/,6X,23('-'),
     .        /,5X,'GROUP',3X,'MATERIAL',3X,'TYPE1',3X,'TYPE2',
     .          3X,'INT. RULE',3X,'NO. PTS.',3X,'AREA/THICKNESS',/)
  901 FORMAT(' ERROR IN OPENING INPUT SET FILE 41 ')
  902 FORMAT(//,6X,'ELEMENT CHARACTERISTICS',/,6X,23('-'),
     .        /,5X,'GROUP',3X,'MATERIAL',3X,'TYPE1',3X,'TYPE2',
     .          3X,'INT. RULE',3X,'NO. PTS.',3X,'CONTACT SET',
     .          3X,'CONTACT INDEX',3X,'AREA/THICKNESS',/)
c 905 FORMAT(6I5,F10.3)
  910 FORMAT(5X,I5,6X,I5,3X,I5,3X,I5,7X,I5,6X,I5,7X,F10.3)
  912 FORMAT(5X,I5,6X,I5,3X,I5,3X,I5,7X,I5,6X,I5,7X,'---')
  913 FORMAT(5X,I5,6X,I5,3X,I5)
  920 FORMAT(5X,I5,6X,I5,3X,I5,3X,I5,7X,I5,6X,I5,7X,I5,6X,I5,7X,F10.3)
  922 FORMAT(5X,I5,6X,I5,3X,I5,3X,I5,7X,I5,6X,I5,9X,I5,11X,I5,7X,'---')
      END
