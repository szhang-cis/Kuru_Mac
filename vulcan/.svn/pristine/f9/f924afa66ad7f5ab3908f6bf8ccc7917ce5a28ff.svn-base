      SUBROUTINE INPFAT(ITAPE,PROEL)
C***********************************************************************
C
C**** THIS ROUTINE READS THE FATIGUE PARAMETERS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'inpo_om.f'
C
      DIMENSION PROEL(NPREL,*)
C
c     IERR1=1
c     IF(ITAPE.NE.LUDAT.AND.ITAPE.NE.LUSET) GOTO 1000
c     IF(ITAPE.EQ.LUSET) THEN       ! external .set file
c      OPEN(UNIT=LUSET,FILE=CB1,STATUS='OLD',ERR=1000)
c      IERR1=0
C
c1000  IF(IERR1.NE.0) THEN
c       IF(IERR1.EQ.1) WRITE(LURES,901)
c       CALL RUNEND('ERROR IN OPENING FILE         ')
c      ENDIF
c     ENDIF
C
      WRITE(LURES,900)
C
C**** INITIALIZATION
C
      DO IPREL=15,18
       DO IGRUP=1,NGRUP
        PROEL(IPREL,IGRUP)=0.0D0
       ENDDO
      ENDDO
C
C**** READ SET CHARACTERISTICS FOR FATIGUE
C
      DO IGRUP=1,NGRUP
C
       NPRIN=1
       CALL LISTEN('INPFAT',NPRIN,ITAPE)
C
C**** LOOK FOR 'SET' CARD
C
       IF(WORDS(1).NE.'SET') GO TO 2000
C
       IF(PARAM(1).EQ.0.0)
     .  CALL RUNEND('INPFAT: WRONG SET NUMBER FOR FATIGUE')
C
       ISETS=INT(PARAM(1))
C
C**** CONTROLS REPEATED SETS WITH FATIGUE
C
       IF(PROEL(15,ISETS).NE.0.0)
     .  CALL RUNEND('ERROR IN SET NUMBER')
C
       ITYPE=INT(PROEL(5,ISETS))
       IF(ITYPE.NE.30)
     .  CALL RUNEND('INPFAT: ERROR ITYPE FOR FATIGUE') 
C
       PROEL(15,ISETS)=1.0                     ! kfati
       CALL LISTEN('INPFAT',NPRIN,ITAPE)
       IF(WORDS(1).EQ.'N_CYC') THEN
        NCYCL=INT(PARAM(1))
        PROEL(16,ISETS)=PARAM(1)
       ELSE
        CALL RUNEND('INPFAT: NUMBER OF CYCLES IS LACKING')
       ENDIF
C
       CALL LISTEN('INPFAT',NPRIN,ITAPE)
       IF(WORDS(1).EQ.'REVER') THEN
        IAUX=0
        IF(WORDS(2).EQ.'CONST') THEN
         IREVF=0
         PROEL(17,ISETS)=0.0
         REVCO=PARAM(1)
         PROEL(18,ISETS)=PARAM(1)
         IAUX=1
        ENDIF
        IF(WORDS(2).EQ.'VARIA') THEN
         IREVF=1
         PROEL(17,ISETS)=1.0
         IAUX=1
        ENDIF
        IF(IAUX.EQ.0)
     .   CALL RUNEND('INPFAT: ERROR IN REVERSIB. DATA')
       ELSE
        CALL RUNEND('INPFAT: NUMBER OF CYCLES IS LACKING')
       ENDIF
C
C**** SOME CONTROLS
C
C     Notes:
C
C     NCYCL=0 is only useful to compute a variable reversibility factor.
C     Warning: when a variable reversibility factor R is used, it is
C     necessary to input at least a first interval with NCYCL=0 (to 
C     compute R) and a second one with NCYCL>0 (to compute the fatigue
C     reduction function).
C
      IF(NCYCL.EQ.0.AND.IREVF.EQ.0)
     . CALL RUNEND('INPFAT: NCYCL=0 & IREVF=0 FACT. ARE NOT COMPATIBLE')
C
C**** WRITES FATIGUE PARAMETERS
C
       WRITE(LURES,911) ISETS
       WRITE(LURES,912) NCYCL
       IF(IREVF.EQ.0) THEN
        WRITE(LURES,913) REVCO
       ELSE
        WRITE(LURES,914)
       ENDIF
C
      ENDDO                                ! igrup=1,ngrup
C
C**** LOOK FOR 'END_FATIGUE' CARD
C
      NPRIN=0
      JTAPE=LUDAT
      CALL LISTEN('INPFAT',NPRIN,JTAPE)
 2000 IF(WORDS(1).NE.'END_F')
     . CALL RUNEND('INPFAT: END_FATIGUE CARD NOT FOUND')
C
      RETURN
C
  900 FORMAT(//,6X,18HFATIGUE PARAMETERS,/,6X,18('-'))
  901 FORMAT(' ERROR IN OPENING PROPERTIES INPUT FILE 42 ')
  911 FORMAT(/,' SET NUMBER FOR FATIGUE ANALYSIS=',I5,/)
  912 FORMAT(/,' NUMBER OF CYCLES               =',I15,/)
  913 FORMAT(/,' REVERSIBILITY FACTOR= CONSTANT= ',E15.6,/)
  914 FORMAT(/,' REVERSIBILITY FACTOR= VARIABLE',/)
      END
