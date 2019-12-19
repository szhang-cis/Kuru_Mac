      SUBROUTINE INPTUN(ITAPE,PROEL)
C***********************************************************************
C
C**** THIS ROUTINE READS THE "TUNING" PARAMETERS FOR THE SETS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'inpo_om.f'
C
      DIMENSION PROEL(NPREL,*)
C
      IERR1=1
      IF(ITAPE.NE.LUDAT) THEN
       ITAPE=LUTUN       ! external .tun file
       OPEN(UNIT=LUTUN,FILE=CH1,STATUS='OLD',ERR=1000)
       IERR1=0
C
 1000  IF(IERR1.NE.0) THEN
        IF(IERR1.EQ.1) WRITE(LURES,901)
        CALL RUNEND('ERROR IN OPENING FILE         ')
       ENDIF
      ENDIF
C
      WRITE(LURES,900)
C
C**** READ ELEMENTS CHARACTERISTICS
C
      DO 10 IGRUP=1,NGRUP
C
      NPRIN=1
      CALL LISTEN('INPTUN',NPRIN,ITAPE)
C
      ISETS=DINT(PARAM(1))
      IF(ISETS.EQ.0) GO TO 20
C
      ITYPE=INT(PROEL(5,ISETS))
C
C**** SOLID ELEMENTS
C
      IF(ITYPE.EQ.30) THEN
       PARAM1=PARAM(2)
       PARAM2=PARAM(3)
       PARAM3=PARAM(4)
       PARAM4=PARAM(5)
C
       PROEL( 8,ISETS)=PARAM1
       PROEL( 9,ISETS)=PARAM2
       PROEL(10,ISETS)=PARAM3
       PROEL(11,ISETS)=PARAM4
      ENDIF 
C
C**** CONTACT & LINKING ELEMENTS
C
      IF(ITYPE.EQ.3.OR.ITYPE.EQ.4.OR.ITYPE.EQ.32.OR.ITYPE.EQ.33) THEN
       PARAM1=PARAM(2)
       PARAM2=PARAM(3)
       PARAM3=PARAM(4)
       PARAM4=PARAM(5)
       PARAM5=PARAM(6)
C
       PROEL( 8,ISETS)=1.0          ! tunings considered
       PROEL( 9,ISETS)=PARAM1
       PROEL(10,ISETS)=PARAM2
       PROEL(11,ISETS)=PARAM3
       PROEL(12,ISETS)=PARAM4
       PROEL(13,ISETS)=PARAM5
      ENDIF
C
C**** PRINT ELEMENTS CHARACTERISTICS
C
      IF(ITYPE.EQ.3.OR.ITYPE.EQ.4.OR.ITYPE.EQ.32.OR.ITYPE.EQ.33) THEN
       WRITE(LURES,800)
       WRITE(LURES,810) ISETS,PARAM1,PARAM2,PARAM3
      ENDIF
      IF(ITYPE.EQ.30) THEN
       WRITE(LURES,801)
       WRITE(LURES,811) ISETS,PARAM1,PARAM2,PARAM3,PARAM4
      ENDIF
C
   10 CONTINUE
C
C**** LOOK FOR THE 'END_TUNING' CARD
C
   20 NPRIN=0
      IF(ITAPE.NE.LUDAT) THEN
       JTAPE=LUDAT
       CALL LISTEN('INPTUN',NPRIN,JTAPE)
      ENDIF
      IF(WORDS(1).NE.'END_T')
     . CALL RUNEND('INPTUN: END_TUNING CARD NOT FOUND')
C
      RETURN
  800 FORMAT(3X,'SET',6X,'PARAM1',6X,'PARAM2',6X,'PARAM3',/)
  801 FORMAT(3X,'SET',6X,'PARAM1',6X,'PARAM2',6X,'PARAM3',6X,'PARAM4',/)
  810 FORMAT(1X,I5,3(2X,E10.3),/)
  811 FORMAT(1X,I5,4(2X,E10.3),/)
  900 FORMAT(//,6X,'TUNING PARAMETERS',/,6X,23('-'),/)
  901 FORMAT(' ERROR IN OPENING INPUT TUNING FILE 47 ')
      END
