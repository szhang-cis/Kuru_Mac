      SUBROUTINE INPCON
C***********************************************************************
C
C**** THIS ROUTINE READS THE CONTACT CONTROL PARAMETERS
C
C     Note: these parameters are the same for every contact set
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'inpo_om.f'
C
      COMMON/TUNING4/RITEN,RITEF,TRUPL,TRUPM,TOLGA,TOLGAM
C
      ITAPE=LUCON                   ! always external file
      IERR1=1
      IF(ITAPE.NE.LUCON) GOTO 1000
      IF(ITAPE.EQ.LUCON) THEN       ! external .con file
       IF(IPRCO.EQ.1) OPEN(UNIT=LUCON,FILE=CI1,STATUS='OLD',ERR=1000)
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
C**** READ CONTACT CONTROL PARAMETERS
C
      NPRIN=1
      CALL LISTEN('INPCON',NPRIN,ITAPE)
C
      IAUXX=DINT(PARAM(1))                ! not used
      TRUPL=PARAM(2)
      TRUPM=PARAM(3)
C
C**** PRINT CONTACT CONTROL PARAMETERS
C
      WRITE(LURES,800)
      WRITE(LURES,810) IAUXX,TRUPL,TRUPM
C
      RETURN
  800 FORMAT(3X,'STEP',6X,'PARAM1',6X,'PARAM2',/)
  810 FORMAT(2X,I5,2(2X,E10.3),/)
  900 FORMAT(//,6X,'CONTACT CONTROL PARAMETERS',/,6X,26('-'),/)
  901 FORMAT(' ERROR IN OPENING INPUT CONTACT CONTROL FILE 48 ')
      END
