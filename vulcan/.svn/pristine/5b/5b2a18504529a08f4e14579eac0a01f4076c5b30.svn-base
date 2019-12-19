      SUBROUTINE COUIND
C***********************************************************************
C
C**** THIS ROUTINE REDEFINES THE DEFAULT CONTROL COUPLING PARAMETERS
C     FOR THIS RUN
C
C     FOR CHANGING THE DEFAULT PARAMETERS
C
C     PARAMETERS         CONTROLLING WORDS
C
C     ICOSTA             'C_COS'
C     ITERME             'C_TYP'
C     ISTAGG             'STAGG'
C     NPOINC             'C_POI'
C     NELEMC             'C_ELE'
C     COUFAC             'C_FAC'
C     CENKEL             'C_FAC'
C     NITERC             'C_FAC'
C
C
C
C     Note: this routine should be finished !!!!!!!
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
      INCLUDE 'addi_oms.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
      INCLUDE 'nued_om.f'
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'inpo_oms.f'
      INCLUDE 'prob_oms.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'auxl_omt.f'
C
c     PARAMETER (MCOMM=7)
c     CHARACTER*5 COMMD(MCOMM)
c     CHARACTER*6 NOTESC
C
c     CHARACTER*8     TITLEC
c     COMMON/TITLESCA/TITLEC(8)
C
c     DATA COMMD/'C_COS','C_TYP','STAGG','C_POI','C_ELE','C_FAC',
c    .           'END_C'/
C
c     IBAND1=0
c     IBAND2=0
C
c     IERR1C=1
c     OPEN(UNIT=LUDATC,FILE=CCOA,STATUS='OLD',ERR=1001)
c     IERR1C=2
c     OPEN(UNIT=LURESC,FILE=CCOB,STATUS='UNKNOWN',ERR=1001)
c     IERR1C=0
C
c1001 IF(IERR1C.NE.0) THEN
c      IF(IERR1C.EQ.1) WRITE(LURESC,901)
c      IF(IERR1C.EQ.2) WRITE(LURESC,902)
c     ENDIF
C
C**** READ HEADING CARD
C
c  10 READ(LUDATC,900) NOTESC,(TITLEC(I),I=1,8)
c     IF(NOTESC.EQ.'VULCAN') GO TO 99
c     GO TO 10
C
c  99 WRITE(LURESC,905) NOTESC,(TITLEC(I),I=1,8)
C
C**** READ CONTROL PARAMETERS
C
c     IPRINC=0
c     ITAPEC=LUDATC
c     CALL LISTENC('COUINP',IPRINC,ITAPEC)
c     IF(WORDS(1).NE.'COUPL') GO TO 2000
c1000 CALL LISTENC('COUINP',IPRINC,ITAPEC)
C
C**** IDENTIFY COMMAND
C
c     DO ICOMM=1,MCOMM
c      IF(WORDS(1).EQ.COMMD(ICOMM)) GOTO 100
c     ENDDO
c     GO TO 2000
C
C**** EXECUTE APROPRIATE COMMAND
C
c 100 CONTINUE
c     GOTO (1,2,3,4,5,6,7), ICOMM
C
c   1 CONTINUE
c     ICOSTA=INT(PARAM(1))
c     GO TO 1000
C
c   2 CONTINUE
c     ITERME=INT(PARAM(1))
C
c     ITERMG=0
c     ITERMP=0
c     ITERMD=0
c     IF(ITERME.EQ.2) THEN                      ! bidirectional coupling
c      ITERMG=1                                 ! defaults
c      ITERMP=0              ! should be 1
c      ITERMD=0              ! should be 1
c      DO I=2,4
c       IF(WORDS(I).EQ.'GAP_D') THEN         ! gap dependency
c        ITERMG=INT(PARAM(I))
c       ENDIF
c       IF(WORDS(I).EQ.'MECHA') THEN         ! mechanical coupling terms
c        ITERMP=INT(PARAM(I))
c       ENDIF
c       IF(WORDS(I).EQ.'DEFOR') THEN         ! deformed shape
c        ITERMD=INT(PARAM(I))
c        IF(LARGE.EQ.0) THEN                    ! ??
c         ITERMD=0
c        ELSE
c         LARGET=LARGE
c        ENDIF
c       ENDIF
c      ENDDO
c     ENDIF
c     GO TO 1000
C
c   3 CONTINUE
c     ISTAGG=INT(PARAM(1))
c     GO TO 1000
C
c   4 CONTINUE
c     NPOINC=INT(PARAM(1))
c     IF(NPOINC.EQ.0) CALL RUNEND('COUINP: NPOINC=0')
c     IBAND1=1
c     GO TO 1000
C
c   5 CONTINUE
c     NELEMC=INT(PARAM(1))
c     IF(NELEMC.EQ.0) CALL RUNEND('COUINP: NELEMC=0')
c     IBAND2=1
c     GO TO 1000
C
c   6 CONTINUE
c     COUFAC=PARAM(1)
c     CENKEL=PARAM(2)
c     NITERC=INT(PARAM(3))
      NITERCS=0
c     NFCOUC=INT(PARAM(4))
c     GO TO 1000
C
c   7 CONTINUE
C
c     IF(IBAND1.NE.1)
c    . CALL RUNEND('NPOINC IN COUPLED FILE HAS NOT BEEN FOUND')
c     IF(IBAND2.NE.1)
c    . CALL RUNEND('NELEMC IN COUPLED FILE HAS NOT BEEN FOUND')
C
C**** LOOK FOR 'STOP' CARD
C
c     CALL LISTENC('COUINP',IPRINC,ITAPEC)
c     IF(WORDS(1).NE.'STOP') GO TO 2001
C
      RETURN 
c2000 CALL RUNEND('COUINP: ERROR IN COUPLED DATA BLOCK')
c2001 CALL RUNEND('COUINP: NO STOP CARD IN COUPLED DATA BLOCK')
C
c 900 FORMAT(A6,1X,8A8)
c 901 FORMAT(' ERROR IN OPENING OUTPUT  FILE  36 (35 LINUX)')
c 902 FORMAT(' ERROR IN OPENING PROCESS FILE 236 (36 LINUX)')
c 905 FORMAT(///5X,A6,1X,8A8//)
      END
