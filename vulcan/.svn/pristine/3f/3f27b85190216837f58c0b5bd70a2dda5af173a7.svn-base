      SUBROUTINE RSOPENS
C***********************************************************************
C
C**** THIS ROUTINE OPENS THE FOLLOWING FILES:
C
C     - FILE 201 (21 LINUX)  PROCESS DATA BASE FILE    ! not implemented
C     - FILE 211 (31 LINUX)  RESTART DATA BASE FILE    ! not implemented
C     - FILE 27 (16 LINUX)   OUTPUT FILE
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL OPENFILE
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_oms.f'
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'auxl_oms.f'
C
      CALL CPUTIMS(TIME1S)
C
      GO TO (101,102,102,101,101,102,101,101), NMACHIS
C
  101 CONTINUE
      LENRCS=512
      GO TO 1001
C
  102 CONTINUE
      LENRCS=512/4
      GO TO 1001
C
 1001 CONTINUE
C
      INQUIRE(UNIT=LUINFS, OPENED=OPENFILE)
      IF(.NOT.OPENFILE) THEN
        OPEN(UNIT=LUINFS,FILE=CTS,STATUS='UNKNOWN',IOSTAT=IERR)
      END IF
      IF(IEVFI.EQ.0) THEN
C
       IF(IRESTS.NE.0) GOTO 200
C***********************************************************************
C
C**** N E W   R U N  *** OPEN NEW RESTART AND OUTPUT FILES . . . . .
C
C***********************************************************************
       IERORT=1
       IF(NMACHIS.NE.8) THEN
        OPEN(UNIT=LURESS,FILE=CGS,STATUS='OLD',ERR=2000)
       ENDIF
       IERORT=0
       GOTO 2000
C***********************************************************************
C
C**** R E S T A R T  R U N *** OPEN OLD RESTART AND OUTPUT FILES . . . .
C
C***********************************************************************
 200   CONTINUE
       IERORT=1
       OPEN(UNIT=LURESS,FILE=CGS,STATUS='UNKNOWN',ACCESS='APPEND',
     .      ERR=2000)
       IERORT=0
C
 2000  IF(IERORT.NE.0)THEN
        IF(IERORT.EQ.2) WRITE(LURESS,902)
        IF(IERORT.EQ.3) WRITE(LURESS,903)
        CALL RUNENDS('  ERROR IN OPENING FILES           ')
       ENDIF
C
      ELSE           ! ievfi=1
C
       OPEN(UNIT=LUPRIS,FILE=CFS)
C
       IF(IRESTS.NE.0) GOTO 100
C***********************************************************************
C
C**** N E W   R U N  *** OPEN NEW RESTART AND OUTPUT FILES . . . . . 
C
C***********************************************************************
       IERORT=1
       IF(NMACHIS.NE.8) THEN
        OPEN(UNIT=LURESS,FILE=CGS,STATUS='OLD',ERR=1000)
       ENDIF
c      IERORT=2
c      IF(NDISKDS.EQ.0)
c    .  OPEN(UNIT=LUDTSS,FILE=CAS,STATUS='UNKNOWN',ACCESS='DIRECT',
c    .       FORM='UNFORMATTED',RECL=LENRCS,ERR=1000)
c      IERORT=3
c      IF(NFURESS.EQ.1) THEN
c       OPEN(UNIT=LUFANS,FILE=COS,STATUS='UNKNOWN',FORM='UNFORMATTED')
c       CLOSE(LUFANS)
c      ENDIF
c      IF(NFURESS.EQ.2)
c    .  OPEN(UNIT=LURSTS,FILE=CKS,STATUS='UNKNOWN',ACCESS='DIRECT',
c    .       FORM='UNFORMATTED',RECL=LENRCS,ERR=1000)
       IERORT=0
       GOTO 1000
C***********************************************************************
C
C**** R E S T A R T  R U N *** OPEN OLD RESTART AND OUTPUT FILES . . . .
C
C***********************************************************************
 100   CONTINUE
       IERORT=1
       OPEN(UNIT=LURESS,FILE=CGS,STATUS='UNKNOWN',ACCESS='APPEND',
     .      ERR=1000)
c      IERORT=2
c      IF(NDISKDS.EQ.0)
c    .  OPEN(UNIT=LUDTSS,FILE=CAS,STATUS='UNKNOWN',ACCESS='DIRECT',
c    .       FORM='UNFORMATTED',RECL=LENRCS,ERR=1000)
c      IERORT=3
c      OPEN(UNIT=LURSTS,FILE=CKS,STATUS='OLD',ACCESS='DIRECT',
c    .      FORM='UNFORMATTED',RECL=LENRCS,ERR=1000)
       IERORT=0
C
 1000  IF(IERORT.NE.0)THEN
        IF(IERORT.EQ.1) WRITE(LUPRIS,901)
        IF(IERORT.EQ.2) WRITE(LURESS,902)
        IF(IERORT.EQ.3) WRITE(LURESS,903)
        CALL RUNENDS('  ERROR IN OPENING FILES           ')
       ENDIF
C
      ENDIF          ! ievfi.eq.0
C
      CALL CPUTIMS(TIME2S)
      CPURSS=CPURSS+(TIME2S-TIME1S)
      RETURN
C
  901 FORMAT(' ERROR IN OPENING OUTPUT FILE 27 (16 LINUX)')
  902 FORMAT(' ERROR IN OPENING PROCESS FILE 201 (21 LINUX)')
  903 FORMAT(' ERROR IN OPENING RESTART FILE 211 (31 LINUX)')
      END
