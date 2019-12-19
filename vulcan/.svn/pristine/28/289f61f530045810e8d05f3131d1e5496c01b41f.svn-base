      SUBROUTINE CHECK0T
C***********************************************************************
C
C**** THIS ROUTINE CHECKS IF THE PRESENT RUN IS A RESTART OR NOT
C     AND DEFINES THE CONTROLLING PARAMETERS IN COMMON 'RESTART' &
C     SOME IN 'CURREN' (THE REMAINING PARAMETERS OF COMMON/CURREN
C     ARE EITHER SET UP DURING EXECUTION OR ,FOR A RESTART, ARE
C     READ BACK IN SUB-RSTAR3)
C
C
C   1)  WORDS(1) = "START"
C
C       IREST = 0 : NEW RUN 
C
C       - WORDS(2) = "NOINItial" - INITI = 0 - No initial conditions
C
C       - WORDS(2) = "INITIal"   - INITI = 1 - Initial conditions to be
C                                              read afterwards
C
C       - WORDS(2) = "PREVIous"  - INITI = 2 - Previous data to be read
C                                              from restart file (*.pan)
C                                              (NFURES=1 in the previous
C                                              run creating *.fan)
C    
C
C   2)  WORDS(1) = "RESTArt" (old form)
C
C       IREST = 2 : RESTART (NFURES=2 in the previous run)
C
C       - WORDS(2) = "CONTInue"  - ISKIP = 0 -
C
C       READ ALL VARIABLES AND PARAMETERS FROM THE LAST TIME STEP WHICH
C       HAS BEEN COMPLETED AND CONTINUE THE ANALYSIS FROM THE LAST
C       STORED RESULTS (AT MOST RESET THE CPU TIME LIMIT AS INPUT BY
C       PARAM(2))
C
C       - WORDS(2) = "SKIP"      - ISKIP = 1 -
C
C       RESTART THE ANALYSIS FROM THE TIME STEP GIVEN IN PARAM(1)
C       PARAMETER AND IN CASE RESET THE CPU TIME LIMIT AS INPUT BY
C       PARAM(2)
C
C
C   FOR ANY IREST:
C
C       WORDS(3) = '     '         - NFURES=0 (no future restart)
C
C       WORDS(3) = FUTURe_analysis - NFURES=1 or 2
C
C       - WORDS(4) = 'FORM1' => NFURES =1
C
C       START (IREST=0) again reading the input data (modified in
C       comparison to the original) with the last stored PREVIous
C       (INITI=2) results as new "initial" data
C
C       - WORDS(4) = 'FORM2' => NFURES =2
C
C       RESTART (IREST=1)
C
C
C
C     NOTES: 
C
C     NFURES=2 is not available (January/1995)
C
C     For coupled problems:
C
C     - the previous analysis must contain the mechanical initial
C       conditions
C
C     - the use of subincrementation could not be compatible for a
C       future analysis
C
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'
C
      PARAMETER (MCOMMT=2)
      CHARACTER*5 COMMDT(MCOMMT)
      CHARACTER*6 NOTEST
C
      DATA COMMDT/'START','RESTA'/
C
C**** OPENING OF FILES
C
C     Note: the option using ILUS1T is more general but does not work
C           on Convex & Power Challenge; why?
C
      IERR1T=1
      if(nmachi.ne.5) then
       OPEN(UNIT=LUDATT,FILE=CET,STATUS='OLD',ERR=1000)
      else
       OPEN(UNIT=LUDATT,FILE=CET(1:ILUS1T),STATUS='OLD',ERR=1000)
      endif
      IERR1T=2
      if(nmachi.ne.5) then
       OPEN(UNIT=LUREST,FILE=CGT,STATUS='UNKNOWN',ERR=1000)
      else
       OPEN(UNIT=LUREST,FILE=CGT(1:ILUS1T),STATUS='UNKNOWN',ERR=1000)
      endif
      IERR1T=0
C
 1000 IF(IERR1T.NE.0) THEN
       IF(IERR1T.EQ.1) WRITE(LUREST,901)
       IF(IERR1T.EQ.2) WRITE(LUREST,902)
      ENDIF
  901 FORMAT(' ERROR IN OPENING OUTPUT  FILE 205 (25 LINUX)')
  902 FORMAT(' ERROR IN OPENING PROCESS FILE 207 (27 LINUX)')
C
C**** READ HEADING CARD
C
   10 READ(LUDATT,900) NOTEST,(TITLET(I),I=1,8)
      IF(NOTEST.EQ.'VULCAN') GO TO 100
      GO TO 10
C
  100 WRITE(LUREST,905) NOTEST,(TITLET(I),I=1,8)
C
C**** READ CARD AND IDENTIFY COMMAND
C
      NPRINT=0
      ITAPET=LUDATT
      CALL LISTENT('CHECK0T',NPRINT,ITAPET)
      DO ICOMMT=1,MCOMMT
       IF(WORDST(1).EQ.COMMDT(ICOMMT)) GO TO 300
      END DO
      CALL RUNENDT('CHECK0T: ERROR READING LAST CARD   ')
C
C**** EXECUTE APPROPRIATE COMMAND
C
  300 GO TO(301,302) ICOMMT
C***********************************************************************
C
C**** THIS IS A NEW RUN
C
C***********************************************************************
  301 CONTINUE
C
      IRESTT = 0
      ISKIPT = 0
      INITIT = 0
      ITIMET = 0
      KTSTET = 0
      TIMSTT = PARAMT(1)
      TTIMET = TIMSTT
      IF(WORDST(2).EQ.'INITI') INITIT=1
      IF(WORDST(2).EQ.'PREVI') INITIT=2
C
      GO TO 2000
C
C***********************************************************************
C
C**** THIS IS A RESTART RUN
C
C***********************************************************************
  302 CONTINUE
C
      IRESTT = 2
      ISKIPT = 0
      IF(WORDST(2).EQ.'SKIP')  ISKIPT=1
      KSTEPT = INT(PARAMT(2))
      IF(KSTEPT.EQ.0) THEN
       TIMSTT = PARAMT(1)
       KTIMET = 0
      ELSE
       TIMSTT = 0.0D+00
       KTIMET = INT(PARAMT(1))
      ENDIF
C
      GO TO 2000
C
 2000 CONTINUE
C
C**** OTHER PARAMETERS
C
C     Defaults values are assumed for (see setdatt.f):
C
C     NMACHI
C     NMEMO,NMEMO1-11
C
      DO IAXWPT=3,MAXWPT
C
       IF(WORDST(IAXWPT).EQ.'DATAB') THEN
        IF(WORDST(IAXWPT+1).EQ.'OUT_O') THEN
         NDISKD=0                                 ! database out of core
        ELSE
         IF(WORDST(IAXWPT+1).EQ.'IN_CO') THEN
          NDISKD=1                                ! database in core
         ELSE
          CALL RUNENDT('ERROR: DATABASE SPECIFICATION IS LACKING')
         ENDIF
        ENDIF
       ENDIF
C
       IF(WORDST(IAXWPT).EQ.'FUTUR') THEN
        IF(WORDST(IAXWPT+1).EQ.'FORM1') THEN
         NFURES=1
        ELSE
         IF(WORDST(IAXWPT+1).EQ.'FORM2') THEN
          NFURES=2
         ELSE
          CALL RUNENDT('ERROR: FUTURE ANALYSIS SPECIFICAT. IS LACKING')
         ENDIF
        ENDIF
       ENDIF
C
      ENDDO
C
      IF(NFURES.EQ.1) THEN
       OPEN(LUFANT,FILE=COT,STATUS='UNKNOWN',ERR=401)
       CLOSE(LUFANT,STATUS='DELETE') !=>OPEN STATUS=UNKNOWN in rsopent.f
  401  CONTINUE
      ENDIF
      IF(NFURES.EQ.2) THEN
       OPEN(LURSTT,FILE=CKT,STATUS='OLD',ERR=402)
       CLOSE(LURSTT,STATUS='KEEP')
       GO TO 302
  402  CONTINUE
      ENDIF
C
      RETURN
 900  FORMAT(A6,1X,8A8)
 905  FORMAT(///5X,A6,1X,8A8//)
      END
