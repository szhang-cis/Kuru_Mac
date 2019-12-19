      SUBROUTINE CHECK0S
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
C     3) WORDS(1): EVOLUTION_EQUATIONS for the microstructural problem
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
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_oms.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nued_om.f'       ! thermal-microstructural
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
      INCLUDE 'inpo_oms.f'
C
      PARAMETER (MCOMMS=3)
      CHARACTER*5 COMMDS(MCOMMS)
      CHARACTER*6 NOTEST
C
      DATA COMMDS/'START','RESTA','EVOLU'/
C
C**** OPENING OF FILES
C
C     Note: the option using ILUS1T is more general but does not work
C           on Convex & Power Challenge; why?
C
      IERR1T=1
      if(nmachis.ne.5) then
       OPEN(UNIT=LUDATS,FILE=CES,STATUS='OLD',ERR=1000)
      else
       OPEN(UNIT=LUDATS,FILE=CES(1:ILUS1S),STATUS='OLD',ERR=1000)
      endif
      IERR1T=2
      if(nmachis.ne.5) then
       OPEN(UNIT=LURESS,FILE=CGS,STATUS='UNKNOWN',ERR=1000)
      else
       OPEN(UNIT=LURESS,FILE=CGS(1:ILUS1S),STATUS='UNKNOWN',ERR=1000)
      endif
      IERR1T=0
C
 1000 IF(IERR1T.NE.0) THEN
       IF(IERR1T.EQ.1) WRITE(LURESS,901)
       IF(IERR1T.EQ.2) WRITE(LURESS,902)
      ENDIF
  901 FORMAT(' ERROR IN OPENING OUTPUT  FILE 25 (15 LINUX)')
  902 FORMAT(' ERROR IN OPENING PROCESS FILE 27 (16 LINUX)')
C
C**** READ HEADING CARD
C
   10 READ(LUDATS,900) NOTEST,(TITLES(I),I=1,8)
      IF(NOTEST.EQ.'VULCAN') GO TO 100
      GO TO 10
C
  100 WRITE(LURESS,905) NOTEST,(TITLES(I),I=1,8)
C
C**** READ CARD AND IDENTIFY COMMAND
C
      NPRINT=0
      ITAPET=LUDATS
      CALL LISTENS('CHECK0S',NPRINT,ITAPET)
      DO ICOMMS=1,MCOMMS
       IF(WORDSS(1).EQ.COMMDS(ICOMMS)) GO TO 300
      END DO
      CALL RUNENDS('CHECK0S: ERROR READING LAST CARD   ')
C
C**** EXECUTE APPROPRIATE COMMAND
C
  300 GO TO(301,302,303) ICOMMS
C***********************************************************************
C
C**** THIS IS A NEW RUN
C
C***********************************************************************
  301 CONTINUE
C
      IEVFI=1
C
      IRESTS = 0
      ISKIPS = 0
      INITIS = 0
      ITIMES = 0
      KTSTES = 0
      TIMSTS = PARAMS(1)
      TTIMES = TIMSTS
      IF(WORDSS(2).EQ.'INITI') INITIS=1
      IF(WORDSS(2).EQ.'PREVI') INITIS=2
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
      IEVFI=1
C
      IRESTS = 2
      ISKIPS = 0
      IF(WORDSS(2).EQ.'SKIP')  ISKIPS=1
      KSTEPS = INT(PARAMS(2))
      IF(KSTEPS.EQ.0) THEN
       TIMSTS = PARAMS(1)
       KTIMES = 0
      ELSE
       TIMSTS = 0.0D+00
       KTIMES = INT(PARAMS(1))
      ENDIF
C
      GO TO 2000
C
C***********************************************************************
C
C**** EVOLUTION EQUATIONS
C
C***********************************************************************
  303 CONTINUE
C
      IEVFI=0
      IF(IMICR.EQ.0)
     .  CALL RUNENDS('ERROR: EVOLUTION EQUATIONS ONLY IN COUPLED PROB.')
      IF(WORDSS(2).EQ.'WEAK_') IMICO=0             ! weak coupling index
      RETURN
C
 2000 CONTINUE
C
C**** OTHER PARAMETERS
C
C     Defaults values are assumed for (see setdatt.f):
C
C     NMACHIS
C     NMEMOS,NMEMO1S-11S
C
      DO IAXWPS=3,MAXWPS
C
       IF(WORDSS(IAXWPS).EQ.'DATAB') THEN
        IF(WORDSS(IAXWPS+1).EQ.'OUT_O') THEN
         NDISKDS=0                                ! database out of core
        ELSE
         IF(WORDSS(IAXWPS+1).EQ.'IN_CO') THEN
          NDISKDS=1                               ! database in core
         ELSE
          CALL RUNENDS('ERROR: DATABASE SPECIFICATION IS LACKING')
         ENDIF
        ENDIF
       ENDIF
C
       IF(WORDSS(IAXWPS).EQ.'FUTUR') THEN
        IF(WORDSS(IAXWPS+1).EQ.'FORM1') THEN
         NFURESS=1
        ELSE
         IF(WORDSS(IAXWPS+1).EQ.'FORM2') THEN
          NFURESS=2
         ELSE
          CALL RUNENDS('ERROR: FUTURE ANALYSIS SPECIFICAT. IS LACKING')
         ENDIF
        ENDIF
       ENDIF
C
      ENDDO
C
      IF(NFURESS.EQ.1) THEN
c      OPEN(LUFANS,FILE=CNS,STATUS='UNKNOWN',ERR=401)
c      CLOSE(LUFANS,STATUS='DELETE') !=>OPEN STATUS=UNKNOWN in rsopens.f
  401  CONTINUE
      ENDIF
      IF(NFURESS.EQ.2) THEN
c      OPEN(LURSTS,FILE=CKS,STATUS='OLD',ERR=402)
c      CLOSE(LURSTS,STATUS='KEEP')
       GO TO 302
  402  CONTINUE
      ENDIF
C
      RETURN
C
 900  FORMAT(A6,1X,8A8)
 905  FORMAT(///5X,A6,1X,8A8//)
      END
