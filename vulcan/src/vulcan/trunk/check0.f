      SUBROUTINE CHECK0
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
C     - the use of incrementation could not be compatible for a future
C       analysis
C
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
      INCLUDE 'inpo_om.f'
C
      PARAMETER (MCOMM=2)
      CHARACTER*5 COMMD(MCOMM)
      PARAMETER (MCOMN=5)
      CHARACTER*5 COMND(MCOMN)
      CHARACTER*6 NOTES
C
      DATA COMMD/'START','RESTA'/
      DATA COMND/'START','RESTA','DATAB','FUTUR','NON_S'/
C
C**** OPENING OF FILES
C
C     Note: the option using ILUS1 is more general but does not work
C           on Convex & Power Challenge; why?
C
      IERR1=1
      if(nmachim.ne.5) then
       OPEN(UNIT=LUDAT,FILE=CE,STATUS='OLD',ERR=1000)
      else
       OPEN(UNIT=LUDAT,FILE=CE(1:ILUS1),STATUS='OLD',ERR=1000)
      endif
      IERR1=2
      if(nmachim.ne.5) then
       OPEN(UNIT=LURES,FILE=CG,STATUS='UNKNOWN',ERR=1000)
      else
       OPEN(UNIT=LURES,FILE=CG(1:ILUS1),STATUS='UNKNOWN',ERR=1000)
      endif
      IERR1=0
C
 1000 IF(IERR1.NE.0) THEN
       IF(IERR1.EQ.1) WRITE(LURES,901)
       IF(IERR1.EQ.2) WRITE(LURES,902)
      ENDIF
  901 FORMAT(' ERROR IN OPENING OUTPUT  FILE 5 ')
  902 FORMAT(' ERROR IN OPENING PROCESS FILE 7 ')
C
C**** READ HEADING CARD
C
   10 READ (LUDAT,900) NOTES,(TITLE(I),I=1,8)
      IF(NOTES.EQ.'VULCAN') GO TO 100
      GO TO 10
C
  100 WRITE(LURES,905) NOTES,(TITLE(I),I=1,8)
C
C**** READ CARD AND IDENTIFY COMMAND
C
      NPRIN=0
      ITAPE=LUDAT
      CALL LISTEN('CHECK0',NPRIN,ITAPE)
      DO ICOMM=1,MCOMM
       IF(WORDS(1).EQ.COMMD(ICOMM)) GO TO 300
      END DO
      CALL RUNEND('CHECK0: ERROR READING LAST CARD    ')
C
C**** EXECUTE APPROPRIATE COMMAND
C
  300 GO TO(301,302) ICOMM
C***********************************************************************
C
C**** THIS IS A NEW RUN
C
C***********************************************************************
  301 CONTINUE
C
      IREST = 0
      ISKIP = 0
      INITI = 0
      INITV = 0
      ITIME = 0
      KTSTE = 0
      TIMST = 0.0
c     TIMST = PARAM(1)
      TTIME = TIMST
      IF(WORDS(2).EQ.'INITI') INITI=1
      IF(WORDS(2).EQ.'PREVI') INITI=2
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
      IREST = 2
      ISKIP = 0
      IF(WORDS(2).EQ.'SKIP')  ISKIP=1
      KSTEP = INT(PARAM(2))
      IF(KSTEP.EQ.0) THEN
       TIMST = PARAM(1)
       KTIME = 0
      ELSE
       TIMST = 0.0D+00
       KTIME = INT(PARAM(1))
      ENDIF
C
      GO TO 2000
C
 2000 CONTINUE
C
C**** OTHER PARAMETERS
C
C     Defaults values are assumed for (see setdat.f):
C
C     NMACHIM
C     NMEMOM,NMEMOM1-11
C
      IAXWO=3
      IAXPP=1
      DO IAXWP=IAXWO,MAXWP
       DO ICOMN=1,MCOMN
        IF(WORDS(IAXWP).EQ.COMND(ICOMN)) THEN
C
         GO TO(501,502,503,504,505) ICOMN
C
  501    TIMST=PARAM(IAXPP)                       ! starting time
         TTIME=TIMST
         IAXPP=IAXPP+1
         GO TO 500
C
  502    IF(INITI.NE.2)                           ! previous analysis
     .    CALL RUNEND('ERROR: INITI MUST BE 1 FOR NON-ST. INIT. COND.')
         TIMST=PARAM(IAXPP)                       ! restarting time
         TTIME=TIMST
         IAXPP=IAXPP+1
         GO TO 500
C
  503    IF(WORDS(IAXWP+1).EQ.'OUT_O') THEN
          NDISKDM=0                               ! database out of core
         ELSE
          IF(WORDS(IAXWP+1).EQ.'IN_CO') THEN
           NDISKDM=1                              ! database in core
          ELSE
           CALL RUNEND('ERROR: DATABASE SPECIFICATION IS LACKING')
          ENDIF
         ENDIF
         GO TO 500
C
  504    IF(WORDS(IAXWP+1).EQ.'FORM1') THEN       ! future restart
          NFURESM=1
         ELSE
          IF(WORDS(IAXWP+1).EQ.'FORM2') THEN
           NFURESM=2
          ELSE
           CALL RUNEND('ERROR: FUTURE ANALYSIS SPECIFICAT. IS LACKING')
          ENDIF
         ENDIF
         GO TO 500
C
  505    IF(INITI.NE.1)                ! non-standard initial conditions
     .    CALL RUNEND('ERROR: INITI MUST BE 1 FOR NON-ST. INIT. COND.')
         INITV=1
         IF(WORDS(IAXWP+1).EQ.'COMPO') THEN       ! number of components
          IF(PARAM(IAXPP).GT.0.0D0) THEN
           NPRE2=INT(PARAM(IAXPP))                ! stress
           IF(NPRE2.GT.6)
     .      CALL RUNEND('ERROR: NPRE2>6-Change dim. STRE1 in frin30.f')
           IAXPP=IAXPP+1
          ENDIF
          IF(PARAM(IAXPP).GT.0.0D0) THEN
           NPRE3=INT(PARAM(IAXPP))                ! internal variables
           IF(NPRE3.GT.8)
     .      CALL RUNEND('ERROR: NPRE3>8-Change dim. VARI1 in frin30.f')
           IAXPP=IAXPP+1
          ENDIF
          NPREA=NPRE1+NPRE2+NPRE3+NPRE4+NPRE5     ! dimension of PREAS
         ENDIF
         GO TO 500
C
  500    CONTINUE
        ENDIF
       END DO
      END DO
C
      IF(NFURESM.EQ.1) THEN
       OPEN(LUFAN,FILE=CO,STATUS='UNKNOWN',ERR=401)
       CLOSE(LUFAN,STATUS='DELETE')   !=>OPEN STATUS=UNKNOWN in rsopen.f
  401  CONTINUE
      ENDIF
      IF(NFURESM.EQ.2) THEN
       OPEN(LURST,FILE=CK,STATUS='OLD',ERR=402)
       CLOSE(LURST,STATUS='KEEP')
       GO TO 302
  402  CONTINUE
      ENDIF
C
      RETURN
 900  FORMAT(A6,1X,8A8)
 905  FORMAT(///5X,A6,1X,8A8//)
      END
