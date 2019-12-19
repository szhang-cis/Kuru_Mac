      SUBROUTINE INTLODT(ELDATT,ELPRET,ELVART,ELMATT,HTLODT,IFFIXT,
     .                   PRESCT,LNODST,MATNOT,PROELT,PROPST,RLOADT,
     .                   RLOAHT,ADVELT,TEMPIT,COORDT,FPCHAT,DISPLT,
     .                   LACTIT,WORK1T)
C***********************************************************************
C
C**** THIS SUBROUTINE SETS UP THE TYPE OF ALGORITHM AND REFERENCE LOAD
C     TO BE USED IN THE CURRENT TIME INTERVAL
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
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'
C
      COMMON/LDFILET/ITAPET
C
      COMMON/CONVE3/ICONV2
C
      DIMENSION MATNOT(NELEMT),          LNODST(NNODET,NELEMT),
     .          PROELT(NPRELT,NGRUPT),   PROPST(NPROPT,NMATST),
     .          ELDATT(NDATAT),          ELPRET(NPREVT),
     .          ELVART(NSTATT),          ELMATT(NMATXT),
     .          WORK1T(*),               RLOADT(*),
     .          HTLODT(NHLODT,NSUBFT,*), IFFIXT(NTOTVT,*)
      DIMENSION PRESCT(NTOTVT,2),        RLOAHT(NTOTVT,NFUNCT)
      DIMENSION ADVELT(NTOTVT*NDIMET),   TEMPIT(NPOINT,2),
     .          COORDT(NDIMET,NPOINT),   FPCHAT(NFPCH,NPOINT),
     .          DISPLT(NTOTVM),          LACTIT(NELEMT)
C
      PARAMETER (MCOMMT=8)
      CHARACTER*5 COMMDT(MCOMMT)
      DATA COMMDT/'HEAT' ,'ADVEC','BOUND','ACTIV','STRAT','FUNCT',
     .            'PROPE','END_I'/
C
      NEWLOT=-1
      NEWBOT=-1
      NEWSTT=-1
      NEWFUT=-1
      NEWADT=-1
      NEWACT=-1
      IF(IRESTT.EQ.1) THEN
       NEWLOT=0
       NEWBOT=0
      ENDIF
C
C**** LOOK FOR INTERVAL CARD KEY WORD WHICH DEFINES IF NEW 
C     HEAT OR/AND NEW INTERVAL DATA HAS TO BE READ.
C
  100 NPRINT=0
      ITAPET=LUDATT
      CALL LISTENT('INTLODT',NPRINT,ITAPET)
C
      DO ICOMMT=1,MCOMMT
       IF(WORDST(1).EQ.COMMDT(ICOMMT)) GO TO 20
      END DO
      GO TO 1000
C
   20 GO TO (1,2,3,4,5,6,7,8),ICOMMT
C
C**** IF REQUIRED SET UP NEW REFERENCE HEAT
C
    1 NEWLOT=0
      IF(WORDST(2).EQ.'NEW_H') NEWLOT=1
      IF(NACTIT.EQ.1) THEN
       IF(NEWACT.EQ.1) THEN
        IF(ISUVOACT.EQ.1) THEN
         IF(NEWLOT.EQ.0)
     .    CALL RUNMENT('WARNING: SURFACE OR VOLUME SHOULD BE REDEFINED')
        ENDIF
       ENDIF
      ENDIF
      IF(NEWLOT.EQ.1) THEN
       IF(PARAMT(1).NE.0.) ITAPET=INT(PARAMT(1))
       CALL LOADPST(ELDATT,ELPRET,ELVART,ELMATT,
     .              LNODST,MATNOT,PROELT,PROPST,RLOADT,RLOAHT,TEMPIT,
     .              WORK1T(ILOSUT(1)),ADVELT,COORDT,
     .              FPCHAT,DISPLT,LACTIT)
C
       NPRINT=0
       IF(ITAPET.NE.LUDATT) THEN
        ITAPET=LUDATT
        CALL LISTENT('INTLODT',NPRINT,ITAPET)
       ENDIF
       IF(WORDST(1).EQ.'END_H') GO TO 100
       WRITE(LUREST,909)
       CALL RUNENDT('ERROR: RUN STOPPED FOR HEAT INPUT')
      ENDIF                  ! newlot.eq.1
C
      GOTO 100
C
C**** IF REQUIRED SET UP ADVECTIVE VELOCITY
C
    2 NEWADT=0
      IF(ICONVT.EQ.0)
     . CALL RUNENDT('ERROR: ADVECT. CARD MUST BE INPUT IN PROBLEM_DATA')
      IF(WORDST(2).EQ.'NEW_A') NEWADT=1
      IF(NEWADT.EQ.1) THEN
       IF(PARAMT(1).NE.0.) ITAPET=INT(PARAMT(1))
       ICONV2=INT(PARAMT(2))                             ! to be revised
       IPRINX=0
       IF(WORDST(3).EQ.'NO_PR') IPRINX=1 ! to avoid a large results file
       CALL ADVECTC(ADVELT,FPCHAT,IPRINX)
C
       NPRINT=0
       ITAPET=LUDATT
       CALL LISTENT('INTLODT',NPRINT,ITAPET)
       IF(WORDST(1).EQ.'END_A') GO TO 100
       WRITE(LUREST,911)
       CALL RUNENDT('RUN STOPPED FOR INPUT ERROR        ')
      ENDIF                    ! newadt.eq.1
C
      GOTO 100
C
C**** IF REQUIRED SET UP NEW BOUNDARY CONDITION
C
    3 NEWBOT=0
      IF(WORDST(2).EQ.'NEW_B') NEWBOT=1
      IF(NEWBOT.EQ.1) THEN
       IF(PARAMT(1).NE.0.) ITAPET=INT(PARAMT(1))
       CALL FIXITYT(IFFIXT,LNODST,PRESCT,WORK1T(ILOSUT(15)),
     .              ELDATT,MATNOT,PROELT,WORK1T(ILOSUT( 1)) )
C
       NPRINT=0
       ITAPET=LUDATT
       CALL LISTENT('INTLODT',NPRINT,ITAPET)
       IF(WORDST(1).EQ.'END_B') GO TO 100
       WRITE(LUREST,910)
       CALL RUNENDT('RUN STOPPED FOR INPUT ERROR        ')
      ENDIF                    ! newbot.eq.1
C
      GOTO 100
C
C**** IF REQUIRED SET UP ACTIVE ELEMENTS
C
    4 NEWACT=0
      IF(NACTIT.EQ.0)
     . CALL RUNENDT('ERROR: ACTIVE CARD MUST BE INPUT IN PROBLEM_DATA')
      IF(WORDST(2).EQ.'NEW_A') NEWACT=1
      IF(NEWBOT.EQ.-1)                         ! control (see activet.f)
     . CALL RUNENDT('ERROR: BOUNDARY DATA BEFORE ACTIVE DATA')
      IF(NEWACT.EQ.1) THEN
       IF(PARAMT(1).NE.0.) ITAPET=INT(PARAMT(1))
       CALL ACTIVET(LACTIT,IFFIXT,LNODST)
C
       NPRINT=0
       ITAPET=LUDATT
       CALL LISTENT('INTLODT',NPRINT,ITAPET)
       IF(WORDST(1).EQ.'END_A') GO TO 100
       WRITE(LUREST,912)
       CALL RUNENDT('RUN STOPPED FOR INPUT ERROR        ')
      ENDIF                    ! newact.eq.1
C
      GOTO 100
C
C**** IF NECESSARY SET UP NEW SOLUTION STRATEGY
C
    5 NEWSTT=0
      IF(WORDST(2).EQ.'NEW_S') NEWSTT=1
      IF(NEWSTT.EQ.1) CALL STRATET
C
      GO TO 100
C
C**** IF REQUIRED SET UP NEW FUNCTION
C
    6 NEWFUT=0
      IF(WORDST(2).EQ.'NEW_F') NEWFUT=1
      IF(NEWFUT.EQ.1) CALL INTFUNT(HTLODT)
C
      GO TO 100
C
C**** IF REQUIRED SET UP NEW MAT. PROPERTIES (at INTERVAL_DATA level)
C
    7 NEWPRT=0
      IF(WORDST(2).EQ.'NEW_P') NEWPRT=1
      IF(NEWPRT.EQ.1) THEN
       ITAPET=LUDATT
       IF(PARAMT(1).NE.0.0) ITAPET=INT(PARAMT(1))
       CALL INPPROT(ITAPET,PROPST,    1)
      ENDIF
C
      GO TO 100
C
C**** CHECK THE CORRECTNESS OF THE INPUT
C
    8 CONTINUE
C
      IF(IRESTT.EQ.0) THEN
       IF(NEWFUT.EQ.-1.OR.NEWLOT.EQ.-1.OR.NEWBOT.EQ.-1.OR.
     .    NEWSTT.EQ.-1)THEN
        WRITE(LUREST,915)
        GO TO 1001
       ENDIF
C
       IF(ITIMET.EQ.1.AND.(NEWFUT.EQ.0.OR.NEWLOT.EQ.0.OR.
     .    NEWBOT.EQ.0.OR.NEWSTT.EQ.0)) THEN
        WRITE(LUREST,920)
        GO TO 1001
       ENDIF
C
       IF(ICONVT.EQ.1.AND.NEWADT.EQ.-1) THEN
        WRITE(LUREST,930)
        GO TO 1001
       ENDIF
C
       IF(NACTIT.EQ.1.AND.NEWACT.EQ.-1) THEN
        WRITE(LUREST,931)
        GO TO 1001
       ENDIF
      ENDIF            ! irestt.eq.0
C
      IF(IRESTT.EQ.1.AND.ISKIPT.EQ.1.AND.NEWSTT.EQ.-1)
     . WRITE(LUREST,925)
C
      RETURN
C
 1000 CALL RUNENDT('INTLOD:ERROR IN INTERVAL DATA BLOCK')
 1001 CALL RUNENDT('INTLOD:STOP FOR LACK OF INFORMATION')
C
  909 FORMAT(
     .'THE LIST OF HEATS DOES NOT END WITH AN        ',/,
     .'"END_HEAT_DATA" COMMAND CARD. VERIFY THE INPUT DATA.   .',/)
  910 FORMAT(
     .'THE LIST OF BOUNDARY CONDITIONS DOES NOT END WITH AN        ',/,
     .'"END_BOUNDARY_DATA" COMMAND CARD. VERIFY THE INPUT DATA.   .',/)
  911 FORMAT(
     .'THE LIST OF ADVECTIVE VELOCITY DOES NOT END WITH AN        ',/,
     .'"END_ADVECTIVE_VELOCITY_DATA" COMMAND CARD. VERIFY THE INPUT',/,
     .' DATA.   .',/)
  912 FORMAT(
     .'THE LIST OF ACTIVE ELEMENTS DOES NOT END WITH AN        ',/,
     .'"END_ACTIVE_ELEMENTS_DATA" COMMAND CARD. VERIFY THE INPUT',/,
     .' DATA.   .',/)
  915 FORMAT(
     .'EACH INTERVAL DATA BLOCK MUST HAVE HEAT, BOUNDARY &
     . STRATEGY ',/,
     .'COMMAND CARDS EVEN IF YOU WANT TO KEEP THE SAME DATA FOR ANY',/,
     .'OF THEM. IN THIS CASE JUST SPECIFY THE QUALIFIER AS "OLD".', /)
  920 FORMAT(
     .' THIS IS SUPPOSED TO BE A NEW RUN THEREFORE IT IS EXPECTED',   /,
     .' TO CONTAIN DATA FOR THE INTERVAL TIME PARAMETERS, THE TYPE',  /,
     .' OF HEAT AND/OR THE FIXITY CONDTIONS.',/)
  925 FORMAT(
     .'           * * * W A R N I N G * * *                        ',//,
     .'THIS IS SUPPOSED TO BE A RESTART TO COMPLETE THE INTERVAL   ',/,
     . I5,'FROM THE TIME ',E12.6,'. IN THIS CASE YOU MOST PROBABLY ',/,
     .'WANTED TO CHANGE THE SOLUTION STRATEGY BUT NO INFORMATION   ',/,
     .'ABOUT IT HAS BEEN GIVEN. CHECK YOUR DATA FILE IN THE RESTART',/,
     .'INTERVAL DATA .',/)
  930 FORMAT(
     .'ADVECTIVE VELOCITY MUST BE INPUT FOR ADVECTIVE PROBLEMS',/)
  931 FORMAT(
     .'ACTIVE ELEMENTS MUST BE INPUT FOR ACTIVE ELEMENTS PROBLEMS',/)
      END
