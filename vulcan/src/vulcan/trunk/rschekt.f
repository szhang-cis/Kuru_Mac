      SUBROUTINE RSCHEKT
C***********************************************************************
C
C     FIND THE STARTING RECORD (NREC0) FOR RESTARTING CONVERGED STEP
C
C***********************************************************************
C
C      POS  1 : KINTX    LAST STORED INTERVAL NUMBER
C      POS  2 : KSTEX    LAST STORED STEP NUMBER
C      POS  3 : TIMMX    LAST STORED TIME
C      POS  4 : NRECX    RECORD NUMBER OF THE START OF RESULTS
C      POS  5 : KTSTE    TOTAL NUMBER OF STORED STEPS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DATA TOLLET/1.0E-07/
C
      CALL CPUTIMT(TIME1T)
C
C**** READ FIRST RECORD TO FIND LAST CONVERGED RESULT
C
      READ(LURSTT,REC=1) KINTXT,KSTEXT,TIMMXT,NRECXT,KTSTET
C
C**** FIND APPROPRIATE NREC0
C
      IF(ISKIPT.EQ.0)  THEN
C
C**** IF THIS IS A CONTINUATION THEN THE STARTING TIME IS THE LAST
C     STORED
C 
         NRECCT=NRECXT
         READ(LURSTT,REC=NRECCT) KINTXT,KSTEXT,TIMMXT,NRECXT,KTSTXT
         GO TO 101
C
      ELSE
C
C**** FIND APPROPIATE RECORD
C
         NRECCT=NRECGT+IDATPT(1,2)*NELEMT   !   START OF CONVERGED DATA
C
         KTSTOT=KTSTET
         DO 100 INDEXT=1,KTSTOT
           READ(LURSTT,REC=NRECCT) KINTXT,KSTEXT,TIMMXT,NRECXT,KTSTXT
           IF(KSTEPT.NE.0) THEN
             IF((KINTXT.EQ.KTIMET).AND.(KSTEXT.EQ.KSTEPT)) GOTO 101
           ELSE
             IF(ABS(TIMSTT-TIMMXT).LT.TOLERT) GOTO 101
           ENDIF
           NRECCT=NRECCT+NLENCT
  100    CONTINUE
C
      ENDIF
C
C**** ERROR CONDITION
C
       CLOSE(LURSTT)
       WRITE(LUREST,905) TIMSTT
       CALL RUNENDT('  ERROR IN RESTART CARD            ')
C
C**** SUCCESFUL SEARCH
C
  101 CONTINUE   
      ITIMET=KINTXT
      TTIMET=TIMMXT
      KTSTET=KTSTXT
      WRITE(LUREST,900) TTIMET,ITIMET,KSTEXT
C
      CALL CPUTIMT(TIME2T)
      CPURST=CPURST+(TIME2T-TIME1T)
C
      RETURN
C
  900 FORMAT('1',///,80('='),//,
     .'         * * * RESTARTS OF COMPUTATION * * *',//,
     .'             STARTING TIME         :',E12.6,/,
     .'             BELONGING TO INTERVAL :',7X,I5,/,
     .'             AT STEP               :',7X,I5,//,80('='))
  905 FORMAT(
     .'   THE DESIRED RESTARTING TIME ',E12.6,'HAS NOT BEEN FOUND',
     .'   INSIDE THE TIME HYSTORY OF THE PREVIOUS RUN')
      END
