      SUBROUTINE RSCHEKS
C***********************************************************************
C
C     FIND THE STARTING RECORD (NREC0) FOR RESTARTING CONVERGED STEP
C
C***********************************************************************
C
C      POS  1 : KINTXS   LAST STORED INTERVAL NUMBER
C      POS  2 : KSTEXS   LAST STORED STEP NUMBER
C      POS  3 : TIMMXS   LAST STORED TIME
C      POS  4 : NRECXS   RECORD NUMBER OF THE START OF RESULTS
C      POS  5 : KTSTES   TOTAL NUMBER OF STORED STEPS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      DATA TOLLES/1.0E-07/
C
      CALL CPUTIMS(TIME1S)
C
C**** READ FIRST RECORD TO FIND LAST CONVERGED RESULT
C
      READ(LURSTS,REC=1) KINTXS,KSTEXS,TIMMXS,NRECXS,KTSTES
C
C**** FIND APPROPRIATE NREC0
C
      IF(ISKIPS.EQ.0) THEN
C
C**** IF THIS IS A CONTINUATION THEN THE STARTING TIME IS THE LAST
C     STORED
C 
         NRECCS=NRECXS
         READ(LURSTS,REC=NRECCS) KINTXS,KSTEXS,TIMMXS,NRECXS,KTSTXS
         GO TO 101
C
      ELSE
C
C**** FIND APPROPIATE RECORD
C
         NRECCS=NRECGS+IDATPS(1,2)*NELEMS   !   START OF CONVERGED DATA
C
         KTSTOS=KTSTES
         DO 100 INDEXS=1,KTSTOS
           READ(LURSTS,REC=NRECCS) KINTXS,KSTEXS,TIMMXS,NRECXS,KTSTXS
           IF(KSTEPS.NE.0) THEN
             IF((KINTXS.EQ.KTIMES).AND.(KSTEXS.EQ.KSTEPS)) GOTO 101
           ELSE
             IF(ABS(TIMSTS-TIMMXS).LT.TOLERS) GOTO 101
           ENDIF
           NRECCS=NRECCS+NLENCS
  100    CONTINUE
C
      ENDIF
C
C**** ERROR CONDITION
C
       CLOSE(LURSTS)
       WRITE(LURESS,905) TIMSTS
       CALL RUNENDS('  ERROR IN RESTART CARD            ')
C
C**** SUCCESFUL SEARCH
C
  101 CONTINUE   
      ITIMES=KINTXS
      TTIMES=TIMMXS
      KTSTES=KTSTXS
      WRITE(LURESS,900) TTIMES,ITIMES,KSTEXS
C
      CALL CPUTIMS(TIME2S)
      CPURSS=CPURSS+(TIME2S-TIME1S)
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