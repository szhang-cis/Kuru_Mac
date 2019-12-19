      SUBROUTINE STORESS(DISPRS,DISTOS,ELPRES,ELVARS,HEADSS,HTLODS,
     .                   IFFIXS,REFORS,RLOADS,TLOADS,TEMPIS,ILATES)
C***********************************************************************
C
C**** THIS ROUTINE STORES THE RESULTS OF THE CURRENT TIME STEP
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_oms.f'
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      DIMENSION HTLODS(NHLODS,NSUBFS,NFUNCS), IFFIXS(NTOTVS,2), 
     .          RLOADS(NTOTVS),               DISPRS(NTOTVS,3),
     .          DISTOS(NTOTVS,3),             HEADSS(NPOINS,4), 
     .          REFORS(NTOTVS),               TLOADS(NTOTVS,2),
     .          ELPRES(NPREVS),               ELVARS(NSTATS),
     .          TEMPIS(NTOTVS,2)
C
C=======================================================================
C
C**** RESTART FILE
C
C=======================================================================
C
      IF(NFURESS.EQ.1) THEN
       CALL RSTAR4S(TEMPIS,DISTOS,ELPRES,ELVARS)
C
c      WRITE(LUPRIS,*)
c      WRITE(LUPRIS,*) 'FUTURE RESTART INFORMATION'
c      WRITE(LUPRIS,*) 'RESULTS STORED FOR TIME:',TTIMES
c      WRITE(LUPRIS,*)
c      CLOSE(LUPRIS)
c      OPEN(UNIT=LUPRIS,FILE=CFS,STATUS='OLD',ACCESS='APPEND')
      ENDIF                      ! nfuress.eq.1
C
C**** DECIDE ON SAVING THESE RESULTS
C
      IF(NFURESS.EQ.2) THEN
C
       ISAVES=0
       IF(MOD(ISTEPS,NOUTPS(1)).EQ.0) ISAVES=1
       IF(ILATES.EQ. 1)               ISAVES=1
       IF(KSAVES.EQ.-1)               ISAVES=0
C
C**** DUMP TO DATA BASE & SET UP APPROPRIATE COUNTERS
C
       IF(ISAVES.EQ.1) THEN
        IF(KTSTES.EQ.0) KTSTES=1
        IF(KSAVES.EQ.1) KTSTES=KTSTES+1
        CALL RSTAR3S(DISPRS,DISTOS,ELPRES,ELVARS,HEADSS,HTLODS,IFFIXS,
     .               REFORS,RLOADS,TLOADS,    1)
        IF(KSAVES.EQ.1) NRECCS=NRECCS+NLENCS
       ENDIF
C
      ENDIF                      ! nfuress.eq.2
C
C=======================================================================
C
C**** PROCESS FILE
C
C=======================================================================
C
C**** INTERCHANGE CURRENT AND CONVERGED VALUES OF ELPRE & ELVAR
C
      DO INDEXS=1,5
       IDUMMS=IDATPS(2,INDEXS)
       IDATPS(2,INDEXS)=IDATPS(4,INDEXS)
       IDATPS(4,INDEXS)=IDUMMS
       IDUMMS=IDATPS(3,INDEXS)
       IDATPS(3,INDEXS)=IDATPS(5,INDEXS)
       IDATPS(5,INDEXS)=IDUMMS
      ENDDO
C
C**** WRITE DISTO TO PROCESS FILE
C
      CALL DATBASS(DISTOS,  12,    1)
C
C=======================================================================
C
C**** RESULTS FILE
C
C=======================================================================
C
      WRITE(LUPRIS,900)
  900 FORMAT(50X,'>>>>>>> STEP WRITTEN UP',/)
ctm
ctm   cambio en el sistema operativo del CONVEX (febrero de 1993)
ctm
c      CLOSE(LURESS)
c      OPEN(UNIT=LURESS,STATUS='OLD',ACCESS='APPEND')
ctm
C
C**** INITIALISE FLAGS IREST,ISKIP
C
      IRESTS=0
      ISKIPS=0
C
      RETURN
      END
