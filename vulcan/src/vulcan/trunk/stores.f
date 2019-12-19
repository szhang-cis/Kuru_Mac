      SUBROUTINE STORES(DISPR,DISTO,ELPRE,ELVAR,HEADS,HTLOD,IFFIX,
     .                  REFOR,RLOAD,TLOAD,LACTI,ILATE)
C***********************************************************************
C
C**** THIS ROUTINE STORES THE RESULTS OF THE CURRENT TIME STEP
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
C
      DIMENSION HTLOD(NHLOD,NSUBF,NFUNC), IFFIX(NTOTV,2),
     .          RLOAD(NTOTV),             DISPR(NTOTV,NDISR),
     .          DISTO(NTOTV,NDISO),       HEADS(NPOIN,4), 
     .          REFOR(NTOTV),             TLOAD(NTOTV,2),
     .          ELPRE(NPREV),             ELVAR(NSTAT)
      DIMENSION LACTI(NELEM)
C
C=======================================================================
C
C**** RESTART FILE
C
C=======================================================================
C
      IF(NFURESM.EQ.1) THEN
       CALL RSTAR4(DISPR,DISTO,ELPRE,ELVAR,TLOAD,LACTI)
C
       WRITE(LUPRI,*)
       WRITE(LUPRI,*) 'FUTURE RESTART INFORMATION'
       WRITE(LUPRI,*) 'RESULTS STORED FOR TIME:',TTIME
       WRITE(LUPRI,*)
       CLOSE(LUPRI)
       OPEN(UNIT=LUPRI,FILE=CF,STATUS='OLD',ACCESS='APPEND')
      ENDIF                      ! nfuresm.eq.1
C
C**** DECIDE ON SAVING THESE RESULTS
C
      IF(NFURESM.EQ.2) THEN
C
       ISAVE=0
       IF(MOD(ISTEP,NOUTP(1)).EQ.0) ISAVE=1
       IF(ILATE.EQ. 1)              ISAVE=1
       IF(KSAVE.EQ.-1)              ISAVE=0
C
C**** DUMP TO DATA BASE & SET UP APPROPRIATE COUNTERS
C
       IF(ISAVE.EQ.1) THEN
        IF(KTSTE.EQ.0) KTSTE=1
        IF(KSAVE.EQ.1) KTSTE=KTSTE+1
        CALL RSTAR3(DISPR,DISTO,ELPRE,ELVAR,HEADS,HTLOD,IFFIX,
     .              REFOR,RLOAD,TLOAD,    1)
        IF(KSAVE.EQ.1) NRECC=NRECC+NLENC
       ENDIF
C
      ENDIF                      ! nfuresm.eq.2
C
C=======================================================================
C
C**** PROCESS FILE
C
C=======================================================================
C
C**** INTERCHANGE CURRENT AND CONVERGED VALUES OF ELPRE & ELVAR
C
      DO INDEX=1,5
       IDUMM=IDATP(2,INDEX)
       IDATP(2,INDEX)=IDATP(4,INDEX)
       IDATP(4,INDEX)=IDUMM
       IDUMM=IDATP(3,INDEX)
       IDATP(3,INDEX)=IDATP(5,INDEX)
       IDATP(5,INDEX)=IDUMM
      ENDDO
C
C**** WRITE DISTO TO PROCESS FILE
C
      CALL DATBAS(DISTO,  12,    1)
C
C=======================================================================
C
C**** RESULTS FILE
C
C=======================================================================
C
      WRITE(LUPRI,900)
  900 FORMAT(50X,'>>>>>>> STEP WRITTEN UP',/)
ctm
ctm   cambio de sistema operativo del CONVEX (febrero de 1993)
ctm
c      CLOSE(LURES)
c      OPEN(UNIT=LURES,STATUS='OLD',ACCESS='APPEND')
ctm
C
C**** INITIALISE FLAGS IREST,ISKIP
C
      IREST=0
      ISKIP=0
C
      RETURN
      END
