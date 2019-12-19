      SUBROUTINE CONINPS
C***********************************************************************
C
C**** THIS ROUTINE REDEFINES THE DEFAULT CONTROL PARAMETERS FOR THIS RUN
C
C     FOR CHANGING THE DEFAULT PARAMETERS
C
C     PARAMETERS         CONTROLLING WORDS
C
C     KPOST              'POSTP'
C     KRENU              'RENUM'
C     KSGAU              'SMOOT'
C     KSOLV              'SOLVE'
C     NBLIM              'DATAL'
C     TLIMT              'CPULI'
C     NHOUR              'HOURG'
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
      INCLUDE 'inpo_oms.f'
C
      PARAMETER (MCOMMS=8)
      CHARACTER*5 COMMDS(MCOMMS)
      DATA COMMDS/'POSTP','RENUM','SMOOT','SOLVE',
     .            'CPULI','DATAL','HOURG','END_C'/
C
C**** READ CONTROL PARAMETERS
C
      NPRINT=0
      ITAPET=LUDATS
      CALL LISTENS('CONINPS',NPRINT,ITAPET)
      IF(WORDSS(1).NE.'CONTR') GO TO 2000
 1000 CALL LISTENS('CONINPS',NPRINT,ITAPET)
C
C**** IDENTIFY COMMAND
C
      DO ICOMMS=1,MCOMMS
       IF(WORDSS(1).EQ.COMMDS(ICOMMS)) GOTO 100
      ENDDO
      GO TO 2000 
C
C**** EXECUTE APROPRIATE COMMAND
C
  100 CONTINUE
      GOTO (1,2,3,4,5,6,7,8), ICOMMS
C
    1 KPOSTS=1
      GO TO 1000
C
    2 KRENUS=1
      GO TO 1000
C
    3 IF(WORDSS(2).EQ.'DISCR') THEN
       KSGAUS=-1
       ISMO1S=INT(PARAMS(1))                             ! see elm005s.f
       ISMO2S=INT(PARAMS(2))
      ENDIF
      IF(WORDSS(2).EQ.'LOCAL') KSGAUS=-2
      CALL RUNENDS('ERROR IN coninps.f: SMOOTHING NOT IMPLEMENTED YET')
      GO TO 1000
C
    4 IF(WORDSS(2).EQ.'PROFI') KSOLVS=0
      IF(WORDSS(2).EQ.'FRONT') KSOLVS=1
      IF(WORDSS(2).EQ.'PCGRA') KSOLVS=2
      IF(WORDSS(2).EQ.'GMRES') KSOLVS=3
      IF(WORDSS(2).EQ.'EXPLI') KSOLVS=4
C
C**** SOLVER'S PARAMETERS (see setdatt for its default values)
C
C     KSOLVS=     skyline, frontal, pcg & gmres
C     KSYMMS=     skyline (params(1)) & frontal (params(1))
C     NWIDTS=     skyline (params(2)) & pcg (params(3))
C     MITCGS=     skyline (params(3)) & pcg (params(1))
C     NBUFAS=     frontal (params(2))
C     TOLCGS=     skyline (params(4)) & pcg (params(2))
C     TOLC1S=     skyline (params(5)) & pcg (params(4))
C     NPRIRS=     pcg (wordss(3))
C
      IF(KSOLVS.EQ.0) THEN        ! Skyline solver (direct & semidirect)
       call runends('error: skyline solver not implemented')
       KSYMMS=INT(PARAMS(1))
       IF(PARAMS(2).NE.0.0) NWIDTS=INT(PARAMS(2))
       NWIDTS=NWIDTS-1            ! Discount Diagonal SKYLINE 
       IF(PARAMS(3).NE.0.0) MITCGS=INT(PARAMS(3))
       IF(PARAMS(4).NE.0.0) TOLCGS=PARAMS(4)
       IF(PARAMS(5).NE.0.0) TOLC1S=PARAMS(5)
      ENDIF
C
      IF(KSOLVS.EQ.1) THEN        ! Frontal Solver 
       KSYMMS=INT(PARAMS(1))
       NBUFAS=INT(PARAMS(2))
      ENDIF
C
      IF(KSOLVS.EQ.2) THEN        ! PCG Solver
       call runends('error: PCG solver not implemented')
       MITCGS=1000
       IF(PARAMS(1).NE.0.0) MITCGS=INT(PARAMS(1))
       IF(PARAMS(2).NE.0.0) TOLCGS=PARAMS(2)
       NWIDTS=0
       IF(PARAMS(3).NE.0.0) NWIDTS=-1
       IF(PARAMS(4).NE.0.0) TOLC1S=PARAMS(4)
       NPRIRS=0
       IF(WORDSS(3).EQ.'PRINT') NPRIRS=1
      ENDIF
C
      IF(KSOLVS.EQ.3) THEN        ! GMRES
       CALL RUNENDS('GMRES SOLVER NOT IMPLEMENTED')
      ENDIF
C
      IF(KSOLVS.EQ.4) THEN        ! Explicit solution (no solver)
       call runends('error: explicit solver not implemented')
       CALL RUNMENS('WARNING: EXPLICIT SOLUTION WILL BE USED            
     .                              ')
       CALL RUNMENS('CHECK INTEGRATION RULE TO OBTAIN A LUMPED MATRIX  
     .                              ')
       CALL RUNMENS('IF A FULL MATRIX IS USED, PUT NCETAT=1 (NOT RECOMME
     .NDED)                         ')
C
       NMEMO7S=1
       NMEMO6S=1
      ENDIF
C
      GOTO 1000
C
    5 IF(PARAMS(1).NE.0.0) TLIMTS=PARAMS(1)
      GO TO 1000
C
    6 NPARAS=INT(PARAMS(1))
      IF(NPARAS.GE.0.AND.NPARAS.LT.NBLIMS) NBLIMS=NPARAS
      GO TO 1000
C
    7 NHOURS=1
      KELASS=0
      HPARAS=1.0
      IF(WORDSS(2).EQ.'ELAST') KELASS=1
      IF(PARAMS(1).NE.0.0)     HPARAS=PARAMS(1)
      GO TO 1000
C
    8 CONTINUE
C
      RETURN 
 2000 CALL RUNENDS('CONINPS:ERROR IN CONTRO DATA BLOCK')
      END
