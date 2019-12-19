      SUBROUTINE CONINPT
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
      INCLUDE 'addi_omt.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'
C

      common/bufapp/nbufp
cmem

C
      PARAMETER (MCOMMT=8)
      CHARACTER*5 COMMDT(MCOMMT)
      DATA COMMDT/'POSTP','RENUM','SMOOT','SOLVE',
     .            'CPULI','DATAL','HOURG','END_C'/
C
C**** READ CONTROL PARAMETERS
C
      NPRINT=0
      ITAPET=LUDATT
      CALL LISTENT('CONINPT',NPRINT,ITAPET)
      IF(WORDST(1).NE.'CONTR') GO TO 2000
 1000 CALL LISTENT('CONINPT',NPRINT,ITAPET)
C
C**** IDENTIFY COMMAND
C
      DO ICOMMT=1,MCOMMT
       IF(WORDST(1).EQ.COMMDT(ICOMMT)) GOTO 100
      ENDDO
      GO TO 2000 
C
C**** EXECUTE APROPRIATE COMMAND
C
  100 CONTINUE
      GOTO (1,2,3,4,5,6,7,8), ICOMMT
C
    1 KPOSTT=1
      GO TO 1000
C
    2 KRENUT=1
      GO TO 1000
C
    3 IF(WORDST(2).EQ.'DISCR') THEN
       KSGAUT=-1
       ISMO1T=INT(PARAMT(1))                             ! see elm005t.f
       ISMO2T=INT(PARAMT(2))
      ENDIF
      IF(WORDST(2).EQ.'LOCAL') KSGAUT=-2
      GO TO 1000
C
    4 IF(WORDST(2).EQ.'PROFI') KSOLVT=0
      IF(WORDST(2).EQ.'FRONT') KSOLVT=1
      IF(WORDST(2).EQ.'PCGRA') KSOLVT=2
      IF(WORDST(2).EQ.'GMRES') KSOLVT=3
      IF(WORDST(2).EQ.'EXPLI') KSOLVT=4
C
C**** SOLVER'S PARAMETERS (see setdatt for its default values)
C
C     KSOLVT=     skyline, frontal, pcg, gmres & explicit
C     KSYMMT=     skyline (paramt(1)) & frontal (paramt(1))
C     NWIDTT=     skyline (paramt(2)) & pcg (paramt(3))
C     MITCGT=     skyline (paramt(3)) & pcg (paramt(1))
C     NBUFAT=     frontal (paramt(2))
C     TOLCGT=     skyline (paramt(4)) & pcg (paramt(2))
C     TOLC1T=     skyline (paramt(5)) & pcg (paramt(4))
C     NPRIRT=     pcg (wordst(3))
C     MITGMT=     gmres (paramt(1))
C     TOLGMT=     gmres (paramt(2))
C     MKRYLT=     gmres (paramt(3))
C     IPGMRT=     gmres (wordst(3))
C
      IF(KSOLVT.EQ.0) THEN        ! Skyline solver (direct & semidirect)
       KSYMMT=INT(PARAMT(1))
       IF(PARAMT(2).NE.0.0D0) NWIDTT=INT(PARAMT(2))
       NWIDTT=NWIDTT-1            ! Discount Diagonal SKYLINE 
       IF(PARAMT(3).NE.0.0D0) MITCGT=INT(PARAMT(3))
       IF(PARAMT(4).NE.0.0D0) TOLCGT=PARAMT(4)
       IF(PARAMT(5).NE.0.0D0) TOLC1T=PARAMT(5)
      ENDIF
C
      IF(KSOLVT.EQ.1) THEN        ! Frontal Solver 
       KSYMMT=INT(PARAMT(1))
       NBUFAT=INT(PARAMT(2))

       nbufp=int(paramt(3))
cmem

      ENDIF
C
      IF(KSOLVT.EQ.2) THEN        ! PCG Solver
       MITCGT=1000
       IF(PARAMT(1).NE.0.0D0) MITCGT=INT(PARAMT(1))
       IF(PARAMT(2).NE.0.0D0) TOLCGT=PARAMT(2)
       NWIDTT=0
       IF(PARAMT(3).NE.0.0D0) NWIDTT=-1
       IF(PARAMT(4).NE.0.0D0) TOLC1T=PARAMT(4)
       NPRIRT=0
       IF(WORDST(3).EQ.'PRINT') NPRIRT=1
      ENDIF
C
      IF(KSOLVT.EQ.3) THEN        ! GMRES
       KSYMMT=0                   ! unsymmetric
       MITGMT=1000
       IF(PARAMT(1).NE.0.0D0) MITGMT=INT(PARAMT(1))
       IF(PARAMT(2).NE.0.0D0) TOLGMT=PARAMT(2)
       MKRYLT=10
       IF(PARAMT(3).NE.0.0D0) MKRYLT=INT(PARAMT(3))
       IPGMRT=0
       IF(WORDST(3).EQ.'PRINT') IPGMRT=1
C
       IF(NMEMO7.EQ.1)
     .  CALL RUNENDT('NMEMO7=1 ARE NOT COMPATIBLE WITH GMRES SOLVER')
      ENDIF
C
      IF(KSOLVT.EQ.4) THEN        ! Explicit solution (no solver)
       CALL RUNMENT('WARNING: EXPLICIT SOLUTION WILL BE USED            
     .                              ')
       CALL RUNMENT('CHECK INTEGRATION RULE TO OBTAIN A LUMPED MATRIX  
     .                              ')
       CALL RUNMENT('IF A FULL MATRIX IS USED, PUT NCETAT=1 (NOT RECOMME
     .NDED)                         ')
C
       NMEMO7=1
       NMEMO6=1
      ENDIF
C
      GOTO 1000
C
    5 IF(PARAMT(1).NE.0.0D0) TLIMTT=PARAMT(1)
      GO TO 1000
C
    6 NPARAT=INT(PARAMT(1))
      IF(NPARAT.GE.0.AND.NPARAT.LT.NBLIMT) NBLIMT=NPARAT
      GO TO 1000
C
    7 NHOURT=1
      KELAST=0
      HPARAT=1.0D0
      IF(WORDST(2).EQ.'ELAST') KELAST=1
      IF(PARAMT(1).NE.0.0D0)   HPARAT=PARAMT(1)
      GO TO 1000
C
    8 CONTINUE
C
      RETURN 
 2000 CALL RUNENDT('CONINPT:ERROR IN CONTRO DATA BLOCK')
      END
