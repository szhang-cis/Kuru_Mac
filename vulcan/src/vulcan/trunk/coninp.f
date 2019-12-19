      SUBROUTINE CONINP
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
C     NHOUR              'HOURG
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
      INCLUDE 'inpo_om.f'
C
      PARAMETER (MCOMM=8)
C
      CHARACTER*5 COMMD(MCOMM)
C
      DATA COMMD/'POSTP','RENUM','SMOOT','SOLVE',
     .           'CPULI','DATAL','HOURG','END_C'/
      KPARAM=2
C
C**** READ CONTROL PARAMETERS
C
      NPRIN=0
      ITAPE=LUDAT
      CALL LISTEN('CONINP',NPRIN,ITAPE)
      IF(WORDS(1).NE.'CONTR') GO TO 2000
 1000 CALL LISTEN('CONINP',NPRIN,ITAPE)
C
C**** IDENTIFY COMMAND
C
      DO ICOMM=1,MCOMM
       IF(WORDS(1).EQ.COMMD(ICOMM)) GOTO 100
      ENDDO
      GO TO 2000 
C
C**** EXECUTE APROPRIATE COMMAND
C
  100 CONTINUE
      GOTO (1,2,3,4,5,6,7,8), ICOMM
C
    1 KPOST=1
      GO TO 1000
C
    2 KRENU=1
      GO TO 1000
C
    3 IF(WORDS(2).EQ.'DISCR') THEN
       KSGAU=-1
       ISMO1=INT(PARAM(1))                                ! see elm030.f
       ISMO2=INT(PARAM(2))
      ENDIF
      IF(WORDS(2).EQ.'LOCAL') KSGAU=-2
      GO TO 1000
C
    4 IF(WORDS(2).EQ.'PROFI') KSOLV=0
      IF(WORDS(2).EQ.'FRONT') KSOLV=1
      IF(WORDS(2).EQ.'PCGRA') KSOLV=2
      IF(WORDS(2).EQ.'GMRES') KSOLV=3
      IF(WORDS(2).EQ.'PARDI') KSOLV=4
C
C**** SOLVER'S PARAMETERS (see setdat for its default values)
C
C     KSOLV=     skyline, frontal, pcg & gmres
C     KSYMM=     skyline (param(1)) & frontal (param(1))
C     NWIDT=     skyline (param(2)) & pcg (param(3))
C     MITCG=     skyline (param(3)) & pcg (param(1))
C     NBUFA=     frontal (param(2))
C     TOLCG=     skyline (param(4)) & pcg (param(2))
C     TOLC1=     skyline (param(5)) & pcg (param(4))
C     NPRIR=     pcg (words(3))
C     MITGM=     gmres (param(1))
C     TOLGM=     gmres (param(2))
C     MKRYL=     gmres (param(3))
C     IPGMR=     gmres (words(3))
C
      IF(KSOLV.EQ.0) THEN         ! Skyline solver (direct & semidirect)
       KSYMM=INT(PARAM(1))
       IF(PARAM(2).NE.0.0D0) NWIDT=INT(PARAM(2))
       NWIDT=NWIDT-1              ! Discount Diagonal SKYLINE 
       IF(PARAM(3).NE.0.0D0) MITCG=INT(PARAM(3))
       IF(PARAM(4).NE.0.0D0) TOLCG=PARAM(4)
       IF(PARAM(5).NE.0.0D0) TOLC1=PARAM(5)
      ENDIF
C
      IF(KSOLV.EQ.1) THEN         ! Frontal Solver 
       KSYMM=INT(PARAM(1))
       NBUFA=INT(PARAM(2))
      ENDIF
C
      IF(KSOLV.EQ.2) THEN         ! PCG Solver
       MITCG=1000
       IF(PARAM(1).NE.0.0D0) MITCG=INT(PARAM(1))
       IF(PARAM(2).NE.0.0D0) TOLCG=PARAM(2)
       NWIDT=0
       IF(PARAM(3).NE.0.0D0) NWIDT=-1
       IF(PARAM(4).NE.0.0D0) TOLC1=PARAM(4)
       NPRIR=0
       IF(WORDS(3).EQ.'PRINT') NPRIR=1
      ENDIF
C
      IF(KSOLV.EQ.3) THEN         ! GMRES
       KSYMM=0                    ! unsymmetric
       MITGM=1000
       IF(PARAM(1).NE.0.0D0) MITGM=INT(PARAM(1))
       IF(PARAM(2).NE.0.0D0) TOLGM=PARAM(2)
       MKRYL=10
       IF(PARAM(3).NE.0.0D0) MKRYL=INT(PARAM(3))
       IPGMR=0
       IF(WORDS(3).EQ.'PRINT') IPGMR=1
C
       IF(NMEMO7.EQ.1)
     .  CALL RUNEND('NMEMO7=1 ARE NOT COMPATIBLE WITH GMRES SOLVER')
      ENDIF
C
      IF(KSOLV.EQ.4) THEN         ! PARDISO solver (direct & semidirect)  !no cambio aparente
       KSYMM=INT(PARAM(1))
       IF(PARAM(2).NE.0.0D0) THEN
        NTHSOM=INT(PARAM(2))
       ELSE
        NTHSOM=1
       ENDIF
      ENDIF
      GOTO 1000
C
    5 IF(PARAM(1).NE.0.0D0) TLIMT=PARAM(1)
      GO TO 1000
C
    6 NPARA=INT(PARAM(1))
      IF(NPARA.GE.0.AND.NPARA.LT.NBLIM) NBLIM=NPARA
      GO TO 1000
C
    7 NHOUR=1
      KELAS=0
      HPARA=1.0D0
      IF(WORDS(2).EQ.'ELAST') KELAS=1
      IF(PARAM(1).NE.0.0D0)   HPARA=PARAM(1)
      GO TO 1000
C
    8 CONTINUE
C
      RETURN 
 2000 CALL RUNEND('CONINP: ERROR IN CONTRO DATA BLOCK')
      END
