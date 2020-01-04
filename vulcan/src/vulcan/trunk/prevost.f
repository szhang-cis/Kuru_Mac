      SUBROUTINE PREVOST(DISTOT,ELDATT,ELPRET,ELVART,ELMATT,HEADST,
     .                   LNODST,MATNOT,PROELT,PROPST,TEMPIT,WORK1T)
C***********************************************************************
C
C**** THIS ROUTINE READS THE INITIAL STATE
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
      COMMON/PRFILET/JTAPET,IFLAGT
C
      DIMENSION MATNOT(NELEMT),        LNODST(NNODET,NELEMT),
     .          PROELT(NPRELT,NGRUPT), PROPST(NPROPT,NMATST),
     .          ELDATT(NDATAT),        ELPRET(NPREVT),
     .          ELVART(NSTATT),        ELMATT(NMATXT)
      DIMENSION DISTOT(NTOTVT,*),      HEADST(NPOINT,*),
     .          TEMPIT(NPOINT)
      DIMENSION WORK1T(*)
C
      PARAMETER   (MCOMMT=7)
      CHARACTER*5 COMMDT(MCOMMT)
C
      SAVE IFLA1T,IFLA2T,ISTOPT
      DATA IFLA1T,IFLA2T,ISTOPT/0,0,0/
C
      DATA COMMDT/'GAUSS','WATER','TEMPI','TEMPE','VELOC','ACCEL',
     .            'END_I'/
C
      CALL CPUTIMT(TIME1T)
      WRITE(LUREST,900)
C
C**** READ INITIAL_DATA CARD
C
      NPRINT=0
      ITAPET=LUDATT
      CALL LISTENT('PREVOST',NPRINT,ITAPET)
      IF(WORDST(1).NE.'INITI') GO TO 1000
C
C**** READ NEW COMMAND
C
  100 CALL LISTENT('PREVOST',NPRINT,ITAPET)
C
C**** IDENTIFY COMMAND
C
      DO ICOMMT=1,MCOMMT
        IF(WORDST(1).EQ.COMMDT(ICOMMT)) GO TO 20
      END DO
      GO TO 1000 
C
C**** EXECUTE APPROPRIATE COMMAND
C
   20 JTAPET=LUDATT
      IF(PARAMT(1).NE.0.0) JTAPET=INT(PARAMT(1))
C
      GO TO (1,2,3,4,4,4,5), ICOMMT
C
C**** READ THE INITIAL STATE VARIABLES AT THE GAUSSIAN POINTS
C
    1 IFLAGT=0
      IF(NMEMO.EQ.0)
     . CALL RUNENDT('ERROR (PREVOST): MEMORY IS NOT DIMENSIONED FOR 
     .               INITIAL GAUSSIAN VARIABLES')
      IF(WORDST(2).EQ.'STRAI') IFLAGT=1
      IF(WORDST(2).EQ.'STRES') IFLAGT=2
      IF(WORDST(2).EQ.'PLAST') IFLAGT=3
      IF(IFLAGT.NE.0)
     .  CALL PREGAUT(ELDATT,ELPRET,ELVART,ELMATT,LNODST,MATNOT,PROELT,
     .               PROPST)
      GO TO 100
C
C**** READ THE INITIAL PORE PRESSURE AT NODES
C
    2 CONTINUE
      IF(KPORET.NE.2) KPORET=1
      GO TO 100
C
C**** READ THE TEMPERATURE DATA AT NODES
C
    3 CONTINUE
      KTEMPT=1
      GO TO 100
C
C**** READ THE INITIAL TEMPERATURE AT NODES
C     (VELOCITIES OR ACCELERATIONS AT NODES)
C
    4 IF(WORDST(1).EQ.'TEMPE') NCOMPX=1
      IF(WORDST(1).EQ.'VELOC') NCOMPX=2
      IF(WORDST(1).EQ.'ACCEL') NCOMPX=3
      CALL PREDIST(DISTOT(1,NCOMPX),NCOMPX)
      DO ITOTVT=1,NTOTVT
       TEMPIT(ITOTVT)=DISTOT(ITOTVT,1)
      ENDDO
      GO TO 100
C
C**** RETURN
C
    5 CALL CPUTIMT(TIME2T)
      CPUDAT=CPUDAT+(TIME2T-TIME1T)
      RETURN
C
 1000 CALL RUNENDT('PREVOST:ERROR IN INITIAL DATA BLOCK')
  900 FORMAT(1H1,///,10X,'ECHO OF INITIAL STATE DATA FOLLOWS :',/
     .               10X,'=================================='/)
      END