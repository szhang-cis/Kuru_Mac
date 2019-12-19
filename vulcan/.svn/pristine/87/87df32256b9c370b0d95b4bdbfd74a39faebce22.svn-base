      SUBROUTINE RSTAR1T(IFLAGT)
C***********************************************************************
C
C**** THIS ROUTINE READS OR WRITES FROM FILE 1 RECORD 2 THE CONTROLLING
C     PARAMETERS SETTED UP IN SUB.CONTRO
C
C   - IFLAG = 1 - NEW RUN    : WRITE INTO RESTART FILE 
C   - IFLAG = 2 - RESTART RUN: READ  FROM RESTART FILE
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
      INCLUDE 'auxl_omt.f'
C
      CALL CPUTIMT(TIME1T)
C
      NREC0T = 2
C
      GOTO (100,200),IFLAGT
C***********************************************************************
C
C*** N E W   R U N  *** . . . . . . . . . . . . . . . . . . . . ..
C
C***********************************************************************
 100  CONTINUE
      IF(NFURES.EQ.2) WRITE(LURSTT,REC=NREC0T)
C---------------------------------------------------------------
C....COMMON/DATBAST/
     .              IDATPT,IDATCT,LENRCT,NLENCT,NLENPT,NRECGT,NRECPT,
C....COMMON/DIMENNT/
     .              NDIMET,NELEMT,NFUNCT,NGRUPT,NHLODT,NHISTT,NMATST,
     .              NPOINT,NPRELT,NPROPT,NSTR1T,NTOTVT 
C
      IF(NFURES.EQ.2) WRITE(LURSTT,REC=NREC0T+1)
C....COMMON/ELDATAT/
     .              IDATAT,IPREVT,ISTATT,IMATXT,
C....COMMON/ELEMNTT/
     .              NDOFCT,NDOFNT,NEVABT,NEVACT,NGAUST,NKONDT,NKOVAT,
     .              NNODET,NDATAT,NPREVT,NSTATT,NMATXT,
C....COMMON/HOURGLT/
     .              NHOURT,KELAST,HPARAT,
C....COMMON/PROBLMT/
     .              KDYNAT,KPORET,KPOSTT,KPROBT,KSGAUT,KSMUST,KTEMPT,
     .              LARGET,
C....COMMON/SOLVERT/
     .              KRENUT,KSOLVT,KSYMMT,NWIDTT,MITCGT,TOLCGT,NBUFAT,
     .              NPRIRT,
     .              TOLC1T,
     .              MITGMT,MKRYLT,IPGMRT,
     .              TOLGMT,
C....COMMON/TITLEST/
     .              TITLET,SUBTIT
C---------------------------------------------------------------
C
      GOTO 300
C***********************************************************************
C
C*** R E S T A R T   R U N  *** . . . . . . . . . . . . . . . . . . . 
C
C***********************************************************************
  200 READ(LURSTT,REC=NREC0T) 
C---------------------------------------------------------------
C....COMMON/DATBAST/
     .              IDATPT,IDATCT,LENRCT,NLENCT,NLENPT,NRECGT,NRECPT,
C....COMMON/DIMENNT/
     .              NDIMET,NELEMT,NFUNCT,NGRUPT,NHLODT,NHISTT,NMATST,
     .              NPOINT,NPRELT,NPROPT,NSTR1T,NTOTVT
C
      READ(LURSTT,REC=NREC0T+1)
C....COMMON/ELDATAT/
     .              IDATAT,IPREVT,ISTATT,IMATXT,
C....COMMON/ELEMNTT/
     .              NDOFCT,NDOFNT,NEVABT,NEVACT,NGAUST,NKONDT,NKOVAT,
     .              NNODET,NDATAT,NPREVT,NSTATT,NMATXT,
C....COMMON/HOURGLT/
     .              NHOURT,KELAST,HPARAT,
C....COMMON/PROBLMT/
     .              KDYNAT,KPORET,KPOSTT,KPROBT,KSGAUT,KSMUST,KTEMPT,
     .              LARGET,
C....COMMON/SOLVERT/
     .              KRENUT,KSOLVT,KSYMMT,NWIDTT,MITCGT,TOLCGT,NBUFAT,
     .              NPRIRT,
     .              TOLC1T,
C....COMMON/TITLEST/
     .              TITLET,SUBTIT
C---------------------------------------------------------------
C
C.....INITIALIZE POINTERS FOR POSITION IN DOUBLE PRECISION WORDS
C
      DO INDEXT=1,11
       IDATPT(INDEXT,5)=-1
      ENDDO
C
 300  CONTINUE
C
      CALL CPUTIMT(TIME2T)
      CPURST=CPURST+(TIME2T-TIME1T)
C
      RETURN
      END
