      SUBROUTINE RSTAR1(IFLAG)
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
      INCLUDE 'addi_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'auxl_om.f'
C
      CALL CPUTIM(TIME1)
C
      NREC0 = 2
C
      GOTO (100,200),IFLAG
C***********************************************************************
C
C*** N E W   R U N  *** . . . . . . . . . . . . . . . . . . . . ..
C
C***********************************************************************
 100  CONTINUE
      IF(NFURESM.EQ.2) WRITE(LURST,REC=NREC0)
C---------------------------------------------------------------
C....COMMON/DATBAS/
     .             IDATP,IDATC,LENRC,NLENC,NLENP,NRECG,NRECP,
C....COMMON/DIMENN/
     .             NDIME,NELEM,NFUNC,NGRUP,NHLOD,NHIST,NMATS,
     .             NPOIN,NPREL,NPROP,NSTR1,NTOTV 
C
      IF(NFURESM.EQ.2) WRITE(LURST,REC=NREC0+1)
C....COMMON/ELDATA/
     .             IDATA,IPREV,ISTAT,IMATX,
C....COMMON/ELEMNT/
     .             NDOFC,NDOFN,NEVAB,NEVAC,NGAUS,NKOND,NKOVA,
     .             NNODE,NDATA,NPREV,NSTAT,NMATX,
C....COMMON/HOURGL/
     .             NHOUR,KELAS,HPARA,
C....COMMON/PROBLM/
     .             KDYNA,KPORE,KPOST,KPROB,KSGAU,KSMUS,KTEMP,
     .             LARGE,
C....COMMON/SOLVER/
     .             KRENU,KSOLV,KSYMM,NWIDT,MITCG,TOLCG,NBUFA,NPRIR,
     .             TOLC1,
C....COMMON/TITLES/
     .             TITLE,SUBTI
C---------------------------------------------------------------
C
      GOTO 300
C***********************************************************************
C
C*** R E S T A R T   R U N  *** . . . . . . . . . . . . . . . . . . . 
C
C***********************************************************************
  200 READ(LURST,REC=NREC0) 
C---------------------------------------------------------------
C....COMMON/DATBAS/
     .             IDATP,IDATC,LENRC,NLENC,NLENP,NRECG,NRECP,
C....COMMON/DIMENN/
     .             NDIME,NELEM,NFUNC,NGRUP,NHLOD,NHIST,NMATS,
     .             NPOIN,NPREL,NPROP,NSTR1,NTOTV
C
      READ(LURST,REC=NREC0+1)
C....COMMON/ELDATA/
     .             IDATA,IPREV,ISTAT,IMATX,
C....COMMON/ELEMNT/
     .             NDOFC,NDOFN,NEVAB,NEVAC,NGAUS,NKOND,NKOVA,
     .             NNODE,NDATA,NPREV,NSTAT,NMATX,
C....COMMON/HOURGL/
     .             NHOUR,KELAS,HPARA,
C....COMMON/PROBLM/
     .             KDYNA,KPORE,KPOST,KPROB,KSGAU,KSMUS,KTEMP,
     .             LARGE,
C....COMMON/SOLVER/
     .             KRENU,KSOLV,KSYMM,NWIDT,MITCG,TOLCG,NBUFA,NPRIR,
     .             TOLC1,
     .             MITGM,MKRYL,IPGMR,
     .             TOLGM,
C....COMMON/TITLES/
     .             TITLE,SUBTI
C---------------------------------------------------------------
C
C.....INITIALIZE POINTERS FOR POSITION IN DOUBLE PRECISION WORDS
C
      DO INDEX=1,11
       IDATP(INDEX,5)=-1
      ENDDO
C
 300  CONTINUE
C
      CALL CPUTIM(TIME2)
      CPURS=CPURS+(TIME2-TIME1)
C
      RETURN
      END
