      SUBROUTINE RSTAR1S(IFLAGS)
C***********************************************************************
C
C**** THIS ROUTINE READS OR WRITES FROM FILE 1 RECORD 2 THE CONTROLLING
C     PARAMETERS SETTED UP IN SUB.CONTRO
C
C   - IFLAGS = 1 - NEW RUN    : WRITE INTO RESTART FILE 
C   - IFLAGS = 2 - RESTART RUN: READ  FROM RESTART FILE
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
      INCLUDE 'auxl_oms.f'
C
      CALL CPUTIMS(TIME1S)
C
      NREC0S = 2
C
      GOTO (100,200),IFLAGS
C***********************************************************************
C
C*** N E W   R U N  *** . . . . . . . . . . . . . . . . . . . . ..
C
C***********************************************************************
 100  CONTINUE
      IF(NFURESS.EQ.2) WRITE(LURSTS,REC=NREC0S)
C---------------------------------------------------------------
C....COMMON/DATBASSA/
     .              IDATPS,IDATCS,LENRCS,NLENCS,NLENPS,NRECGS,NRECPS,
C....COMMON/DIMENNSA/
     .              NDIMES,NELEMS,NFUNCS,NGRUPS,NHLODS,NHISTS,NMATSS,
     .              NPOINS,NPRELS,NPROPS,NSTR1S,NTOTVS 
C
      IF(NFURESS.EQ.2) WRITE(LURSTS,REC=NREC0S+1)
C....COMMON/ELDATASA/
     .              IDATAS,IPREVS,ISTATS,IMATXS,
C....COMMON/ELEMNTSA/
     .              NDOFCS,NDOFNS,NEVABS,NEVACS,NGAUSS,NKONDS,NKOVAS,
     .              NNODES,NDATAS,NPREVS,NSTATS,NMATXS,
C....COMMON/HOURGLS/
     .              NHOURS,KELASS,HPARAS,
C....COMMON/PROBLMSA/
     .              KDYNAS,KPORES,KPOSTS,KPROBS,KSGAUS,KSMUSS,KTEMPS,
     .              LARGES,
C....COMMON/SOLVERS/
     .              KRENUS,KSOLVS,KSYMMS,NWIDTS,MITCGS,TOLCGS,NBUFAS,
     .              NPRIRS,
     .              TOLC1S,
C....COMMON/TITLESSA/
     .              TITLES,SUBTIS
C---------------------------------------------------------------
C
      GOTO 300
C***********************************************************************
C
C*** R E S T A R T   R U N  *** . . . . . . . . . . . . . . . . . . . 
C
C***********************************************************************
  200 READ(LURSTS,REC=NREC0S) 
C---------------------------------------------------------------
C....COMMON/DATBASSA/
     .              IDATPS,IDATCS,LENRCS,NLENCS,NLENPS,NRECGS,NRECPS,
C....COMMON/DIMENNSA/
     .              NDIMES,NELEMS,NFUNCS,NGRUPS,NHLODS,NHISTS,NMATSS,
     .              NPOINS,NPRELS,NPROPS,NSTR1S,NTOTVS
C
      READ(LURSTS,REC=NREC0S+1)
C....COMMON/ELDATASA/
     .              IDATAS,IPREVS,ISTATS,IMATXS,
C....COMMON/ELEMNTSA/
     .              NDOFCS,NDOFNS,NEVABS,NEVACS,NGAUSS,NKONDS,NKOVAS,
     .              NNODES,NDATAS,NPREVS,NSTATS,NMATXS,
C....COMMON/HOURGLS/
     .              NHOURS,KELASS,HPARAS,
C....COMMON/PROBLMSA/
     .              KDYNAS,KPORES,KPOSTS,KPROBS,KSGAUS,KSMUSS,KTEMPS,
     .              LARGES,
C....COMMON/SOLVERS/
     .              KRENUS,KSOLVS,KSYMMS,NWIDTS,MITCGS,TOLCGS,NBUFAS,
     .              NPRIRS,
     .              TOLC1S,
C....COMMON/TITLESSA/
     .              TITLES,SUBTIS
C---------------------------------------------------------------
C
C.....INITIALIZE POINTERS FOR POSITION IN DOUBLE PRECISION WORDS
C
      DO INDEXS=1,11
       IDATPS(INDEXS,5)=-1
      ENDDO
C
 300  CONTINUE
C
      CALL CPUTIMS(TIME2S)
      CPURSS=CPURSS+(TIME2S-TIME1S)
C
      RETURN
      END
