      SUBROUTINE EXPLICT(DISIT,IFFIX,REFOR,LNODS,
     .                   NDOFN,NELEM,NEVAB,NNODE,NPOIN,NTOTV,
     .                   INRHS,IPASS,KRESL,
     .                   CSTIF,DISIM,FOREL,
     .                   CSTII,CSTIT,
     .                   NPREL,NGRUP,NPROP,NMATS,
     .                   NDATA,NPREV,NSTAT,NMATX,NDIME,NTOTVM,NFPCH,
     .                   MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .                   ELVAR,ELMAT,WORK1,DISTO,COORD,
     .                   ADVEL,TEMPI,PREAS,TGAPS,DISPL,FPCHA,
     .                   NACTI,LACTI)
C***********************************************************************
C
C**** THIS ROUTINE SOLVES THE SET OF LINEAR EQUATIONS USING THE
C     "EXPLICIT_SOLUTION" PROCEDURE
C
C     EXPLIC
C            - PCGASS 
C            - PCGTRI
C            - PCGITE  - PCGPRE - PCGMUL
C                      - SKYRHS
C            - SKYRHS
C
C
C     Notes:
C
C     CSTII: elemental diagonal jacobian matrix
C     CSTIT: global diagonal jacobian matrix
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON/SOLVERTA/KRENUT,KSOLVT,KSYMMT,NWIDTT,MITCGT,NBUFAT,NPRIRT,
     .                MITGMT,MKRYLT,IPGMRT
      COMMON/SOLVERTB/TOLCGT,TOLC1T,TOLGMT
C
      DIMENSION DISIT(NTOTV),   IFFIX(NTOTV),   REFOR(NTOTV),
     .          LNODS(NNODE,NELEM)
      DIMENSION CSTIF(NEVAB,*), DISIM(*),
     .          FOREL(*),
     .          CSTII(NEVAB),   CSTIT(NTOTV)
C
      DIMENSION MATNO(NELEM),
     .          PROEL(NPREL,NGRUP), PROPS(NPROP,NMATS),
     .          ELDAT(NDATA),       ELPRE(NPREV),
     .          ELVAR(NSTAT),       ELMAT(NMATX)
      DIMENSION DISTO(NTOTV,3),     COORD(NDIME,NPOIN)
      DIMENSION ADVEL(NTOTV*NDIME), TEMPI(NPOIN,2)
      DIMENSION PREAS(NPOIN),       TGAPS(NPOIN)
      DIMENSION DISPL(NTOTVM),      FPCHA(NFPCH,NPOIN),
     .          LACTI(NELEM)
C
C**** ASSEMBLE & FACTORISE MATRIX IF NECESSARY
C
      IF(KRESL.EQ.1) THEN
       CALL EXPASST(CSTIF,LNODS,
     .              NDOFN,NELEM,NEVAB,NNODE,NPOIN,
     .              NPREL,NGRUP,NPROP,NMATS,
     .              NDATA,NPREV,NSTAT,NMATX,NTOTV,NDIME,NTOTVM,NFPCH,
     .              MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .              ELVAR,ELMAT,WORK1,DISTO,COORD,DISIT,
     .              ADVEL,TEMPI,PREAS,TGAPS,DISPL,FPCHA,LACTI,
     .              CSTII,CSTIT)
C
      ENDIF
C
C**** IF NECESSARY MODIFY RHS FOR NON-ZERO BOUNDARY CONDITIONS
C
      IF(INRHS.EQ.1)
     . CALL SKYRHST(CSTIF,DISIT,REFOR,LNODS,NDOFN,NELEM,NEVAB,
     .              NNODE,NPOIN,DISIM,FOREL,    1,
     .              NPREL,NGRUP,NPROP,NMATS,
     .              NDATA,NPREV,NSTAT,NMATX,NTOTV,NDIME,NTOTVM,NFPCH,
     .              MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .              ELVAR,ELMAT,WORK1,DISTO,COORD,
     .              ADVEL,TEMPI,PREAS,TGAPS,DISPL,FPCHA,LACTI)
C
C**** SOLVE THE EQUATIONS 
C
      CALL EXPITET(DISIT,IFFIX,REFOR,
     .             NTOTV,
     .             CSTIT)
C
C**** COMPUTE REACTIONS
C
      IF(IPASS.EQ.2) THEN
       CALL SKYRHST(CSTIF,DISIT,REFOR,LNODS,NDOFN,NELEM,NEVAB,
     .              NNODE,NPOIN,DISIM,FOREL,    0,
     .              NPREL,NGRUP,NPROP,NMATS,
     .              NDATA,NPREV,NSTAT,NMATX,NTOTV,NDIME,NTOTVM,NFPCH,
     .              MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .              ELVAR,ELMAT,WORK1,DISTO,COORD,
     .              ADVEL,TEMPI,PREAS,TGAPS,DISPL,FPCHA,LACTI)
      ENDIF
C
      RETURN
      END
