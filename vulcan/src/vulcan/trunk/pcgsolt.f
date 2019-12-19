      SUBROUTINE PCGSOLT(DISIT,IFFIX,REFOR,LNODS,
     .                   NDOFN,NELEM,NEVAB,NNODE,NPOIN,NTOTV,
     .                   INRHS,IPASS,KRESL,
     .                   GSTDI,CSTIF,DISIM,FOREL,ALOAD,DELTA,
     .                   NPREL,NGRUP,NPROP,NMATS,
     .                   NDATA,NPREV,NSTAT,NMATX,NDIME,NTOTVM,NFPCH,
     .                   MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .                   ELVAR,ELMAT,WORK1,DISTO,COORD,
     .                   ADVEL,TEMPI,PREAS,TGAPS,DISPL,FPCHA,
     .                   NACTI,LACTI)
C***********************************************************************
C
C**** THIS ROUTINE SOLVES THE SET OF LINEAR EQUATIONS USING THE
C         PRECONDITIONED CONJUGATE GRADIENT ITERATIVE SOLVER 
C
C     PCGSOL
C            - PCGASS 
C            - PCGTRI
C            - PCGITE  - PCGPRE - PCGMUL
C                      - SKYRHS
C            - SKYRHS
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
      DIMENSION GSTDI(*),       CSTIF(NEVAB,*), DISIM(*),
     .          FOREL(*),       ALOAD(*),       DELTA(*)
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
       CALL PCGASST(CSTIF,GSTDI,LNODS,
     .              NDOFN,NELEM,NEVAB,NNODE,NPOIN,NWIDTT,
     .              NPREL,NGRUP,NPROP,NMATS,
     .              NDATA,NPREV,NSTAT,NMATX,NTOTV,NDIME,NTOTVM,NFPCH,
     .              MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .              ELVAR,ELMAT,WORK1,DISTO,COORD,DISIT,
     .              ADVEL,TEMPI,PREAS,TGAPS,DISPL,FPCHA,LACTI)
C
       CALL PCGTRI(GSTDI,IFFIX,NDOFN,NPOIN,NWIDTT,TOLC1T)
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
      CALL PCGITET(DISIT,IFFIX,LNODS,REFOR,
     .             NDOFN,NELEM,NEVAB,NNODE,NPOIN,NTOTV,
     .             ALOAD,CSTIF,DELTA,GSTDI,DISIM,FOREL,
     .             NPREL,NGRUP,NPROP,NMATS,
     .             NDATA,NPREV,NSTAT,NMATX,NDIME,NTOTVM,NFPCH,
     .             MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .             ELVAR,ELMAT,WORK1,DISTO,COORD,
     .             ADVEL,TEMPI,PREAS,TGAPS,DISPL,FPCHA,LACTI)
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
