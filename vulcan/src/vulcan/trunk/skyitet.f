      SUBROUTINE SKYITET(ELOAD,LNODS,NDOFN,NELEM,NEQNS,NEVAB,NNODE,
     .                   NPOIN,NTOTV,CSTIF,LNUEQ,LPONT,DISIM,FOREL,
     .                   ALOAD,DELTA,DISIC,GSTDI,GSTLO,GSTUP,LOCAL,
     .                   NPREL,NGRUP,NPROP,NMATS,
     .                   NDATA,NPREV,NSTAT,NMATX,NDIME,NTOTVM,NFPCH,
     .                   MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .                   ELVAR,ELMAT,WORK1,DISTO,COORD,DISIT,
     .                   ADVEL,TEMPI,PREAS,TGAPS,DISPL,FPCHA,LACTI)
C***********************************************************************
C
C**** THIS ROUTINE SOLVES THE SYSTEM OF EQUATIONS USING  
C     A PRECONDITIONED CONJUGATED GRADIENT ITERATIVE SOLVER
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON/SOLVERTA/KRENUT,KSOLVT,KSYMMT,NWIDTT,MITCGT,NBUFAT,NPRIRT,
     .                MITGMT,MKRYLT,IPGMRT
      COMMON/SOLVERTB/TOLCGT,TOLC1T,TOLGMT
C
      COMMON/LOGUNT/LUDTST,LUSOLT,LUFROT,LUFRHT,LUDATT,LUPRIT,LUREST,
     .              LUSO2T,LUFR2T,LUPOST,LURSTT,LUBFGT,LUPIPT,LUPANT,
     .              LUGEOT,LUSETT,LUMATT,LUINIT,LULOAT,LUFIXT,LUADVT,
     .              LUACTT,LUFANT,LUSTRT,
     .              LUCU1T,LUCU2T,LUCU3T,LUCU4T,LUCU5T,LUCU6T,LUCU7T,
     .              LUCU8T,LUCU9T,LUC10T
C
      DIMENSION ALOAD(*), CSTIF(*), DELTA(*), DISIC(*),
     .          DISIM(*), ELOAD(*), FOREL(*), LNODS(NNODE,*),
     .          LNUEQ(*), LPONT(*), GSTDI(*), GSTLO(*), GSTUP(*),
     .          LOCAL(*)
C
      DIMENSION MATNO(NELEM),
     .          PROEL(NPREL,NGRUP), PROPS(NPROP,NMATS),
     .          ELDAT(NDATA),       ELPRE(NPREV),
     .          ELVAR(NSTAT),       ELMAT(NMATX)
      DIMENSION DISTO(NTOTV,3),     COORD(NDIME,NPOIN),
     .          DISIT(NTOTV)
      DIMENSION ADVEL(NTOTV*NDIME), TEMPI(NPOIN,2)
      DIMENSION PREAS(NPOIN),       TGAPS(NPOIN)
      DIMENSION DISPL(NTOTVM),      FPCHA(NFPCH,NPOIN),
     .          LACTI(NELEM)
C
      ANORM=VECDOT(ELOAD,ELOAD,NEQNS)
      IF(ANORM.LE.TOLC1T) RETURN
C
      BETAB=1.0D00
C
C**** INITIALISE SOME ARRAYS
C
      DO 10 IEQNS=1,NEQNS
      ALOAD(IEQNS)=ELOAD(IEQNS)
      DELTA(IEQNS)=0.0
      DISIC(IEQNS)=0.0
   10 CONTINUE
C
C**** ITERATIVE SOLUTION
C
      DO 100 IITER=1,MITCGT
C
C**** COMPUTE FIRST PART OF ITERATIVE SOLUTION
C
      CALL SKYBAK(GSTDI,GSTLO,GSTUP,LPONT,NEQNS,ALOAD,LUREST)
C
C**** EVALUATE CG PARAMETER
C
      BETAT=VECDOT(ALOAD,ELOAD,NEQNS)
      BETAK=BETAT/BETAB
C
C**** EVALUATE ADVANCING DIRECTION
C
      DO 30 IEQNS=1,NEQNS
   30 DELTA(IEQNS)=ALOAD(IEQNS)+BETAK*DELTA(IEQNS)
C
C**** EVALUATE PRODUCT K*d
C
      CALL SKYPROT(CSTIF,DELTA,LNODS,LNUEQ,NDOFN,NELEM,NEQNS,
     .             NEVAB,NNODE,NPOIN,ALOAD,DISIM,FOREL,LOCAL,
     .             NPREL,NGRUP,NPROP,NMATS,
     .             NDATA,NPREV,NSTAT,NMATX,NDIME,NTOTV,NTOTVM,NFPCH,
     .             MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .             ELVAR,ELMAT,WORK1,DISTO,COORD,
     .             ADVEL,TEMPI,PREAS,TGAPS,DISPL,FPCHA,LACTI)
C
C**** COMPUTE LINE SEARCH PARAMETER
C
      EETAT=VECDOT(DELTA,ELOAD,NEQNS)
      EETAB=VECDOT(DELTA,ALOAD,NEQNS)
      EETAK=EETAT/EETAB
      BETAB=BETAT
C
C**** UPDATE SOLUTION
C
      DO 40 IEQNS=1,NEQNS
       DISIC(IEQNS)=DISIC(IEQNS)+EETAK*DELTA(IEQNS)
       ELOAD(IEQNS)=ELOAD(IEQNS)-EETAK*ALOAD(IEQNS)
       ALOAD(IEQNS)=ELOAD(IEQNS)
   40 CONTINUE
C
C**** CHECK FOR CONVERGENCE ( NORM OF RESIDUAL FORCES )
C
      BNORM=VECDOT(ELOAD,ELOAD,NEQNS)
      RATIO=DSQRT(BNORM/ANORM)
C
      IF(RATIO.LT.TOLCGT) THEN
       WRITE(LUPRIT,900) IITER,RATIO
       DO 50 IEQNS=1,NEQNS
       ELOAD(IEQNS)=DISIC(IEQNS)
   50  CONTINUE
       RETURN
      ENDIF
C
  100 CONTINUE
C
C**** NOT CONVERGED
C
      CALL RUNENDT(' "SEDATIVE" SOLVER NOT CONVERGED   ')
C
      RETURN
  900 FORMAT(5X,'PCG ITER. = ',I5,3X,'PCG RATIO = ',E15.6)
      END   
