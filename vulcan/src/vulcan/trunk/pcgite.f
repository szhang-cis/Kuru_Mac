      SUBROUTINE PCGITE(DISIT,IFFIX,LNODS,REFOR,
     .                  NDOFN,NELEM,NEVAB,NNODE,NPOIN,NTOTV,
     .                  ALOAD,CSTIF,DELTA,GSTDI,DISIM,FOREL,
     .                  NPREL,NGRUP,NPROP,NMATS,
     .                  NDATA,NPREV,NSTAT,NMATX,NDIME,
     .                  NDISR,NDISO,
     .                  MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .                  ELVAR,ELMAT,WORK1,TEMPN,DISPR,DTEMP,
     .                  VNORM,DISTO,COORD,INFRI,COFRI,LACTI)
C***********************************************************************
C
C**** THIS ROUTINE SOLVES THE SYSTEM OF EQUATIONS USING  
C     A PRECONDITIONED CONJUGATED GRADIENT ITERATIVE SOLVER
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON/SOLVERA/KRENU,KSOLV,KSYMM,NWIDT,MITCG,NBUFA,NPRIR,
     .               MITGM,MKRYL,IPGMR,KPARAM, NTHMCM, NTHSOM
      COMMON/SOLVERB/TOLCG,TOLC1,TOLGM
C
      COMMON/LOGUN/LUDTS,LUSOL,LUFRO,LUFRH,LUDAT,LUPRI,LURES,
     .             LUSO2,LUFR2,LUPOS,LURST,LUBFG,LUPIP,LUPAN,LUINF,
     .             LUGEO,LUSET,LUMAT,LUINI,LULOA,LUFIX,LUIN1,
     .             LUTUN,LUCON,LUACT,LUFAN,
     .             LUCU1,LUCU2,LUCU3,LUCU4,LUCU5,LUCU6,LUCU7,
     .             LUCU8,LUCU9,LUC10
C
      DIMENSION DISIT(*), IFFIX(*), LNODS(NNODE,*), REFOR(*)
      DIMENSION ALOAD(*), CSTIF(*), DELTA(*), GSTDI(*), 
     .          DISIM(*), FOREL(*)
C
      DIMENSION MATNO(NELEM),
     .          PROEL(NPREL,NGRUP), PROPS(NPROP,NMATS),
     .          ELDAT(NDATA),       ELPRE(NPREV),
     .          ELVAR(NSTAT),       ELMAT(NMATX),
     .          WORK1(*)
      DIMENSION TEMPN(NPOIN,2),     DISPR(NTOTV,NDISR),
     .          DTEMP(NPOIN),       VNORM(NTOTV)
      DIMENSION DISTO(NTOTV,NDISO), COORD(NDIME,NPOIN)
      DIMENSION INFRI(*),           COFRI(*),
     .          LACTI(NELEM)
C
C**** INITIALISE SOME ARRAYS
C
      DO ITOTV=1,NTOTV
       IF(IFFIX(ITOTV).NE.0) REFOR(ITOTV)=0.0
       ALOAD(ITOTV)=REFOR(ITOTV)
       DELTA(ITOTV)=0.0
      ENDDO
C
      ANORM=VECDOT(REFOR,REFOR,NTOTV)
      IF(ANORM.LE.TOLC1) RETURN
C
C**** ITERATIVE SOLUTION
C
      BETAB=1.0D00
C
      DO 100 IITER=1,MITCG
C
C**** COMPUTE FIRST PART OF ITERATIVE SOLUTION
C
      CALL PCGPRE(ALOAD,GSTDI,FOREL,NWIDT,NDOFN,NPOIN)
C
C**** EVALUATE CG PARAMETER
C
      BETAT=VECDOT(ALOAD,REFOR,NTOTV)
      BETAK=BETAT/BETAB
      BETAB=BETAT
C
C**** EVALUATE ADVANCING DIRECTION
C
      DO ITOTV=1,NTOTV
       DELTA(ITOTV)=ALOAD(ITOTV)+BETAK*DELTA(ITOTV)
      ENDDO
C
C**** EVALUATE PRODUCT K*d
C
      CALL SKYRHS(CSTIF,DELTA,ALOAD,LNODS,NDOFN,NELEM,NEVAB,
     .            NNODE,NPOIN,DISIM,FOREL,    0,
     .            NPREL,NGRUP,NPROP,NMATS,
     .            NDATA,NPREV,NSTAT,NMATX,NTOTV,NDIME,
     .            NDISR,NDISO,
     .            MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .            ELVAR,ELMAT,WORK1,TEMPN,DISPR,DTEMP,
     .            VNORM,DISTO,COORD,INFRI,COFRI,LACTI)
C
      DO ITOTV=1,NTOTV
       IF(IFFIX(ITOTV).NE.0) ALOAD(ITOTV)=0.0
      ENDDO
C
C**** COMPUTE LINE SEARCH PARAMETER
C
      EETAT=VECDOT(DELTA,REFOR,NTOTV)
      EETAB=VECDOT(DELTA,ALOAD,NTOTV)
      EETAK=EETAT/EETAB
C
C**** UPDATE SOLUTION
C
      DO ITOTV=1,NTOTV
       DISIT(ITOTV)=DISIT(ITOTV)+EETAK*DELTA(ITOTV)
       REFOR(ITOTV)=REFOR(ITOTV)-EETAK*ALOAD(ITOTV)
       ALOAD(ITOTV)=REFOR(ITOTV)
      ENDDO
C
C**** CHECK FOR CONVERGENCE ( NORM OF RESIDUAL FORCES )
C
      BNORM=VECDOT(REFOR,REFOR,NTOTV)
      RATIO=DSQRT(BNORM/ANORM)
C
      IF(NPRIR.EQ.1) WRITE(LUPRI,900) IITER,RATIO
      IF(RATIO.LT.TOLCG) THEN
       RETURN
      ENDIF
C
  100 CONTINUE
C
C**** NOT CONVERGED
C
      CALL RUNEND('PCGITE: "PCG" SOLVER NOT CONVERGED ')
C
      RETURN
  900 FORMAT(5X,'PCG ITER. = ',I5,3X,'DPCG RATIO = ',E15.6)
      END   