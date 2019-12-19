      SUBROUTINE SKYGRA(DISIT,REFOR,LNODS,
     .                  NDOFN,NELEM,NEVAB,NNODE,NPOIN,NTOTV,
     .                  IDSK2,INRHS,IPASS,KPASS,KPORE,KRESL,NEQNS,NLAST,
     .                  IFFIX,
     .                  GSTDI,GSTLO,GSTUP,CSTIF,ELOAD,LNUEQ,LPONT,
     .                  DISIM,FOREL,ALOAD,DELTA,DISIC,LOCAL,
     .                  iiter,miter,
     .                  NPREL,NGRUP,NPROP,NMATS,
     .                  NDATA,NPREV,NSTAT,NMATX,NDIME,
     .                  NDISR,NDISO,
     .                  MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .                  ELVAR,ELMAT,WORK1,TEMPN,DISPR,DTEMP,
     .                  VNORM,DISTO,COORD,NPOIC,IAUGM,INFRI,COFRI,
     .                  NACTI,LACTI)
C***********************************************************************
C
C**** THIS ROUTINE SOLVES THE SET OF LINEAR EQUATIONS USING THE
C         SEmi Direct And iteraTIVE SOLVER ( "SEDATIVE" ) 
C
C     SKYDES
C     SKYGRA
C            - SKYASS 
C            - SKYTRI
C                      - SKYCEK
C                      - SKYDIA
C            - SKYTAP
C            - SKYRHS 
C         !  - SKYBAK              ! Direct Skyline Solver
C         !  - SKYITE  - SKYBAK    ! Iterative PCG Solver
C                      - SKYPRO
C            - SKYRHS
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
      DIMENSION DISIT(NTOTV),   REFOR(NTOTV),  LNODS(NNODE,NELEM)
      DIMENSION LNUEQ(NTOTV),   LPONT(NTOTV)      
      DIMENSION GSTDI(*),       GSTLO(*),      GSTUP(*), 
     .          CSTIF(NEVAB,*), ELOAD(*),      DISIM(*),
     .          FOREL(*),       ALOAD(*),      DELTA(*),
     .          DISIC(*),       LOCAL(*)
      DIMENSION IFFIX(*)
C
      DIMENSION MATNO(NELEM),
     .          PROEL(NPREL,NGRUP), PROPS(NPROP,NMATS),
     .          ELDAT(NDATA),       ELPRE(NPREV),
     .          ELVAR(NSTAT),       ELMAT(NMATX)
      DIMENSION TEMPN(NPOIN,2),     DISPR(NTOTV,NDISR),
     .          DTEMP(NPOIN),       VNORM(NTOTV)
      DIMENSION DISTO(NTOTV,NDISO), COORD(NDIME,NPOIN)
      DIMENSION INFRI(*),           COFRI(*),
     .          LACTI(NELEM)
C
C**** IDENTIFY THE DESTINATION OF THE EQUATIONS INTO PROFILE
C     NEQNS,NLAST,LNUEQ,LPONT
C
      IF((KPASS.EQ.0).OR.((KPORE.EQ.2).AND.(NEQNS.GT.0))) THEN
       REWIND IDSK2
       READ(IDSK2) NEQNS,NLAST,LNUEQ,LPONT
       KPASS=1
      ENDIF
C
C**** ASSEMBLE & FACTORISE MATRIX IF NECESSARY
C
      IF(KRESL.EQ.1) THEN
       CALL SKYASS(CSTIF,GSTDI,GSTLO,GSTUP,LNUEQ,LPONT,KSYMM,
     .             LNODS,NDOFN,NELEM,NEQNS,NEVAB,NLAST,NNODE,
     .             NWIDT,NTOTV,NPOIN,
     .             NPREL,NGRUP,NPROP,NMATS,
     .             NDATA,NPREV,NSTAT,NMATX,NDIME,
     .             NDISR,NDISO,
     .             MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .             ELVAR,ELMAT,WORK1,TEMPN,DISPR,DTEMP,
     .             VNORM,DISTO,COORD,INFRI,COFRI,LACTI)
C
       CALL SKYTRI(GSTDI,GSTLO,GSTUP,KSYMM,LPONT,KPORE,NEQNS,LURES,
     .             NPOIC,IAUGM)
      ENDIF
C
C**** STORE OR READ BACK THE FACTORISED GLOBAL STIFFNESS MATRICES
C
      IF(KPORE.EQ.2) THEN
       IF(NEQNS.GT.0)
     .  CALL SKYTAP(KRESL,KSYMM,GSTDI,GSTLO,GSTUP,IDSK2,NLAST,
     .              NEQNS)
      ENDIF
C
C**** IF NECESSARY MODIFY RHS FOR NON-ZERO BOUNDARY CONDITIONS
C
      IF(INRHS.EQ.1)
     . CALL SKYRHS(CSTIF,DISIT,REFOR,LNODS,NDOFN,NELEM,NEVAB,
     .             NNODE,NPOIN,DISIM,FOREL,    1,
     .             NPREL,NGRUP,NPROP,NMATS,
     .             NDATA,NPREV,NSTAT,NMATX,NTOTV,NDIME,
     .             NDISR,NDISO,
     .             MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .             ELVAR,ELMAT,WORK1,TEMPN,DISPR,DTEMP,
     .             VNORM,DISTO,COORD,INFRI,COFRI,LACTI)
C
C**** COMPRESS RHS
C
      DO ITOTV=1,NTOTV
       IEQNS=LNUEQ(ITOTV)
       IF(IEQNS.GT.0) ELOAD(IEQNS)=REFOR(ITOTV)
      ENDDO
C
C**** SOLVE THE EQUATIONS 
C
      IF(MITCG.EQ.0) THEN                                ! Skyline
       CALL SKYBAK(GSTDI,GSTLO,GSTUP,LPONT,NEQNS,ELOAD,LURES)
      ELSE                                               ! PCG
       CALL SKYITE(ELOAD,LNODS,NDOFN,NELEM,NEQNS,NEVAB,NNODE,
     .             NPOIN,NTOTV,CSTIF,LNUEQ,LPONT,DISIM,FOREL,
     .             ALOAD,DELTA,DISIC,GSTDI,GSTLO,GSTUP,LOCAL,
     .             NPREL,NGRUP,NPROP,NMATS,
     .             NDATA,NPREV,NSTAT,NMATX,NDIME,
     .             NDISR,NDISO,
     .             MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .             ELVAR,ELMAT,WORK1,TEMPN,DISPR,DTEMP,
     .             VNORM,DISTO,COORD,INFRI,COFRI,LACTI)
      ENDIF
C
C**** REALLOCATE SOLUTION ARRAY INTO TOTAL NTOTV
C
      DO ITOTV=1,NTOTV
       IEQNS=LNUEQ(ITOTV)
       IF(IEQNS.GT.0) DISIT(ITOTV)=ELOAD(IEQNS)
      ENDDO
C
C**** COMPUTE REACTIONS (only necessary for coupled prob. with kpore=2)
C
      IF(IPASS.EQ.2) THEN
       CALL SKYRHS(CSTIF,DISIT,REFOR,LNODS,NDOFN,NELEM,NEVAB,
     .             NNODE,NPOIN,DISIM,FOREL,    0,
     .             NPREL,NGRUP,NPROP,NMATS,
     .             NDATA,NPREV,NSTAT,NMATX,NTOTV,NDIME,
     .             NDISR,NDISO,
     .             MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .             ELVAR,ELMAT,WORK1,TEMPN,DISPR,DTEMP,
     .             VNORM,DISTO,COORD,INFRI,COFRI,LACTI)
      ENDIF
C
      RETURN
      END
