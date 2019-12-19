      SUBROUTINE SOLVER(DISIT,REFOR,IFFIX,IMPOS,KPORE,KRESL,LNODS,
     .                  NDOFN,NELEM,NEVAB,NNODE,NPOIN,NTOTV,IPASS,
     .                  IITER,WORK1,miter,
     .                  NPREL,NGRUP,NPROP,NMATS,
     .                  NDATA,NPREV,NSTAT,NMATX,NDIME,
     .                  NDISR,NDISO,
     .                  MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .                  ELVAR,ELMAT,TEMPN,DISPR,DTEMP,
     .                  VNORM,DISTO,COORD,NPOIC,IAUGM,INFRI,COFRI,
     .                  NACTI,LACTI,
     .                  NEQNS)
C***********************************************************************
C
C**** THIS ROUTINE SOLVES THE SYSTEM OF EQUATIONS SELECTING THE
C     THE APPROPRIATE SOLVER
C
C     KRESL=1  ASSEMBLE & FACTORISE THE MATRIX (see algors.f)
C     IPASS=1  (with no arc-length or displ. control; see esepas.f)
C     KPASS=0  (see below)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'auxl_om.f'
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
      DIMENSION DISIT(*), IFFIX(*), LNODS(NNODE,*),  REFOR(*),
     .          WORK1(*)
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
      CALL CPUTIM(TIME1)
C
      INRHS=0
      IF(IITER.EQ.1) INRHS=IMPOS
ctm   IF(NEWBO.EQ.1) KPASS=0
      KPASS=0
ctm
      IDSK2=LUSOL
      IF(IPASS.EQ.2) THEN
       INRHS=1
       IDSK2=LUSO2
      ENDIF
C
C**** CALL THE PROFILE SOLVER
C
      IF(KSOLV.EQ.0)
     .  CALL SKYGRA(DISIT,REFOR,LNODS,
     .              NDOFN,NELEM,NEVAB,NNODE,NPOIN,NTOTV,
     .              IDSK2,INRHS,IPASS,KPASS,KPORE,KRESL,NEQNS,NLAST,
     .              IFFIX,
     .              WORK1(ISOLV( 1)),WORK1(ISOLV( 2)),WORK1(ISOLV( 3)),
     .              WORK1(ISOLV( 4)),WORK1(ISOLV( 5)),WORK1(ISOLV( 6)),
     .              WORK1(ISOLV( 7)),WORK1(ISOLV( 8)),WORK1(ISOLV( 9)),
     .              WORK1(ISOLV(10)),WORK1(ISOLV(11)),WORK1(ISOLV(12)),
     .              WORK1(ISOLV(13)),
     .              iiter,miter,
     .              NPREL,NGRUP,NPROP,NMATS,
     .              NDATA,NPREV,NSTAT,NMATX,NDIME,
     .              NDISR,NDISO,
     .              MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .              ELVAR,ELMAT,WORK1,TEMPN,DISPR,DTEMP,
     .              VNORM,DISTO,COORD,NPOIC,IAUGM,INFRI,COFRI,
     .              NACTI,LACTI)

C
C**** CALL THE FRONTAL SOLVER
C
      IF(KSOLV.EQ.1) THEN
       IDSK3=LUFRO
       IF(IPASS.EQ.2) IDSK3=LUFR2
       CALL FRONTS(DISIT,IDSK2,IDSK3,IFFIX,IPASS,LNODS,KRESL,
     .             NDOFN,NELEM,NEVAB,NNODE,NPOIN,KPORE,REFOR,
     .             WORK1(ISOLV(1)), WORK1(ISOLV(2)), WORK1(ISOLV(3)),
     .             WORK1(ISOLV(4)), WORK1(ISOLV(5)), WORK1(ISOLV(6)),
     .             WORK1(ISOLV(7)), WORK1(ISOLV(8)), WORK1(ISOLV(9)),
     .             WORK1(ISOLV(10)),WORK1(ISOLV(11)),WORK1(ISOLV(12)),
     .             KSYMM,
     .             NPREL,NGRUP,NPROP,NMATS,
     .             NDATA,NPREV,NSTAT,NMATX,NTOTV,NDIME,
     .             NDISR,NDISO,
     .             MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .             ELVAR,ELMAT,WORK1,TEMPN,DISPR,DTEMP,
     .             VNORM,DISTO,COORD,NPOIC,IAUGM,INFRI,COFRI,
     .             NACTI,LACTI)
      ENDIF
C
C**** CALL THE PCG SOLVER
C
      IF(KSOLV.EQ.2)
     . CALL PCGSOL(DISIT,IFFIX,REFOR,LNODS,
     .             NDOFN,NELEM,NEVAB,NNODE,NPOIN,NTOTV,
     .             INRHS,IPASS,KRESL,
     .             WORK1(ISOLV(1)),WORK1(ISOLV(2)),WORK1(ISOLV(3)),
     .             WORK1(ISOLV(4)),WORK1(ISOLV(5)),WORK1(ISOLV(6)),
     .             NPREL,NGRUP,NPROP,NMATS,
     .             NDATA,NPREV,NSTAT,NMATX,NDIME,
     .             NDISR,NDISO,
     .             MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .             ELVAR,ELMAT,WORK1,TEMPN,DISPR,DTEMP,
     .             VNORM,DISTO,COORD,INFRI,COFRI,
     .             NACTI,LACTI)
C
C**** CALL THE GMRES SOLVER
C
      IF(KSOLV.EQ.3)
     . CALL GMRESS(WORK1(ISOLV( 4)),NEQNS,WORK1(ISOLV(15)),
     .             WORK1(ISOLV(16)),WORK1(ISOLV(17)),
     .             DISIT,REFOR,
     .             WORK1(ISOLV( 1)),WORK1(ISOLV( 6)),
     .             LNODS,NDOFN,NELEM,NEVAB,NNODE,IPASS,NTOTV,
     .             NPOIN,
     .             NLAST,IDSK2,WORK1(ISOLV( 7)),       ! added variables
     .             MITGM,MKRYL,IPGMR,TOLGM,
     .             LURES)
     
     
      IF(KSOLV.EQ.4)
     . CALL PARDISM(DISIT,REFOR,LNODS, !---------
     .              NDOFN,NELEM,NEVAB,NNODE,NPOIN,NTOTV,!---------
     .              IDSK2,INRHS,IPASS,KPASS,KPORE,KRESL,NEQNS,NLAST,!-----
     .              IFFIX,!---------
     .              WORK1(ISOLV( 1)),WORK1(ISOLV( 2)),WORK1(ISOLV( 3)),
     .              WORK1(ISOLV( 4)),WORK1(ISOLV( 5)),WORK1(ISOLV( 6)),
     .              WORK1(ISOLV( 7)),WORK1(ISOLV( 8)),WORK1(ISOLV( 9)),
     .              WORK1(ISOLV(10)),WORK1(ISOLV(11)),WORK1(ISOLV(12)),
     .              WORK1(ISOLV(13)),!---------
     .              iiter,miter,!---------
     .              NPREL,NGRUP,NPROP,NMATS,!---------
     .              NDATA,NPREV,NSTAT,NMATX,NDIME,!---------
     .              NDISR,NDISO,!---------
     .              MATNO,PROEL,PROPS,ELDAT,ELPRE,!---------
     .              ELVAR,ELMAT,WORK1,TEMPN,DISPR,DTEMP,!---------
     .              VNORM,DISTO,COORD,NPOIC,IAUGM,INFRI,COFRI,!---------
     .              NACTI,LACTI)!---------
C
      CALL CPUTIM(TIME2)
      CPUSO=CPUSO+(TIME2-TIME1)
C
      RETURN
      END
