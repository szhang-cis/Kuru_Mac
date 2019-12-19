      SUBROUTINE ITERAT(DISIT,DISPR,DISTO,ELDAT,ELPRE,ELVAR,ELMAT,
     .                  HEADS,IFFIX,PRESC,LNODS,MATNO,PROEL,PROPS,
     .                  REFOR,RLOAD,FICTO,TLOAD,TEMPN,DTEMP,INFRI,
     .                  COFRI,PWORK,PREAS,TGAPS,VNORM,LPNTN,COORD,
     .                  FPCHA,LACTI,NOPRF,PRESF,VANIS,WORK1,
     .                  NEQNS)
C***********************************************************************
C
C**** THIS ROUTINE PERFORMS THE ITERATIVE CORRECTION ON THE RESPONSE
C
C     Note: Dimensions of PREAS are inverted in order to properly
C           transfer it to stifmx.f, setmtx.f & output.f
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
      INCLUDE 'inpo_om.f'
C
      DIMENSION DISIT(NTOTV,2),     DISPR(NTOTV,NDISR),
     .          DISTO(NTOTV,NDISO),
     .          ELDAT(NDATA),       ELPRE(NPREV),       ELVAR(NSTAT),
     .          ELMAT(NMATX),       HEADS(NPOIN,4),
     .          IFFIX(NTOTV,2),     LNODS(NNODE,NELEM), MATNO(NELEM),
     .          PROEL(NPREL,NGRUP), PROPS(NPROP,NMATS), REFOR(NTOTV,2),
     .          RLOAD(NTOTV),       TLOAD(NTOTV,2),
     .          TEMPN(NPOIN,2),     DTEMP(NPOIN),LPNTN(*)
      DIMENSION INFRI(*),           COFRI(NSKEW,NDIME,*),
     .          PWORK(NPOIN,2),     PREAS(NPOIN,NPREA),
     .          TGAPS(NPOIN),       VNORM(NTOTV),
     .          COORD(NDIME,NPOIN), FPCHA(NFPCH,NPOIN)
      DIMENSION PRESC(NTOTV,2),     FICTO(NFUNC),
     .          LACTI(NELEM)
      DIMENSION NOPRF(NNODE,NELEM), PRESF(NNODE,NDIME,NELEM),
     .          VANIS(NANIV,NANIC,NELEM)
      REAL*8, ALLOCATABLE, INTENT(INOUT) :: WORK1(:)
      INCLUDE 'soladd.inc'
C
C**** UPDATE ITERATION COUNTER
C
      IITER=IITER+1
C
C**** ASSEMBLE GLOBAL COUPLED STIFFNESS MATRIX, IF NECESSARY
C
      CALL ALGORS      ! select solution flags KSTIF, KRESL
C
      IF(NMEMO7M.EQ.0) THEN
       NPRX5=NPRE1+NPRE2+NPRE3+NPRE4+1
       IF(KSTIF.EQ.1)  ! evaluate a new stiffness matrix
     .  CALL STIFMX(ELDAT,ELPRE,ELVAR,ELMAT,LNODS,MATNO,PROEL,
     .              PROPS,WORK1,TEMPN(1,1),DISPR(1,1),DTEMP,VNORM,COORD,
     .              DISTO(1,1),INFRI,COFRI,LACTI,PREAS(1,NPRX5))
      ENDIF
C
      IF(NMEMO6M.EQ.0) THEN
       IF(KRESL.EQ.1)  ! assemble global coupled stiffness matrix
     .  CALL ASELMT(MATNO,PROEL,PROPS,LNODS,INFRI,COFRI,
     .              ELMAT(IMATX(1)),ELMAT(IMATX(2)),ELMAT(IMATX(3)),
     .              ELMAT(IMATX(4)),ELMAT(IMATX(5)),ELMAT(IMATX(6)),
     .              WORK1)
      ENDIF
C
C**** INITIALISE THE ITER. DISPL. TAKING INTO ACCOUNT FIXITY CONDITIONS
C
      CALL INIDIS(DISIT,IFFIX,PRESC,FICTO,RLOAD)
C
C**** APPLY SKEW SYSTEM CONDITIONS ON REFOR (GLOBAL > LOCAL)
C
C     REFOR: RESIDUAL VECTOR
C
C     Note: DISIT (PRESCRIBED DISPLACEMENTS) does not need to be
C           transformed since the prescribed displacements at the nodes
C           with skew systems must be input in these local systems.
C           The skew (local) systems are only useful for prescribed
C           displacements since the loads are input either in the global
C           system (point, surface and volume loads) or in surface local
C           systems that are different from the skew systems (surface
C           loads).
C
      IF(NSKEW.GT.0) THEN
       CALL TRSKEW(NPOIN,NTOTV,NDOFN,REFOR,INFRI,COFRI,    1,
     .             NSKEW)
      ENDIF
C
C**** SOLVE EQUATIONS FOR ITERATIVE 'DISPLACEMENTS'
C
      CALL ESEPAS(DISIT,DISPR,IFFIX,LNODS,REFOR,RLOAD,TLOAD,WORK1,
     .            MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .            ELVAR,ELMAT,TEMPN,DTEMP,
     .            VNORM,DISTO,COORD,INFRI,COFRI,LACTI,
     .            NEQNS)
C
C**** APPLY SKEW SYSTEM CONDITIONS ON DISIT & REFOR (LOCAL > GLOBAL)
C
C     DISIT: ITERATIVE DISPLACEMENTS
C     REFOR: RESIDUAL VECTOR
C
      IF(NSKEW.GT.0) THEN
       CALL TRSKEW(NPOIN,NTOTV,NDOFN,DISIT,INFRI,COFRI,    2,
     .             NSKEW)
       CALL TRSKEW(NPOIN,NTOTV,NDOFN,REFOR,INFRI,COFRI,    2,
     .             NSKEW)
      ENDIF
C
C**** DEALS WITH CONTACT NON-COINCIDENT MESH
C     (only for non-linearized computation of normal vector n)
C
      IF(NOCOI.GT.0.AND.LICOI.EQ.1) THEN
       IF(LARGC.NE.0) THEN
        IPRIX=0
        IF(IITER.EQ.1) IPRIX=1
        CALL SETMTX(ELDAT,ELPRE,ELVAR,ELMAT,LNODS,MATNO,PROEL,PROPS,
     .              COORD,VNORM,DISTO,PREAS(1,NPRE1+1),
     .              PREAS(1,NPRE1+NPRE2+1),WORK1,IPRIX)
C
        CALL SOLADD(IFFIX,LNODS,LPNTN,PRESC,
     .              WORK1(IFIXY( 1)),WORK1(IFIXY( 2)),WORK1(IFIXY( 3)),
     .              WORK1(IFIXY( 4)),WORK1(IFIXY( 5)),WORK1(IFIXY( 6)),
     .              WORK1(IFIXY( 7)),
     .              NEQNS,WORK1)
       ENDIF
      ENDIF
C
C**** CALCULATE RESIDUAL FORCES
C
      CALL RESIDU(DISIT,DISPR,DISTO,ELDAT,ELPRE,ELVAR,ELMAT,
     .            HEADS,IFFIX,LNODS,MATNO,PROEL,PROPS,REFOR,
     .            TLOAD,TEMPN,DTEMP,INFRI,COFRI,PWORK,PREAS,
     .            TGAPS,VNORM,COORD,FPCHA,LACTI,NOPRF,PRESF,
     .            VANIS,WORK1)
C
C**** APPLY SKEW SYSTEM CONDITIONS ON REFOR & TLOAD (GLOBAL > LOCAL)
C
      IF(NSKEW.GT.0) THEN
       CALL TRSKEW(NPOIN,NTOTV,NDOFN,REFOR,INFRI,COFRI,    1,
     .             NSKEW)
       CALL TRSKEW(NPOIN,NTOTV,NDOFN,TLOAD,INFRI,COFRI,    1,
     .             NSKEW)
      ENDIF
C
C**** CHECK FOR CONVERGENCE
C
      CALL CONVER(DISIT,DISPR,DISTO,REFOR,TLOAD,IFFIX)
C
C**** APPLY SKEW SYSTEM CONDITIONS ON REFOR & TLOAD (LOCAL > GLOBAL)
C
      IF(NSKEW.GT.0) THEN
       CALL TRSKEW(NPOIN,NTOTV,NDOFN,REFOR,INFRI,COFRI,    2,
     .             NSKEW)
       CALL TRSKEW(NPOIN,NTOTV,NDOFN,TLOAD,INFRI,COFRI,    2,
     .             NSKEW)
      ENDIF
C
C**** OUTPUT RESULTS
C
      CALL OUTPUT(DISTO,ELDAT,ELPRE,ELVAR,ELMAT,HEADS,IFFIX,
     .            LNODS,MATNO,PROEL,PROPS,TLOAD,REFOR,COORD,
     .            PWORK(1,2),TGAPS,PREAS(1,1),DISTO,INFRI,LACTI,
     .            PREAS(1,NPRE1+1),PREAS(1,NPRE1+NPRE2+1),WORK1)
C
      RETURN
      END
