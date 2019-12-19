      SUBROUTINE ESEPAS(DISIT,DISPR,IFFIX,LNODS,REFOR,RLOAD,TLOAD,WORK1,
     .                  MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .                  ELVAR,ELMAT,TEMPN,DTEMP,
     .                  VNORM,DISTO,COORD,INFRI,COFRI,LACTI,
     .                  NEQNS)
C***********************************************************************
C
C**** THIS ROUTINE PERFOMS THE SOLUTION FOR UNCOUPLED PROBLEM
C     (CONSTRAINED METHODS)
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
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION DISIT(NTOTV,*),   DISPR(NTOTV,*),  IFFIX(NTOTV,*),
     .          LNODS(NNODE,*),   REFOR(NTOTV,*),  RLOAD(NTOTV),
     .          TLOAD(NTOTV,*)
      DIMENSION WORK1(*)
C
      DIMENSION MATNO(NELEM),
     .          PROEL(NPREL,NGRUP), PROPS(NPROP,NMATS),
     .          ELDAT(NDATA),       ELPRE(NPREV),
     .          ELVAR(NSTAT),       ELMAT(NMATX)
      DIMENSION TEMPN(NPOIN,2),
     .          DTEMP(NPOIN),       VNORM(NTOTV)
      DIMENSION DISTO(NTOTV,NDISO), COORD(NDIME,NPOIN)
      DIMENSION INFRI(*),           COFRI(*),
     .          LACTI(NELEM)
C
C**** SOLVE EQUATIONS  ( FIRST PASS )
C
      CALL ANTESO(REFOR(1,1),REFOR(1,2),DISIT(1,2),
     .            WORK1(IQUAS(1)),WORK1(IQUAS(2)))
C
      CALL SOLVER(DISIT(1,1),REFOR(1,1),IFFIX(1,1),
     .            IMPOS,KPORE,KRESL,
     .            LNODS,NDOFN,NELEM,NEVAB,NNODE,NPOIN,
     .            NTOTV,    1,IITER,WORK1,miter,
     .            NPREL,NGRUP,NPROP,NMATS,
     .            NDATA,NPREV,NSTAT,NMATX,NDIME,
     .            NDISR,NDISO,
     .            MATNO,PROEL,PROPS,ELDAT,ELPRE,
     .            ELVAR,ELMAT,TEMPN,DISPR,DTEMP,
     .            VNORM,DISTO,COORD,NPOIC,IAUGM,INFRI,COFRI,
     .            NACTI,LACTI,
     .            NEQNS)
C
      CALL POSTSO(DISIT(1,1),DISIT(1,2),REFOR(1,1),REFOR(1,2),
     .            WORK1(IQUAS(1)),WORK1(IQUAS(2)))
C
C***********************************************************************
C
      IF(KARCL.EQ.0.OR.ISTEP.LE.2) GOTO 100
C
C**** FOR ARC-LENGTH AND DISPLACEMENT CONTROL : (to be revised !!!)
C
        IF(KRESL.EQ.1.OR.(ISTEP.EQ.3.AND.IITER.EQ.1)) THEN
C
C - PREPARE FOR SECOND PASS
C
c            DO 10 ITOTV=1,NTOTV
c   10       IF(IFFIX(ITOTV).NE.0) DISPR(ITOTV,2)=RLOAD(ITOTV,1)
C
C - SOLVE EQUATIONS  ( SECOND PASS )
C
c            CALL SOLVER(DISPR(1,2),RLOAD(1,2),IFFIX,
c    .                   IMPOS,KPORE,KRESL,
c    .                   LNODS,NDOFN,NELEM,NEVAB,NNODE,NPOIN,
c    .                   NTOTV,    1,IITER,WORK1,miter,
c    .                   NPREL,NGRUP,NPROP,NMATS,
c    .                   NDATA,NPREV,NSTAT,NMATX,NDIME,
c    .                   NDISR,NDISO,
c    .                   MATNO,PROEL,PROPS,ELDAT,ELPRE,
c    .                   ELVAR,ELMAT,TEMPN,DISPR,DTEMP,
c    .                   VNORM,DISTO,COORD,NPOIC,IAUGM,INFRI,COFRI,
c    .                   NACTI,LACTI)
C
        ENDIF
C
C - UPDATE LOAD LEVEL AND ITERATIVE DISPLACEMENT VECTOR
C
c        CALL UPDATE(DISPR(1,2),DISPR(1,3),DISPR(1,1),DISIT,
c     .              NELEM,NEVAB,NTOTV,RLOAD,TLOAD)
C
 100  CONTINUE
C
      RETURN
      END
