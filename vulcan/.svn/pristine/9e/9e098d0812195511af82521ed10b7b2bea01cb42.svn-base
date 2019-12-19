      SUBROUTINE TMSTEP(DISIT,DISPR,DISTO,ELDAT,ELPRE,ELVAR,ELMAT,
     .                  HEADS,HTLOD,IFFIX,PRESC,LNODS,MATNO,PROEL,
     .                  PROPS,REFOR,RLOAD,RLOAH,FICTO,TFICT,TLOAD,
     .                  TEMPN,DTEMP,INFRI,COFRI,PWORK,PREAS,TGAPS,
     .                  VNORM,COORD,FPCHA,LACTI,LPNTN,NOPRF,PRESF,
     .                  PREHF,VANIS,WORK1,
     .                  NEQNS)
C***********************************************************************
C
C**** THIS ROUTINE ADVANCES IN TIME, UPDATE TIME AND PRESCR. VARIABLES,
C     AND PERFORMS THE PREDICTOR STAGE
C
C     Note: Dimensions of PREAS are inverted in order to properly
C           transfer it to setmtx.f
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
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
C
      DIMENSION DISIT(NTOTV),             DISPR(NTOTV,NDISR),
     .          DISTO(NTOTV,NDISO),       ELDAT(NDATA),
     .          ELPRE(NPREV),             ELVAR(NSTAT),
     .          ELMAT(NMATX),             HEADS(NPOIN,4),
     .          HTLOD(NHLOD,NSUBF,NFUNC), IFFIX(NTOTV,2),
     .          LNODS(NNODE,NELEM),       MATNO(NELEM),
     .          PROEL(NPREL,NGRUP),       PROPS(NPROP,NMATS),
     .          REFOR(NTOTV),             RLOAD(NTOTV),
     .          TLOAD(NTOTV,2),
     .          TEMPN(NPOIN,2),           DTEMP(NPOIN)
      DIMENSION INFRI(*),                 COFRI(NSKEW,NDIME,*),
     .          PWORK(NPOIN,2),           PREAS(NPOIN,NPREA),
     .          TGAPS(NPOIN),             VNORM(NTOTV),
     .          COORD(NDIME,NPOIN),       FPCHA(NFPCH,NPOIN)
      DIMENSION PRESC(NTOTV,2),           RLOAH(NTOTV,NFUNC),
     .          FICTO(NFUNC),             TFICT(NFUNC),
     .          LACTI(NELEM),             LPNTN(*)
      DIMENSION NOPRF(NNODE,NELEM),       PRESF(NNODE,NDIME,NELEM),
     .          PREHF(NNODE,NDIME,NELEM,NFUNC),
     .          VANIS(NANIV,NANIC,NELEM)
      REAL*8, ALLOCATABLE, INTENT(INOUT) :: WORK1(:)
      INCLUDE 'soladd.inc'
C
      IF(KARCL.EQ.0.OR.ISTEP.LE.2) THEN
C
C**** UPDATE DTIME AND TTIME
C
       CALL TMARCH
C
C**** APPLY SKEW SYSTEM CONDITIONS ON RLOAH (GLOBAL > LOCAL)
C
       IF(NSKEW.GT.0) THEN
        DO IFUNC=1,NFUNC
         CALL TRSKEW(NPOIN,NTOTV,NDOFN,RLOAH(1,IFUNC),INFRI,COFRI,    1,
     .               NSKEW)
        ENDDO
       ENDIF
C
C**** INCREMENT LOAD FACTORS, APPLIED LOADS AND NON-TENSIONAL STRAINS
C
       CALL INCREM(ELDAT,      ELPRE,ELVAR,ELMAT,HEADS,HTLOD,
     .             IFFIX(1,1), PRESC,LNODS,MATNO,PROEL,PROPS,
     .             RLOAD,      RLOAH,PRESF,PREHF,FICTO,TFICT,TLOAD)
C
C**** APPLY SKEW SYSTEM CONDITIONS ON TLOAD & RLOAH (LOCAL > GLOBAL)
C
       IF(NSKEW.GT.0) THEN
        CALL TRSKEW(NPOIN,NTOTV,NDOFN,TLOAD,INFRI,COFRI,    2,
     .              NSKEW)
C
        DO IFUNC=1,NFUNC
         CALL TRSKEW(NPOIN,NTOTV,NDOFN,RLOAH(1,IFUNC),INFRI,COFRI,    2,
     .               NSKEW)
        ENDDO
       ENDIF
C
C**** PREDICTOR PHASE
C
       CALL PREDIC(DISIT,DISPR,DISTO)
C
C**** DEALS WITH CONTACT NON-COINCIDENT MESH
C
       IF(NOCOI.GT.0.AND.LICOI.EQ.0.AND.MOD(ISTEP,NSKIC).EQ.0) THEN
        IF(LARGC.NE.0) THEN
         CALL SETMTX(ELDAT,ELPRE,ELVAR,ELMAT,LNODS,MATNO,PROEL,PROPS,
     .               COORD,VNORM,DISTO,PREAS(1,NPRE1+1),
     .               PREAS(1,NPRE1+NPRE2+1),WORK1,    1)
C
         CALL SOLADD(IFFIX,LNODS,LPNTN,PRESC,
     .               WORK1(IFIXY( 1)),WORK1(IFIXY( 2)),WORK1(IFIXY( 3)),
     .               WORK1(IFIXY( 4)),WORK1(IFIXY( 5)),WORK1(IFIXY( 6)),
     .               WORK1(IFIXY( 7)),
     .               NEQNS,WORK1)
        ENDIF
       ENDIF
C
C**** CALCULATE RESIDUAL FORCES
C
       CALL RESIDU(DISIT,DISPR,DISTO,ELDAT,ELPRE,ELVAR,ELMAT,
     .             HEADS,IFFIX,LNODS,MATNO,PROEL,PROPS,REFOR,
     .             TLOAD,TEMPN,DTEMP,INFRI,COFRI,PWORK,PREAS,
     .             TGAPS,VNORM,COORD,FPCHA,LACTI,NOPRF,PRESF,
     .             VANIS,WORK1)
      ENDIF
C
      RETURN
      END
