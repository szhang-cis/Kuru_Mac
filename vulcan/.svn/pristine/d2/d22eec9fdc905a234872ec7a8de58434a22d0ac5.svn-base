        SUBROUTINE ITERAT(DISIT,DISPR,DISTO,ELDAT,ELPRE,ELVAR,ELMAT,
     .                    HEADS,IFFIX,PRESC,LNODS,MATNO,PROEL,PROPS,
     .                    REFOR,RLOAD,FICTO,TLOAD,TEMPN,DTEMP,INFRI,
     .                    COFRI,PWORK,PREAS,TGAPS,VNORM,LPNTN,COORD,
     .                    FPCHA,LACTI,NOPRF,PRESF,VANIS,WORK1,!)!,
     .                    NEQNS)
          IMPLICIT REAL*8(A-H,O-Z)
          INCLUDE 'nuec_om.f'
          INCLUDE 'prob_om.f'
          INCLUDE 'inte_om.f'
          INCLUDE 'auxl_om.f'
          INCLUDE 'inpo_om.f'

          DIMENSION LNODS(NNODE,NELEM), MATNO(NELEM),
     .          PROEL(NPREL,NGRUP), PROPS(NPROP,NMATS),
     .          COORD(NDIME,NPOIN), HTLOD(NHLOD,NSUBF,NFUNC)
       DIMENSION IFFIX(NTOTV,2),     PRESC(NTOTV,2),
     .          RLOAD(NTOTV),       RLOAH(NTOTV,NFUNC),
     .          FICTO(NFUNC),       TFICT(NFUNC)
      DIMENSION DISIT(NTOTV,2),     DISPR(NTOTV,NDISR),
     .          DISTO(NTOTV,NDISO), HEADS(NPOIN,4),
     .          REFOR(NTOTV,2),     TLOAD(NTOTV,2)
      DIMENSION LPNTN(NPOIN),       ELDAT(NDATA),
     .          ELPRE(NPREV),       ELVAR(NSTAT),
     .          ELMAT(NMATX),       TEMPN(NPOIN,2),
     .          DTEMP(NPOIN)
      DIMENSION INFRI(NPOIN),       COFRI(NSKEW,NDIME,*),
     .          PWORK(NPOIN,2),     PREAS(NPOIN,NPREA),
     .          TGAPS(NPOIN),       VNORM(NTOTV),
     .          FPCHA(NFPCH,NPOIN), LACTI(NELEM)
      DIMENSION NOPRF(NNODE,NELEM), PRESF(NNODE,NDIME,NELEM),
     .          PREHF(NNODE,NDIME,NELEM,NFUNC),
     .          VANIS(2,3,NELEM)
        REAL*8, ALLOCATABLE, INTENT(INOUT) :: WORK1(:)
        END SUBROUTINE

