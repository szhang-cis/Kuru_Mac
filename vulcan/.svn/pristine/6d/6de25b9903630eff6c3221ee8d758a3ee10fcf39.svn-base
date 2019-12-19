      SUBROUTINE COUINI(LNODS,MATNO,PROEL,PROPS,COORD,HTLOD,
     .                  IFFIX,PRESC,RLOAD,RLOAH,FICTO,TFICT,
     .                  DISIT,DISPR,DISTO,HEADS,REFOR,TLOAD,
     .                  LPNTN,ELDAT,ELPRE,ELVAR,ELMAT,TEMPN,
     .                  DTEMP,INFRI,COFRI,PWORK,PREAS,TGAPS,
     .                  VNORM,FPCHA,LACTI,NOPRF,PRESF,PREHF,
     .                  VANIS,WORK1,NEQNS,
     .                  LNODST,MATNOT,PROELT,PROPST,COORDT,HTLODT,
     .                  IFFIXT,PRESCT,RLOADT,RLOAHT,FICTOT,TFICTT,
     .                  DISITT,DISPRT,DISTOT,HEADST,REFORT,TLOADT,
     .                  LPNTNT,ELDATT,ELPRET,ELVART,ELMATT,DISPLT,
     .                  PWORKT,PREAST,TGAPST,TEMPIT,ADVELT,FPCHAT,
     .                  LACTIT,
     .                  WORK1T)
C***********************************************************************
C
C**** THIS ROUTINE ESTABLISH THE INITIAL CONDITIONS FOR THE MECHANICAL
C     PROBLEM
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8, ALLOCATABLE, INTENT(INOUT) :: WORK1(:), WORK1T(:)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
      INCLUDE 'nued_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
      INCLUDE 'inpo_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'

C
C**** MECHANICAL VARIABLES
C
      DIMENSION LNODS(NNODE,NELEM), MATNO(NELEM),       
     .          PROEL(NPREL,NGRUP), PROPS(NPROP,NMATS),
     .          COORD(NDIME,NPOIN), HTLOD(NHLOD,NSUBF,NFUNC)
      DIMENSION IFFIX(NTOTV,2),     PRESC(NTOTV,2),
     .          RLOAD(NTOTV),       RLOAH(NTOTV,NFUNC),
     .          FICTO(NFUNC),       TFICT(NFUNC)
      DIMENSION DISIT(NTOTV,2),     DISPR(NTOTV,NDISR),
     .          DISTO(NTOTV,NDISO), HEADS(NPOIN,4),
     .          REFOR(NTOTV,2),     TLOAD(NTOTV,2)
      DIMENSION LPNTN(NPOIN)
      DIMENSION ELDAT(NDATA),       ELPRE(NPREV),
     .          ELVAR(NSTAT),       ELMAT(NMATX)
      DIMENSION TEMPN(NPOIN,2),     DTEMP(NPOIN),
     .          INFRI(NPOIN),       COFRI(NSKEW,NDIME,*),
     .          PWORK(NPOIN,2),     PREAS(NPREA,NPOIN),
     .          TGAPS(NPOIN),       VNORM(NTOTV),
     .          FPCHA(NFPCH,NPOIN), LACTI(NELEM)
      DIMENSION NOPRF(NNODE,NELEM), PRESF(NNODE,NDIME,NELEM),
     .          PREHF(NNODE,NDIME,NELEM,NFUNC),
     .          VANIS(NANIV,NANIC,NELEM)
C
C**** THERMAL VARIABLES
C
      DIMENSION LNODST(NNODET,NELEMT), MATNOT(NELEMT),
     .          PROELT(NPRELT,NGRUPT), PROPST(NPROPT,NMATST),
     .          COORDT(NDIMET,NPOINT), HTLODT(NHLODT,NSUBFT,NFUNCT)
      DIMENSION IFFIXT(NTOTVT,2),      PRESCT(NTOTVT,2),
     .          RLOADT(NTOTVT),        RLOAHT(NTOTVT,NFUNCT),
     .          FICTOT(NFUNCT),        TFICTT(NFUNCT)
      DIMENSION DISITT(NTOTVT,2),  
     .          DISPRT(NTOTVT,3),      DISTOT(NTOTVT,3),
     .          HEADST(NPOINT,4),
     .          REFORT(NTOTVT,2),      TLOADT(NTOTVT,2),
     .          LPNTNT(NPOINT)
      DIMENSION ELDATT(NDATAT),        ELPRET(NPREVT),
     .          ELVART(NSTATT),        ELMATT(NMATXT)
      DIMENSION DISPLT(NTOTVM),        PWORKT(NPOINT,3),
     .          PREAST(NPREAT,NPOINT), TGAPST(NPOINT),
     .          TEMPIT(NPOINT,2),      ADVELT(NTOTVT*NDIMET),
     .          FPCHAT(NFPCH,NPOINT),  LACTIT(NELEMT)
      END SUBROUTINE
