      SUBROUTINE COUPLE0S(LNODSS,MATNOS,PROELS,PROPSS,COORDS,HTLODS,
     .                    IFFIXS,PRESCS,RLOADS,RLOAHS,FICTOS,TFICTS,
     .                    DISITS,DISPRS,DISTOS,HEADSS,REFORS,TLOADS,
     .                    LPNTNS,ELDATS,ELPRES,ELVARS,ELMATS,DISPLS,
     .                    PWORKS,PREASS,TGAPSS,TEMPIS,ADVELS,FPCHAS,
     .                    LACTIS,
     .                    WORK1S,
     .                    LNODST,MATNOT,PROELT,PROPST,COORDT,HTLODT,
     .                    IFFIXT,PRESCT,RLOADT,RLOAHT,FICTOT,TFICTT,
     .                    DISITT,DISPRT,DISTOT,HEADST,REFORT,TLOADT,
     .                    LPNTNT,ELDATT,ELPRET,ELVART,ELMATT,DISPLT,
     .                    PWORKT,PREAST,TGAPST,TEMPIT,ADVELT,FPCHAT,
     .                    LACTIT,
     .                    WORK1T)
C***********************************************************************
C
C**** THIS ROUTINE SOLVES THE THERMO-MICROSTRUCTURAL COUPLED PROBLEM
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8, ALLOCATABLE, INTENT(INOUT) :: WORK1S(:),WORK1T(:)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
      INCLUDE 'nued_om.f'
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
C**** MICROSTRUCTURAL VARIABLES
C
      DIMENSION LNODSS(NNODES,NELEMS),  MATNOS(NELEMS),
     .          PROELS(NPRELS,NGRUPS),  PROPSS(NPROPS,NMATSS),
     .          COORDS(NDIMES,NPOINS),  HTLODS(NHLODS,NSUBFS,NFUNCS)
      DIMENSION IFFIXS(NTOTVS,2),       PRESCS(NTOTVS,2),
     .          RLOADS(NTOTVS),         RLOAHS(NTOTVS,NFUNCS),
     .          FICTOS(NFUNCS),         TFICTS(NFUNCS)
      DIMENSION DISITS(NTOTVS,2),
     .          DISPRS(NTOTVS,3),       DISTOS(NTOTVS,3),
     .          HEADSS(NPOINS,4),
     .          REFORS(NTOTVS,2),       TLOADS(NTOTVS,2),
     .          LPNTNS(NPOINS)
      DIMENSION ELDATS(NDATAS),         ELPRES(NPREVS),
     .          ELVARS(NSTATS),         ELMATS(NMATXS)
      DIMENSION DISPLS(NTOTVMS),        PWORKS(NPOINS,3),
     .          PREASS(NPREAS,NPOINS),  TGAPSS(NPOINS),
     .          TEMPIS(NPOINS,2),       ADVELS(NTOTVS*NDIMES),
     .          FPCHAS(NFPCH,NPOINS),   LACTIS(NELEMS)
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
