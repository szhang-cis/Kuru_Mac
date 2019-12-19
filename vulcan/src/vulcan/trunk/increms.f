      SUBROUTINE INCREMS(ELDATS,ELPRES,ELVARS,ELMATS,HEADSS,HTLODS,
     .                   IFFIXS,PRESCS,LNODSS,MATNOS,PROELS,PROPSS,
     .                   RLOADS,RLOAHS,FICTOS,TFICTS,TLOADS)
C***********************************************************************
C
C**** THIS ROUTINE INCREMENTS THE APPLIED LOADING
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      DIMENSION MATNOS(NELEMS),          LNODSS(NNODES,NELEMS),
     .          PROELS(NPRELS,NGRUPS),   PROPSS(NPROPS,NMATSS),
     .          ELDATS(NDATAS),          ELPRES(NPREVS),
     .          ELVARS(NSTATS),          ELMATS(NMATXS)
      DIMENSION HEADSS(NPOINS,*), 
     .          HTLODS(NHLODS,NSUBFS,*), IFFIXS(*),       
     .          RLOADS(*),               TLOADS(*)
      DIMENSION PRESCS(NTOTVS,2),        RLOAHS(NTOTVS,NFUNCS),
     .          FICTOS(NFUNCS),          TFICTS(NFUNCS)
C
C**** DETERMINE THE LOAD/DISPLACEMENT INCREMENT FACTOR
C
      CALL INCFACS(HTLODS,FICTOS,TFICTS)
C
C**** INCREMENT THE APPLIED LOAD
C
      CALL INCFORS(HEADSS,IFFIXS,RLOADS,RLOAHS,FICTOS,TLOADS)
C
      RETURN
      END
