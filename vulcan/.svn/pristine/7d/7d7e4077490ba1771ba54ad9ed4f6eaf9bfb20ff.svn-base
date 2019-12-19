       SUBROUTINE INTERVS(DISTOS,ELDATS,ELPRES,ELVARS,ELMATS,HTLODS,
     .                    IFFIXS,PRESCS,LNODSS,LPNTNS,MATNOS,PROELS,
     .                    PROPSS,RLOADS,RLOAHS,FICTOS,TFICTS,ADVELS,
     .                    HEADSS,TLOADS,COORDS,TEMPIS,FPCHAS,PREASS,
     .                    DISPLS,LACTIS,WORK1S)
        IMPLICIT REAL*8(A-H,O-Z)
        INCLUDE 'nuec_om.f'
        INCLUDE 'nued_om.f'
        INCLUDE 'prob_oms.f'
        INCLUDE 'inte_oms.f'
        INCLUDE 'auxl_oms.f'
        INCLUDE 'inpo_oms.f'


        DIMENSION MATNOS(NELEMS),          LNODSS(NNODES,NELEMS),
     .            PROELS(NPRELS,NGRUPS),   PROPSS(NPROPS,NMATSS),
     .            ELDATS(NDATAS),          ELPRES(NPREVS),
     .            ELVARS(NSTATS),          ELMATS(NMATXS)
        DIMENSION DISTOS(*),
     .            HTLODS(NHLODS,NSUBFS,*), IFFIXS(NTOTVS,*),
     .            LPNTNS(*),               RLOADS(*)
        DIMENSION PRESCS(NTOTVS,2),        RLOAHS(NTOTVS,NFUNCS),
     .            FICTOS(NFUNCS),          TFICTS(NFUNCS)
        DIMENSION ADVELS(NTOTVS*NDIMES)
        DIMENSION HEADSS(NPOINS,4),        TLOADS(NTOTVS,2),
     .            COORDS(NDIMES,NPOINS),   TEMPIS(NPOINS,2),
     .            FPCHAS(NFPCH,NPOINS),    PREASS(NPREAS,NPOINS),
     .            DISPLS(NTOTVMS),         LACTIS(NELEMS)

        REAL*8, ALLOCATABLE, INTENT(INOUT) :: WORK1S(:)
       END SUBROUTINE
