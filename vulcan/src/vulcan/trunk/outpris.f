      SUBROUTINE OUTPRIS(DISTOS,ELDATS,ELPRES,ELVARS,ELMATS,HEADSS,
     .                   LNODSS,MATNOS,PROELS,PROPSS,
     .                   SMSTSS,SMSTNS,SMSTPS,SMSPPS,
     .                   COORDS,FPCHAS,DISPLS,
     .                   WORK1S)
C***********************************************************************
C
C**** THIS ROUTINE OUTPUTS TEMPERATURES & HEAT FLUXES
C
C       KPRI0=0  DO NOT WRITE TO OUTPUT FILE
C             1  WRITE TO OUTPUT FILES
C       KPRI1=0  DO NOT WRITE TEMPERATURES TO OUTPUT FILE
C             1  WRITE TEMPERATURES TO OUTPUT FILES
C       KPRI2=0  DO NOT WRITE GAUSSIAN HEAT FLUXES TO OUTPUT FILE
C             1  WRITE GAUSSIAN HEAT FLUXES TO OUTPUT FILES
C       KPRI3=0  ---
C             1  ---
C       KPRI4=0  DO NOT WRITE GAUSS INTERNAL VARIABLES TO OUTPUT FILE
C             1  WRITE GAUSS INTERNAL VARIABLES TO OUTPUT FILES
C       KPRI5=0  DO NOT WRITE NODAL HEAT FLUXES
C             1  WRITE NODAL HEAT FLUXES
C       KPRI6=0  ---
C             1  ---
C       KPRI7=0  DO NOT WRITE NODAL INTERNAL VARIABLES TO OUTPUT FILE
C             1  WRITE NODAL INTERNAL VARIABLES TO OUTPUT FILES
C       KPRI8=0 DO NOT WRITE NODAL POROSITY CRITERIA
C            =1 WRITE NODAL POROSITY CRITERIA
C
C       KPRIN=1  WRITE TO OUTPUT FILES TEMPERATURES & HEAT FLUXES
C
C       KFEMV=0  DO NOT WRITE FOR POSTPROCESOR
C             1  WRITE FOR POSTPROCESOR
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
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
      INCLUDE 'inpo_oms.f'
C
      DIMENSION MATNOS(NELEMS),        LNODSS(NNODES,NELEMS),
     .          PROELS(NPRELS,NGRUPS), PROPSS(NPROPS,NMATSS),
     .          ELDATS(NDATAS),        ELPRES(NPREVS),
     .          ELVARS(NSTATS),        ELMATS(NMATXS)
      DIMENSION DISTOS(*),             HEADSS(NPOINS,*),
     .          SMSTSS(NSTR1S,NPOINS), SMSTNS(NSTR1S,NPOINS)
      DIMENSION SMSTPS(NNUINS,NPOINS), SMSPPS(NPOROS,NPOINS)
      DIMENSION STRSPS(3)
      DIMENSION COORDS(NDIMES,NPOINS), FPCHAS(NFPCH,NPOINS),
     .          DISPLS(NTOTVMS)
      DIMENSION WORK1S(*)
C
      IF(KFEMVS.NE.0)
     . OPEN(UNIT=LUPOSS,FILE=CJS,STATUS='OLD',FORM='UNFORMATTED',
     .      ACCESS='APPEND')
C
C**** OUTPUT TEMPERATURES
C
      CALL OUTDISS(DISTOS,HEADSS,DISPLS)
C
C**** OUTPUT GAUSSIAN VARIABLES
C
      IF((KPRI2S+KPRI4S).NE.0)
     . call runends('outgaus not implemented')
c    . CALL OUTGAUT(ELDATT,ELPRET,ELVART,ELMATT,
c    .              LNODST,MATNOT,PROELT,PROPST,COORDT,DISPLT,
c    .              WORK1T)
C
C**** OUTPUTS MACROSCOPICAL PHASE-CHANGES
C
c     CALL OUTMICS(FPCHAS)
C
C**** OUTPUT NODAL HEAT FLUXES & INTERNAL VARIABLES
C
      IF(((KPRI5S+KPRI7S+KPRI8S+KFEMVS).GT.0).AND.(KSGAUS.NE.0))
     . call runends('outnods not implemented')
c    . CALL OUTNODT(SMSTST,SMSTNT,SMSTPT,SMSPPT)
C
C**** CLOSE POSTPROCESS FILE
C
C     Warning: in CONVEX C-120, comment "close(lupost)" (february 1993)
C
      IF(KFEMVS.NE.0) CLOSE(LUPOSS)
C
      RETURN
      END
