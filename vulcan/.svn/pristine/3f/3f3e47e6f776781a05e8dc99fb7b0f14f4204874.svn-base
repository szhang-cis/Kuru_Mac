      SUBROUTINE OUTPUTT(DISTOT,ELDATT,ELPRET,ELVART,ELMATT,HEADST,
     .                   IFFIXT,LNODST,MATNOT,PROELT,PROPST,TLOADT,
     .                   COORDT,FPCHAT,PREAST,DISPLT,LACTIT,WORK1T)
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
C       KPRI4=0  DO NOT WRITE GAUSSIAN INTERNAL VARIABLES TO OUTPUT FILE
C             1  WRITE GAUSSIAN INTERNAL VARIABLES TO OUTPUT FILES
C       KPRI5=0  DO NOT WRITE NODAL HEAT FLUXES
C             1  WRITE NODAL HEAT FLUXES
C       KPRI6=0  ---
C             1  ---
C       KPRI7=0 DO NOT WRITE NODAL INTERNAL VARIABLES TO OUTPUT FILE
C             1 WRITE NODAL INTERNAL VARIABLES TO OUTPUT FILES
C       KPRI8=0 DO NOT WRITE NODAL POROSITY CRITERIA
C            =1 WRITE NODAL POROSITY CRITERIA
C
C       KPRIN=1 WRITE TO OUTPUT FILES TEMPERATURES & HEAT FLUXES
C
C       KFEMV=0 DO NOT WRITE FOR POSTPROCESOR
C             1 WRITE FOR POSTPROCESOR
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
      INCLUDE 'nued_om.f'
      INCLUDE 'nuef_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'
C
      DIMENSION MATNOT(NELEMT),        LNODST(NNODET,NELEMT),
     .          PROELT(NPRELT,NGRUPT), PROPST(NPROPT,NMATST),
     .          ELDATT(NDATAT),        ELPRET(NPREVT),
     .          ELVART(NSTATT),        ELMATT(NMATXT),
     .          WORK1T(*)
      DIMENSION DISTOT(NTOTVT,3),      HEADST(NPOINT,*), 
     .          IFFIXT(*),             TLOADT(*)
      DIMENSION COORDT(NDIMET,NPOINT), FPCHAT(NFPCH,NPOINT),
     .          PREAST(NPOROT,NPOINT), DISPLT(NTOTVM),
     .          LACTIT(NELEMT)
C
      CALL CPUTIMT(TIME1T)
C
C**** IF NECESSARY PRINT OUT X-Y CURVES
C
      IF(NCHEKT.EQ.0.AND.NCURVT.NE.0.AND.NCKGLO.EQ.0)
     . CALL PLOOUTT(DISTOT,ELVART,TLOADT)
C
C**** CHECK IF THE RESULTS SHOULD BE PRINTOUT
C
      CALL OUTCEKT(IITERT,ISTEPT,KPRINT,NCHEKT,NOUTPT,NSTEPT,TLIMTT,
     .             NCKGLO)
C
C**** PRINTOUT CONVERGENCE INFORMATION
C
      IF(IPRCOT.GT.0) THEN       ! no print of initial conditions
       WRITE(LUREST,897) ITIMET,ISTEPT,IITERT,NCHEKT,NCKGLO
C
C**** PRINTOUT GENERAL TITLE
C
       IF(KPRI0T.EQ.0.AND.KFEMVT.EQ.1.AND.NCHEKT.EQ.0.AND.NCKGLO.EQ.0) 
     .                          WRITE(LUREST,885) ITIMET,ISTEPT,IITERT
       IF(KPRI0T.EQ.1.AND.NCHEKT.EQ.0.AND.NCKGLO.EQ.0)
     .                          WRITE(LUREST,890) ITIMET,ISTEPT,IITERT
       IF(KPRI0T.EQ.1.AND.NCHEKT.NE.0)
     .                          WRITE(LUREST,895) ITIMET,ISTEPT,IITERT
C
C**** TO DEAL WITH COUPLED PROBLEMS
C
       IF(KPRI0T.EQ.1.AND.NCHEKT.EQ.0.AND.NCKGLO.EQ.1)
     .                          WRITE(LUREST,896) ITIMET,ISTEPT,IITERT
      ENDIF
C
C**** PERFORMS SMOOTHING OVER GAUSSIAN HEAT FLUXES, INTERNAL VARIABLES &
C     POROSITY, IF NECESSARY
C
C     Notes:
C
C     - Smoothing is always necessary when porosity is used.
C     - Porosity is only computed and output for converged results.
C     However, note that when not-converged results are required for the
C     postprocess file (determined at the input data), zero porosity
C     variables are written in the .pos file in order to not break
C     the output variables order for the output interfaces: intfgvt,
C     intergplos, etc. (see outsmot.f & outnodt.f).
C
C     - Smoothing is always necessary for microstructural thermally
C     coupled flow problems (warning: the option implemented is only
C     valid for the locally converged coupling strategy; the option
C     IF(ITERMEF.GT.0.AND.IMICR.EQ.1) ITHCFP=1 is valid for the locally
C     non-converged strategy.
C
      ICALPO=0
      IF(NPOROT.NE.0.AND.NCHEKT.EQ.0.AND.NCKGLO.EQ.0) ICALPO=1
C
c     ITHCFP=0
      ITHCFP=1
      IF(ITERMEF.GT.0.AND.IMICR.EQ.1.AND.NCHEKT.EQ.0.AND.NCKGLO.EQ.0)
     .                                                          ITHCFP=1
C
      IF(KSGAUT.NE.0) THEN
       IF((KPRI5T+KPRI7T+KPRI8T).GT.0.OR.(KFEMVT.EQ.1).OR.
     .                              (ICALPO.EQ.1).OR.(ITHCFP.EQ.1)) THEN
        CALL OUTSMOT(ELDATT,ELPRET,ELVART,ELMATT,LNODST,MATNOT,PROELT,
     .               PROPST,COORDT,PREAST,FPCHAT,DISTOT(1,1),
     .               DISTOT(1,2),DISPLT,LACTIT,
     .               WORK1T(IGSMOT(1)),WORK1T(IGSMOT(2)),
     .               WORK1T(IGSMOT(3)),
     .               WORK1T(IGSMOT(8)),WORK1T(IGSMOT(22)),  !smstp,smspp
     .               WORK1T)
       ENDIF
      ENDIF
C
C**** PRINTOUT RESULTS TO RESULTS FILES
C
      IF(KPRI0T.EQ.1.OR.KFEMVT.EQ.1)
     . CALL OUTPRIT(DISTOT,ELDATT,ELPRET,ELVART,ELMATT,HEADST,
     .              LNODST,MATNOT,PROELT,PROPST,
     .              WORK1T(IGSMOT(1)),WORK1T(IGSMOT( 2)),
     .              WORK1T(IGSMOT(8)),WORK1T(IGSMOT(22)),
     .              COORDT,FPCHAT,DISPLT,WORK1T)
C
      CALL CPUTIMT(TIME2T)
      CPUOUT=CPUOUT+(TIME2T-TIME1T)
C
      RETURN
  885 FORMAT(1H1,///,
     .10X,'* * *  CONVERGED RESULTS STORED ONLY ON FILES * * *',//,
     .        5X,'ITIME =',I5,   5X,'ISTEP =',I5,   5X,'IITER =',I5)
  890 FORMAT(1H1,///,
     .10X,'* * * CONVERGED RESULTS * * *',//,
     .        5X,'ITIME =',I5,   5X,'ISTEP =',I5,   5X,'IITER =',I5)
  895 FORMAT(//,
     .10X,'* * * NOT-CONVERGED RESULTS * * *',//,
     .        5X,'ITIME =',I5,   5X,'ISTEP =',I5,   5X,'IITER =',I5)
C
  896 FORMAT(//,
     .10X,'* * * CONVERGED RESULTS - NOT-CONVERGED GLOBAL 
     .           RESULTS * * *',//,
     .        5X,'ITIME =',I5,   5X,'ISTEP =',I5,   5X,'IITER =',I5)
C
  897 FORMAT(//,
     .'* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *',
     .' * * * * * * *',//,
     .10X,'* * * CONVERGENCE INFORMATION * * *',//,
     .        3X,'ITIME =',I5,   3X,'ISTEP =',I5,   3X,'IITER =',I5,
     .        3X,'NCHEK =',I5,   3X,'NCKGLO =',I5,//,
     .'* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *',
     .' * * * * * * *')
C
      END
