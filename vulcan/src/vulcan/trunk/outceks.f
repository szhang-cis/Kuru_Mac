      SUBROUTINE OUTCEKS(IITERS,ISTEPS,KPRINS,NCHEKS,NOUTPS,NSTEPS,
     .                   TLIMTS,NCKGLO)
C***********************************************************************
C
C**** THIS ROUTINE CHECKS IF THE RESULTS SHOULD BE PRINTOUT
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
C       KPRI6=0  DO NOT WRITE NODAL MACROSCOPIC PHASE-CHANGE FUNCTION
C             1  WRITE NODAL MACROSCOPIC PHASE-CHANGE FUNCTION
C       KPRI7=0  DO NOT WRITE NODAL INTERNAL VARIABLES TO OUTPUT FILE
C             1  WRITE NODAL INTERNAL VARIABLES TO OUTPUT FILES
C       KPRI8=0  DO NOT NODAL WRITE POROSITY CRITERIA
C            =1  WRITE NODAL POROSITY CRITERIA
C
C       KPRIN=1  WRITE TO OUTPUT FILES TEMPERATURES AND HEAT FLUXES
C
C       KFEMV=0  DO NOT WRITE FOR POSTPROCESOR
C             1  WRITE FOR POSTPROCESOR
C
C       NOUTP( 1) RESTART INDICATOR
C       NOUTP( 2) POSTPROCESOR         INDICATOR FOR CONVERGED
C                                      ITERATIONS
C
C       NOUTP( 4) TEMPERATURES     OUTPUT  "      "      "         "
C       NOUTP( 5) GAUSS HEAT FLUXES   "    "      "      "         "
C       NOUTP( 6) 
C       NOUTP( 7) GAUSSIAN INT. VAR.  "    "      "      "         "
C       NOUTP( 8) NODAL HEAT FLUXES   "    "      "      "         "
C       NOUTP( 9) NODAL MACROSCOPIC PHASE-CHANGE FUNCTION
C       NOUTP(10) TOTAL               "    "      "      "         "
C       NOUTP(11) POSTPROCESOR
C       NOUTP(20) NODAL INT. VAR.     "    "      "      "         "
C       NOUTP(22) NODAL POROSITY CRITERIA
C
C                                      INDICATOR FOR NON CONVERGED
C                                      ITERATIONS
C
C       NOUTP(13) TEMPERATURES     OUTPUT  "      "   "      "        "
C       NOUTP(14) GAUSS HEAT FLUXES   "    "      "   "      "        "
C       NOUTP(15) 
C       NOUTP(16) GAUSSIAN INT. VAR.  "    "      "   "      "        "
C       NOUTP(17) NODAL HEAT FLUXES   "    "      "   "      "        "
C       NOUTP(18) NODAL MACROSCOPIC PHASE-CHANGE FUNCTION
C       NOUTP(19) TOTAL               "    "      "   "      "        "
C       NOUTP(21) NODAL INT. VAR.     "    "      "   "      "        "
C       NOUTP(23) NODAL POROSITY CRITERIA
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'inpo_oms.f'
C
      DIMENSION NOUTPS(50)
C
      KFEMVS=0
      KPRINS=0
      KPRI0S=0
      KPRI1S=0
      KPRI2S=0
      KPRI3S=0
      KPRI4S=0
      KPRI5S=0
      KPRI6S=0
      KPRI7S=0
      KPRI8S=0
      KPRI9S=0
      KPRI10S=0
C
      IF(NCHEKS.EQ.0.AND.NCKGLO.EQ.0) THEN        ! CONVERGED RESULTS
       IF(NOUTPS( 2).NE.0) THEN
        IF(MOD(ISTEPS,NOUTPS( 2)).EQ.0.OR.ISTEPS.EQ.NSTEPS) KFEMVS=1
       ENDIF
       IF(NOUTPS( 4).NE.0) THEN
        IF(MOD(ISTEPS,NOUTPS( 4)).EQ.0.OR.ISTEPS.EQ.NSTEPS) KPRI1S=1
       END IF
       IF(NOUTPS( 5).NE.0) THEN
        IF(MOD(ISTEPS,NOUTPS( 5)).EQ.0.OR.ISTEPS.EQ.NSTEPS) KPRI2S=1
       END IF
       IF(NOUTPS( 7).NE.0) THEN
        IF(MOD(ISTEPS,NOUTPS( 7)).EQ.0.OR.ISTEPS.EQ.NSTEPS) KPRI4S=1
       END IF
       IF(NOUTPS( 8).NE.0) THEN
        IF(MOD(ISTEPS,NOUTPS( 8)).EQ.0.OR.ISTEPS.EQ.NSTEPS) KPRI5S=1
       END IF
       IF(NOUTPS( 9).NE.0) THEN
        IF(MOD(ISTEPS,NOUTPS( 9)).EQ.0.OR.ISTEPS.EQ.NSTEPS) KPRI6S=1
       END IF
       IF(NOUTPS(20).NE.0) THEN
        IF(MOD(ISTEPS,NOUTPS(20)).EQ.0.OR.ISTEPS.EQ.NSTEPS) KPRI7S=1
       END IF
       IF(NOUTPS(22).NE.0) THEN
        IF(MOD(ISTEPS,NOUTPS(22)).EQ.0.OR.ISTEPS.EQ.NSTEPS) KPRI8S=1
       END IF
       IF(NOUTPS(10).NE.0) THEN
        IF(MOD(ISTEPS,NOUTPS(10)).EQ.0.OR.ISTEPS.EQ.NSTEPS) KPRINS=1
       END IF
      ELSE                                ! NON CONVERGED RESULTS
       IF(NOUTPS(11).NE.0) THEN
        IF(MOD(IITERS,NOUTPS(11)).EQ.0) KFEMVS=1
       ENDIF
       IF(NOUTPS(13).NE.0) THEN
        IF(MOD(IITERS,NOUTPS(13)).EQ.0) KPRI1S=1
       END IF
       IF(NOUTPS(14).NE.0) THEN
        IF(MOD(IITERS,NOUTPS(14)).EQ.0) KPRI2S=1
       END IF
       IF(NOUTPS(16).NE.0) THEN
        IF(MOD(IITERS,NOUTPS(16)).EQ.0) KPRI4S=1
       END IF
       IF(NOUTPS(17).NE.0) THEN
        IF(MOD(IITERS,NOUTPS(17)).EQ.0) KPRI5S=1
       END IF
       IF(NOUTPS(18).NE.0) THEN
        IF(MOD(IITERS,NOUTPS( 9)).EQ.0) KPRI6S=1
       END IF
       IF(NOUTPS(21).NE.0) THEN
        IF(MOD(IITERS,NOUTPS(21)).EQ.0) KPRI7S=1
       END IF
       IF(NOUTPS(23).NE.0) THEN
        IF(MOD(IITERS,NOUTPS(23)).EQ.0) KPRI8S=1
       END IF
       IF(NOUTPS(19).NE.0) THEN
        IF(MOD(IITERS,NOUTPS(19)).EQ.0) KPRINS=1
       END IF
      ENDIF
C
      IF(KPRINS.EQ.1) THEN
       KPRI0S=1
       KPRI1S=1
       KPRI2S=1
       KPRI4S=1
       KPRI5S=1
       KPRI6S=1
       KPRI7S=1
      ENDIF
      IF((KPRI1S+KPRI2S+KPRI4S+KPRI5S+KPRI6S+KPRI7S).GT.0) KPRI0S=1
C
      CALL CPUTIMS(CUCPUS)
      IF((CUCPUS.GE.TLIMTS).AND.(NCHEKS.EQ.0).AND.(NCKGLO.EQ.0)) THEN
       IF((NOUTPS(2)+NOUTPS(11)).GT.0) KFEMVS=1
      ENDIF
C
C**** SIMILAR TREATMENT FOR THE POSTPROCESS OPTIONS (not implem. yet)
C
cpost
cpost       only nodal values for postprocessing
cpost       see: strinpt, outpost, outdist, outgaut & outnodt
C
      RETURN
      END