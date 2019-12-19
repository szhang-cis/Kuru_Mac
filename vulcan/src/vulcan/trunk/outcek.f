      SUBROUTINE OUTCEK(IITER,ISTEP,KPRIN,NCHEK,NOUTP,NSTEP,TLIMT,
     .                  NCKGLO)
C***********************************************************************
C
C**** THIS ROUTINE CHECKS IF THE RESULTS SHOULD BE PRINTOUT
C
C      KPRI0= 0  DO NOT WRITE TO OUTPUT FILE
C             1  WRITE TO OUTPUT FILES
C      KPRI1= 0  DO NOT WRITE DISPLACEMENTS TO OUTPUT FILE
C             1  WRITE DISPLACEMENTS TO OUTPUT FILES
C      KPRI2= 0  DO NOT WRITE GAUSSIAN STRESSES TO OUTPUT FILE
C             1  WRITE GAUSSIAN STRESSES TO OUTPUT FILES
C      KPRI3= 0  DO NOT WRITE GAUSSIAN PRINCIPAL STRESSES TO OUTPUT FILE
C             1  WRITE GAUSSIAN PRINCIPAL STRESSES TO OUTPUT FILES
C      KPRI4= 0  DO NOT WRITE GAUSSIAN INTERNAL VARIABLES TO OUTPUT FILE
C             1  WRITE GAUSSIAN INTERNAL VARIABLES TO OUTPUT FILES
C      KPRI5= 0  DO NOT WRITE NODAL STRESSES
C             1  WRITE NODAL STRESSES
C      KPRI6= 0  DO NOT WRITE NODAL PRINCIPAL STRESSES
C             1  WRITE NODAL PRINCIPAL STRESSES
C      KPRI7= 0  DO NOT WRITE NODAL INTERNAL VARIABLES TO OUTPUT FILE
C             1  WRITE NODAL INTERNAL VARIABLES TO OUTPUT FILES
C      KPRI8= 0  DO NOT WRITE REACTIONS TO OUTPUT FILE
C             1  WRITE REACTIONS TO OUTPUT FILES
C      KPRI9= 0  DO NOT WRITE NODAL STRAINS
C             1  WRITE NODAL STRAINS
C      KPRI10=0  DO NOT WRITE NODAL PRINCIPAL STRAINS
C             1  WRITE NODAL PRINCIPAL STRAINS
C      KPRI11=0  DO NOT WRITE NORMAL GAP
C             1  WRITE NORMAL GAP
C      KPRI12=0  DO NOT WRITE HOMOGENIZED STRESS
C             1  WRITE HOMOGENIZED STRESS
C      KFEMV= 0  DO NOT WRITE FOR POSTPROCESOR
C             1  WRITE FOR POSTPROCESOR
C      KPRIN= 1  WRITE TO OUTPUT FILES DISPLACEMENTS AND STRESSES
C
C      NOUTP( 1) RESTART INDICATOR
C      NOUTP( 2) POSTPROCESOR
C
C                                     INDICATOR FOR CONVERGED
C                                     ITERATIONS
C
C      NOUTP( 4) DISPLACEMENTS    OUTPUT  "      "      "         "
C      NOUTP( 5) GAUSS STRESSES      "    "      "      "         "
C      NOUTP( 6) PRINCIPAL G. STRESS "    "      "      "         "
C      NOUTP( 7) GAUSSIAN INT. VAR.  "    "      "      "         "
C      NOUTP( 8) NODAL STRESSES      "    "      "      "         "
C      NOUTP( 9) PRINCIPAL N. STRESS."    "      "      "         "
C      NOUTP(10) TOTAL               "    "      "      "         "
C      NOUTP(11) POSTPROCESOR
C      NOUTP(20) NODAL INT. VAR.     "    "      "      "         "
C      NOUTP(22) REACTIONS           "    "      "      "         "
C      NOUTP(24) NODAL STRAINS       "    "      "      "         "
C      NOUTP(26) PRINCIPAL N. STRAIN "    "      "      "         "
C      NOUTP(28) NORMAL GAP          "    "      "      "         "
C      NOUTP(30) HOMOGENIZED STRESS  "    "      "      "         "
C
C                                     INDICATOR FOR NON CONVERGED
C                                     ITERATIONS
C
C      NOUTP(13) DISPLACEMENTS    OUTPUT  "      "      "         "
C      NOUTP(14) GAUSS STRESSES      "    "      "      "         "
C      NOUTP(15) PRINCIPAL G. STRESS."    "      "      "         "
C      NOUTP(16) GAUSSIAN INT. VAR.  "    "      "      "         "
C      NOUTP(17) NODAL STRESSES      "    "      "      "         "
C      NOUTP(18) PRINCIPAL N. STRESS."    "      "      "         "
C      NOUTP(19) TOTAL               "    "      "      "         "
C      NOUTP(21) NODAL INT. VAR.     "    "      "      "         "
C      NOUTP(23) REACTIONS           "    "      "      "         "
C      NOUTP(25) NODAL STRAINS       "    "      "      "         "
C      NOUTP(27) PRINCIPAL N. STRAIN "    "      "      "         "
C      NOUTP(29) NORMAL GAP          "    "      "      "         "
C      NOUTP(31) HOMOGENIZED STRESS  "    "      "      "         "
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'inpo_om.f'
C
      DIMENSION NOUTP(50)
C
      KFEMV=0
      KPRIN=0
      KPRI0=0
      KPRI1=0
      KPRI2=0
      KPRI3=0
      KPRI4=0
      KPRI5=0
      KPRI6=0
      KPRI7=0
      KPRI8=0
      KPRI9=0
      KPRI10=0
      KPRI11=0
      KPRI12=0
C
      IF(NCHEK.EQ.0.AND.NCKGLO.EQ.0) THEN           ! CONVERGED RESULTS
       IF(NOUTP( 2).NE.0) THEN
        IF(MOD(ISTEP,NOUTP( 2)).EQ.0.OR.ISTEP.EQ.NSTEP) KFEMV=1
       ENDIF
       IF(NOUTP( 4).NE.0) THEN
        IF(MOD(ISTEP,NOUTP( 4)).EQ.0.OR.ISTEP.EQ.NSTEP) KPRI1=1
       END IF
       IF(NOUTP( 5).NE.0) THEN
        IF(MOD(ISTEP,NOUTP( 5)).EQ.0.OR.ISTEP.EQ.NSTEP) KPRI2=1
       END IF
       IF(NOUTP( 6).NE.0) THEN
        IF(MOD(ISTEP,NOUTP( 6)).EQ.0.OR.ISTEP.EQ.NSTEP) KPRI3=1
       END IF
       IF(NOUTP( 7).NE.0) THEN
        IF(MOD(ISTEP,NOUTP( 7)).EQ.0.OR.ISTEP.EQ.NSTEP) KPRI4=1
       END IF
       IF(NOUTP( 8).NE.0) THEN
        IF(MOD(ISTEP,NOUTP( 8)).EQ.0.OR.ISTEP.EQ.NSTEP) KPRI5=1
       END IF
       IF(NOUTP( 9).NE.0) THEN
        IF(MOD(ISTEP,NOUTP( 9)).EQ.0.OR.ISTEP.EQ.NSTEP) KPRI6=1
       END IF
       IF(NOUTP(20).NE.0) THEN
        IF(MOD(ISTEP,NOUTP(20)).EQ.0.OR.ISTEP.EQ.NSTEP) KPRI7=1
       END IF
       IF(NOUTP(22).NE.0) THEN
        IF(MOD(ISTEP,NOUTP(22)).EQ.0.OR.ISTEP.EQ.NSTEP) KPRI8=1
       END IF
       IF(NOUTP(24).NE.0) THEN
        IF(MOD(ISTEP,NOUTP(24)).EQ.0.OR.ISTEP.EQ.NSTEP) KPRI9=1
       END IF
       IF(NOUTP(26).NE.0) THEN
        IF(MOD(ISTEP,NOUTP(26)).EQ.0.OR.ISTEP.EQ.NSTEP) KPRI10=1
       END IF
       IF(NOUTP(28).NE.0) THEN
        IF(MOD(ISTEP,NOUTP(28)).EQ.0.OR.ISTEP.EQ.NSTEP) KPRI11=1
       END IF
       IF(NOUTP(30).NE.0) THEN
        IF(MOD(ISTEP,NOUTP(30)).EQ.0.OR.ISTEP.EQ.NSTEP) KPRI12=1
       END IF
       IF(NOUTP(10).NE.0) THEN
        IF(MOD(ISTEP,NOUTP(10)).EQ.0.OR.ISTEP.EQ.NSTEP) KPRIN=1
       END IF
      ELSE                                ! NON CONVERGED RESULTS
       IF(NOUTP(11).NE.0) THEN
        IF(MOD(IITER,NOUTP(11)).EQ.0) KFEMV=1
       ENDIF
       IF(NOUTP(13).NE.0) THEN
        IF(MOD(IITER,NOUTP(13)).EQ.0) KPRI1=1
       END IF
       IF(NOUTP(14).NE.0) THEN
        IF(MOD(IITER,NOUTP(14)).EQ.0) KPRI2=1
       END IF
       IF(NOUTP(15).NE.0) THEN
        IF(MOD(IITER,NOUTP(15)).EQ.0) KPRI3=1
       END IF
       IF(NOUTP(16).NE.0) THEN
        IF(MOD(IITER,NOUTP(16)).EQ.0) KPRI4=1
       END IF
       IF(NOUTP(17).NE.0) THEN
        IF(MOD(IITER,NOUTP(17)).EQ.0) KPRI5=1
       END IF
       IF(NOUTP(18).NE.0) THEN
        IF(MOD(IITER,NOUTP(18)).EQ.0) KPRI6=1
       END IF
       IF(NOUTP(21).NE.0) THEN
        IF(MOD(IITER,NOUTP(21)).EQ.0) KPRI7=1
       END IF
       IF(NOUTP(23).NE.0) THEN
        IF(MOD(IITER,NOUTP(23)).EQ.0) KPRI8=1
       END IF
       IF(NOUTP(25).NE.0) THEN
        IF(MOD(ISTEP,NOUTP(25)).EQ.0) KPRI9=1
       END IF
       IF(NOUTP(27).NE.0) THEN
        IF(MOD(ISTEP,NOUTP(27)).EQ.0) KPRI10=1
       END IF
       IF(NOUTP(29).NE.0) THEN
        IF(MOD(ISTEP,NOUTP(29)).EQ.0.OR.ISTEP.EQ.NSTEP) KPRI11=1
       END IF
       IF(NOUTP(31).NE.0) THEN
        IF(MOD(ISTEP,NOUTP(31)).EQ.0.OR.ISTEP.EQ.NSTEP) KPRI12=1
       END IF
       IF(NOUTP(19).NE.0) THEN
        IF(MOD(IITER,NOUTP(19)).EQ.0) KPRIN=1
       END IF
      ENDIF
C
      IF(KPRIN.EQ.1) THEN
       KPRI0=1
       KPRI1=1
       KPRI2=1
       KPRI3=1
       KPRI4=1
       KPRI5=1
       KPRI6=1
       KPRI7=1
       KPRI8=1
       KPRI9=1
       KPRI10=1
       KPRI11=1
       KPRI12=1
      ENDIF
      IF((KPRI1+KPRI2+KPRI3+KPRI4+KPRI5+KPRI6+KPRI7+KPRI8+KPRI9+
     . KPRI10+KPRI11+KPRI12).GT.0) KPRI0=1
C
      CALL CPUTIM(CUCPU)
      IF((CUCPU.GE.TLIMT).AND.(NCHEK.EQ.0).AND.(NCKGLO.EQ.0)) THEN
       IF((NOUTP(2)+NOUTP(11)).GT.0) KFEMV=1
      ENDIF
C
C**** SIMILAR TREATMENT FOR THE POSTPROCESS OPTIONS (not implem. yet)
C
cpost
cpost       only nodal values for postprocessing
cpost       see: strinp, outpos, outdis, outgau & outnod
C
      RETURN
      END
