      SUBROUTINE TMARCHS
C***********************************************************************
C
C**** THIS ROUTINE CHOOSES THE TIME MARCHING SCHEME
C
C     KDTIM=1 EITHER THETA-PARAMETER OR TIME-INCREMENT HAS CHANGED
C     KDTIM=2 BOTH THETA-PARAMETER AND TIME INCREMENT ARE UNCHANGED
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'inte_oms.f'
      INCLUDE 'prob_oms.f'
C
      IITERS=0
      NCHEKS=1
      KDTIMS=2
      IF(KINTES.EQ.1)TIMEXS=1./TBETAS
C
      IF(KOPTIS.NE.0)THEN
        IF(KOPTIS.EQ.3)THEN
          INTERS=XTIMES                            ! LINEAR SCALE
          IF(MOD(ISTEPS-1,INTERS).EQ.0)THEN
            KDTIMS=1
            DTIMES=TTIMES*9.0/XTIMES
            TIMEXS=1.0+1.0/DTIMES-1.0/DLOG(1+DTIMES)
          ENDIF
        ELSE
          KDTIMS=1
          DTIMES=TTIMES*(10.0**(1.0/XTIMES)-1.0)    ! LOG. SCALE
          IF(KOPTIS.EQ.2)                         ! SANDHU
     .     TIMEXS=1.0+1.0/DTIMES-1.0/DLOG(1+DTIMES)
        ENDIF
      ELSE
        IF(ISTEPS.EQ.1) KDTIMS=1
      ENDIF
C
C***PRINTOUT TIME PARAMETERS
C
      TTIMES=TTIMES+DTIMES
      IF(KINTES.EQ.1) TBETAS=1./TIMEXS
C
      WRITE(LURESS,900) ISTEPS,TTIMES,DTIMES
C
  900 FORMAT(1H1,//,132('*'),//,
     .       10X,'TIME STEP NO.',I5,//,
     .       15X,'CURRENT TIME        =',E15.5,/,
     .       15X,'TIME STEP           =',E15.5,/)
      RETURN
      END
