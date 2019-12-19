      SUBROUTINE TMARCHT
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
      INCLUDE 'inte_omt.f'
      INCLUDE 'prob_omt.f'
C
      IITERT=0
      NCHEKT=1
      KDTIMT=2
      IF(KINTET.EQ.1)TIMEXT=1./TBETAT
C
      IF(KOPTIT.NE.0)THEN
        IF(KOPTIT.EQ.3)THEN
          INTERT=XTIMET                            ! LINEAR SCALE
          IF(MOD(ISTEPT-1,INTERT).EQ.0)THEN
            KDTIMT=1
            DTIMET=TTIMET*9.0/XTIMET
            TIMEXT=1.0+1.0/DTIMET-1.0/DLOG(1+DTIMET)
          ENDIF
        ELSE
          KDTIMT=1
          DTIMET=TTIMET*(10.0**(1.0/XTIMET)-1.0)    ! LOG. SCALE
          IF(KOPTIT.EQ.2)                         ! SANDHU
     .     TIMEXT=1.0+1.0/DTIMET-1.0/DLOG(1+DTIMET)
        ENDIF
      ELSE
        IF(ISTEPT.EQ.1) KDTIMT=1
      ENDIF
C
C***PRINTOUT TIME PARAMETERS
C
      TTIMET=TTIMET+DTIMET
      IF(KINTET.EQ.1) TBETAT=1./TIMEXT
C
      WRITE(LUREST,900) ISTEPT,TTIMET,DTIMET
C
  900 FORMAT(1H1,//,132('*'),//,
     .       10X,'TIME STEP NO.',I5,//,
     .       15X,'CURRENT TIME        =',E15.5,/,
     .       15X,'TIME STEP           =',E15.5,/)
      RETURN
      END
