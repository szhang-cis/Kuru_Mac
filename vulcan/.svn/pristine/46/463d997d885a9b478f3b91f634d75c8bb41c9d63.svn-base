      SUBROUTINE TMARCH
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
      INCLUDE 'inte_om.f'
      INCLUDE 'prob_om.f'
C
      IITER=0
      NCHEK=1
      KDTIM=2
      IF(KINTE.EQ.1) TIMEX=1.0D0/TBETA
C
      IF(KOPTI.NE.0) THEN
       IF(KOPTI.EQ.3) THEN
        INTER=XTIME                                  ! LINEAR SCALE
        IF(MOD(ISTEP-1,INTER).EQ.0) THEN
         KDTIM=1
         DTIME=TTIME*9.0D0/XTIME
         TIMEX=1.0D0+1.0D0/DTIME-1.0D0/DLOG(1.0D0+DTIME)
        ENDIF
       ELSE
        KDTIM=1
        DTIME=TTIME*(10.0D0**(1.0D0/XTIME)-1.0D0)    ! LOG. SCALE
        IF(KOPTI.EQ.2)                               ! SANDHU
     .   TIMEX=1.0D0+1.0D0/DTIME-1.0D0/DLOG(1.0D0+DTIME)
       ENDIF
      ELSE
       IF(ISTEP.EQ.1) KDTIM=1
      ENDIF
C
C**** PRINTOUT TIME PARAMETERS
C
      TTIME=TTIME+DTIME
      IF(KINTE.EQ.1) TBETA=1.0D0/TIMEX
C
      WRITE(LURES,900) ISTEP,TTIME,DTIME
C
  900 FORMAT(1H1,//,132('*'),//,
     .       10X,'TIME STEP NO.',I5,//,
     .       15X,'CURRENT TIME        =',E15.5,/,
     .       15X,'TIME STEP           =',E15.5,/)
      RETURN
      END
