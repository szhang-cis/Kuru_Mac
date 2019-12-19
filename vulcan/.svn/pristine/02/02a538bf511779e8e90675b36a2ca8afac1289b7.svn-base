      SUBROUTINE ALGORS
C***********************************************************************
C
C**** THIS ROUTINE SETS STIFFNESS REFORMING INDEXES: KSTIF , KRESL
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      KRESL=0
      KSTIF=0
C
C**** DECIDE ON REFORMING STIFFNESS
C
      IF(ITIME*ISTEP*IITER.EQ.1)        KSTIF=1 ! at the beginning
      IF(IREST.NE.0.AND.IITER.EQ.1)     KSTIF=1 ! after restart
      IF(NEWBO*ISTEP*IITER.EQ.1)        KSTIF=1 ! new boundary
C
      IF(NALGO.EQ.2)                    KSTIF=1 ! always
      IF((NALGO.EQ.3).AND.(IITER.EQ.KALGO))  KSTIF=1 ! at K iteration
      IF((NALGO.EQ.4).AND.(MOD(ISTEP,KALGO).EQ.0).
     .                AND.(IITER.EQ.1))      KSTIF=1 ! every K steps
      IF((NALGO.GT.2).AND.(KUNLD.EQ.1))      KSTIF=1 ! unloading
C
C**** DECIDE IF ASSEMBLE A NEW TOTAL STIFFNESS AND SOLVE
C
      KRESL=KSTIF                                       ! new stiffness
      IF(   (KDYNA.EQ.1.OR.KPORE.EQ.2).AND.
     .      (IITER.EQ.1.AND.KDTIM.EQ.1)    )   KRESL=1  ! new time inc.
C
C**** OTHER CONTROLS
C
      IF(NMEMO7M.EQ.1) THEN
       IF(KSTIF.NE.1.OR.KRESL.NE.1)
     .  CALL RUNMEN('WARNING: KSTIF=KRESL SHOULD BE 1 FOR NMEMO7M=1')
       KSTIF=1
       KRESL=1
      ENDIF
C
      RETURN
      END
