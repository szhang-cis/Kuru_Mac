      SUBROUTINE ALGORSS
C***********************************************************************
C
C**** THIS ROUTINE SETS STIFFNESS REFORMING INDEXES: KSTIF , KRESL
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_oms.f'
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      KRESLS=0
      KSTIFS=0
C
C**** DECIDE ON REFORMING STIFFNESS
C
      IF(ITIMES*ISTEPS*IITERS.EQ.1)         KSTIFS=1 ! at the beginning
      IF(IRESTS.NE.0.AND.IITERS.EQ.1)       KSTIFS=1 ! after restart
      IF(NEWBOS*ISTEPS*IITERS.EQ.1)         KSTIFS=1 ! new boundary
C
      IF(NALGOS.EQ.2)                          KSTIFS=1 ! always
      IF((NALGOS.EQ.3).AND.(IITERS.EQ.KALGOS)) KSTIFS=1 ! at K iteration
      IF((NALGOS.EQ.4).AND.(MOD(ISTEPS,KALGOS).EQ.0).
     .                AND.(IITERS.EQ.1))       KSTIFS=1 ! every K steps
      IF((NALGOS.GT.2).AND.(KUNLDS.EQ.1))      KSTIFS=1 ! unloading
C
C**** DECIDE IF ASSEMBLE A NEW TOTAL STIFFNESS AND SOLVE
C
      KRESLS=KSTIFS                                      ! new stiffness
      IF(   (KDYNAS.EQ.1.OR.KPORES.EQ.2).AND.
     .      (IITERS.EQ.1.AND.KDTIMS.EQ.1)    ) KRESLS=1  ! new time inc.
C
C**** PHASE-CHANGE MATRIX CONTROL
C
      IF(KPROBS.EQ.4) THEN
       IF(IITERS.EQ.1) NSOL1S=1 !K is always updated for first iteration
       IF(IITERS.GT.NSOL1S) THEN
        KSTIFS=0
        KRESLS=1
       ENDIF
C
       IF(NMEMO7S.EQ.1) THEN
        IF(KSTIFS.NE.1.OR.KRESLS.NE.1)
     .   CALL RUNMENS('WARNING: KSTIF=KRESL SHOULD BE 1 FOR NMEMO7S=1')
        KSTIFS=1
        KRESLS=1
       ENDIF
      ENDIF          ! kprobs.eq.4
C
      RETURN
      END
