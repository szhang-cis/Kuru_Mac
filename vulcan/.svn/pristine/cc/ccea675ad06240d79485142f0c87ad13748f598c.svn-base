      SUBROUTINE ALGORST
C***********************************************************************
C
C**** THIS ROUTINE SETS STIFFNESS REFORMING INDEXES: KSTIF , KRESL
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      KRESLT=0
      KSTIFT=0
C
C**** DECIDE ON REFORMING STIFFNESS
C
      IF(ITIMET*ISTEPT*IITERT.EQ.1)         KSTIFT=1 ! at the beginning
      IF(IRESTT.NE.0.AND.IITERT.EQ.1)       KSTIFT=1 ! after restart
      IF(NEWBOT*ISTEPT*IITERT.EQ.1)         KSTIFT=1 ! new boundary
C
      IF(NALGOT.EQ.2)                          KSTIFT=1 ! always
      IF((NALGOT.EQ.3).AND.(IITERT.EQ.KALGOT)) KSTIFT=1 ! at K iteration
      IF((NALGOT.EQ.4).AND.(MOD(ISTEPT,KALGOT).EQ.0).
     .                AND.(IITERT.EQ.1))       KSTIFT=1 ! every K steps
      IF((NALGOT.GT.2).AND.(KUNLDT.EQ.1))      KSTIFT=1 ! unloading
C
C**** DECIDE IF ASSEMBLE A NEW TOTAL STIFFNESS AND SOLVE
C
      KRESLT=KSTIFT                                      ! new stiffness
      IF(   (KDYNAT.EQ.1.OR.KPORET.EQ.2).AND.
     .      (IITERT.EQ.1.AND.KDTIMT.EQ.1)    ) KRESLT=1  ! new time inc.
C
C**** PHASE-CHANGE MATRIX CONTROL
C
      IF(KPROBT.EQ.4) THEN
       IF(IITERT.EQ.1) NSOL1=1 ! K is always updated for first iteration
       IF(IITERT.GT.NSOL1) THEN
        KSTIFT=0
        KRESLT=1
       ENDIF
C
       IF(NMEMO7.EQ.1) THEN
        IF(KSTIFT.NE.1.OR.KRESLT.NE.1)
     .   CALL RUNMENT('WARNING: KSTIF=KRESL SHOULD BE 1 FOR NMEMO7=1')
        KSTIFT=1
        KRESLT=1
       ENDIF
      ENDIF          ! kprobt.eq.4
C
      RETURN
      END
