      SUBROUTINE ASSIFIC
C***********************************************************************
C
C**** THIS ROUTINE ASSIGNS THERMAL FILES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
      INCLUDE 'addi_omt.f'
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'nuec_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
C
C**** COUPLING FILES ALREADY DEFINED IN nuec_om.f
C
      IF(NMACHI.NE.NMACHIM)
     . CALL RUNEND('ERROR: NMACHI NE NMACHIM')
C
      GO TO (10,10,30,10,50,10,10,80), NMACHI
C
   10 CALL ASSIFIC1
      RETURN
C
   30 CALL ASSIFIC3
      RETURN
C
   50 CALL ASSIFIC5
      RETURN
C
   80 CALL ASSIFIC8
      RETURN
C
      END
