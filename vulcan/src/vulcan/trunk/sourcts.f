      SUBROUTINE SOURCTS(DENSES,SOURCS,QINTES)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE SOURCE TERM
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
C**** CONSTANT PART OF THE SOURCE TERM
C
      QINTES=DENSES*SOURCS
C
C**** ADDS POSSIBLY VARIABLE PART OF THE SOURCE TERM
C
      RETURN
      END
