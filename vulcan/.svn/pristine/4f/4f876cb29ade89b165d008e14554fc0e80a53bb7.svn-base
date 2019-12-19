      SUBROUTINE ASSIFIS
C***********************************************************************
C
C**** THIS ROUTINE ASSIGNS MICROSTRUCTURAL FILES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_oms.f'
C
C**** MICROSTRUCTURAL FILES ALREADY DEFINED IN prob_oms.f
C
      GO TO (10,10,30,10,50,10,10,80), NMACHIS
C
   10 CALL ASSIFIS1
      RETURN
C
   30 CALL ASSIFIS3
      RETURN
C
   50 CALL ASSIFIS5
      RETURN
C
   80 CALL ASSIFIS8
      RETURN
C
      END
