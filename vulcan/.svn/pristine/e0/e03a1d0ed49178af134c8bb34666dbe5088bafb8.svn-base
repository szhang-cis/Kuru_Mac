      SUBROUTINE ASSIFI
C***********************************************************************
C
C**** THIS ROUTINE ASSIGNS MECHANICAL FILES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** MECHANICAL FILES ALREADY DEFINED IN prob_om.f
C
      GO TO (10,10,30,10,50,10,10,80), NMACHIM
C
   10 CALL ASSIFI1
      RETURN
C
   30 CALL ASSIFI3
      RETURN
C
   50 CALL ASSIFI5
      RETURN
C
   80 CALL ASSIFI8
      RETURN
C
      END
