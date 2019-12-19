      SUBROUTINE CPUTIM(RTIME)
C***********************************************************************
C
C**** THIS ROUTINE FINDS OUT THE CPU TIME IN SEC.
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
      GO TO (10,20,30,20,50,20,70,80), NMACHIM
C
   10 CALL CPUTIM1(RTIME)
      RETURN
C
   20 CALL CPUTIM2(RTIME)
      RETURN
C
   30 CALL CPUTIM3(RTIME)
      RETURN
C
   50 CALL CPUTIM5(RTIME)
      RETURN
C
   70 CALL CPUTIM7(RTIME)
      RETURN
C
   80 CALL CPUTIM8(RTIME)
      RETURN
C
      END
