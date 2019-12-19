      SUBROUTINE CPUTIMT(RTIMET)
C***********************************************************************
C
C**** THIS ROUTINE FINDS OUT THE CPU TIME IN SEC.
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
      GO TO (10,20,30,20,50,20,70,80), NMACHI
C
   10 CALL CPUTIMT1(RTIMET)
      RETURN
C
   20 CALL CPUTIMT2(RTIMET)
      RETURN
C
   30 CALL CPUTIMT3(RTIMET)
      RETURN
C
   50 CALL CPUTIMT5(RTIMET)
      RETURN
C
   70 CALL CPUTIMT7(RTIMET)
      RETURN
C
   80 CALL CPUTIMT8(RTIMET)
      RETURN
C
      END
