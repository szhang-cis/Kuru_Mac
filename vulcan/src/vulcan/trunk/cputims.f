      SUBROUTINE CPUTIMS(RTIMES)
C***********************************************************************
C
C**** THIS ROUTINE FINDS OUT THE CPU TIME IN SEC.
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_oms.f'
C

      rtimes=0.0


c     GO TO (10,20,30,20,50,20,70,80), NMACHI
C
c  10 CALL CPUTIMT1(RTIMET)
c     RETURN
C
c  20 CALL CPUTIMT2(RTIMET)
c     RETURN
C
c  30 CALL CPUTIMT3(RTIMET)
c     RETURN
C
c  50 CALL CPUTIMT5(RTIMET)
c     RETURN
C
c  70 CALL CPUTIMT7(RTIMET)
c     RETURN
C
c  80 CALL CPUTIMT8(RTIMET)
      RETURN
C
      END
