      SUBROUTINE CONTROS
C***********************************************************************
C
C**** CONTROL ROUTINE FOR VULCAN
C
C***********************************************************************
C
C**** ESTABLISH THE CONTROL PARAMETERS
C
      CALL CONINPS
C
C**** ESTABLISH THE PROBLEM PARAMETERS
C
      CALL PROINPS
C
C**** ESTABLISH THE REMAINING PARAMETERS FOR PROBLEM
C
      CALL CONSETS
C
C**** CHECK ON MAXIMUM DIMENSION AND COMPATIBILITY OF INPUT PROBLEM
C     OPTIONS
C
      CALL CONCEKS
C
      END
