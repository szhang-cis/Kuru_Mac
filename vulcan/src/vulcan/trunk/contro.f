      SUBROUTINE CONTRO
C***********************************************************************
C
C**** CONTROL ROUTINE FOR VULCAN
C
C***********************************************************************
C
C**** ESTABLISH THE CONTROL PARAMETERS
C
      CALL CONINP
C
C**** ESTABLISH THE PROBLEM PARAMETERS
C
      CALL PROINP
C
C**** ESTABLISH THE REMAINING PARAMETERS FOR PROBLEM
C
      CALL CONSET
C
C**** CHECK ON MAXIMUM DIMENSION AND COMPATIBILITY OF INPUT PROBLEM
C     OPTIONS
C
      CALL CONCEK
C
      END
