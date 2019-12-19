      SUBROUTINE STRATET
C***********************************************************************
C
C**** THIS ROUTINE SETS UP THE CONTROLLING PARAMETERS FOR THE
C     CURRENT TIME INTERVAL
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** READ NEW PARAMETERS 
C
      CALL STRINPT
C
C**** CHECK INPUT PARAMETERS AND DETECT INCOMPATIBILITIES
C
      CALL STRCEKT
C
C**** PRINT OUT THE DECIDED STRATEGY
C
      CALL STROUTT
C
      RETURN
      END
