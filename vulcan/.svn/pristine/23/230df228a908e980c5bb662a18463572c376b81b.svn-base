      SUBROUTINE ASSIFIT
C***********************************************************************
C
C**** THIS ROUTINE ASSIGNS THERMAL FILES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** THERMAL FILES ALREADY DEFINED IN prob_omt.f
C
      GO TO (10,10,30,10,50,10,10,80), NMACHI
C
   10 CALL ASSIFIT1
      RETURN
C
   30 CALL ASSIFIT3
      RETURN
C
   50 CALL ASSIFIT5
      RETURN
C
   80 CALL ASSIFIT8
      RETURN
C
      END
