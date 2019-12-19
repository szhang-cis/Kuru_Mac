      SUBROUTINE SOLPOR(PROPS,TEMPG,CPOEF,CPOE1,CPOE2)
#ifndef restricted
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES A VARIABLE ASSOCIATED TO THE POROSITY
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'prob_om.f'
C
      DIMENSION PROPS(*)
C
C**** MATERIAL PROPERTIES (EXPLICIT TEMPERATURE FUNCTIONS)
C
      CALL CPOEFT(CPOEF,CPOE1,CPOE2,TEMPG,PROPS)
C
#endif
      RETURN
      END
