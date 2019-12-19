      SUBROUTINE SOLDAM(PROPS,TEMPG,
     .                  CDAMA,CFRAC)
#ifndef restricted
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES OTHER VARIABLES OF THE DAMAGE MODEL
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
      CALL CDAMAT(CDAMA,TEMPG,PROPS)
      CALL CFRACT(CFRAC,TEMPG,PROPS)
C
#endif
      RETURN
      END
