      SUBROUTINE SOLISO(PROPS,TEMPG,
     .                  CCOEF,CCOEB,CCOEQ,
     .                  CCOEFF,CCOEFA,CCOEBF,CCOEBA,CCOEQF,CCOEQA,
     .                  CKOEF,CKOEB,CKOEQ,
     .                  CKOEFF,CKOEFA,CKOEBF,CKOEBA,CKOEQF,CKOEQA)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE VARIABLES ASSOCIATED WITH THE ISOTROPIC
C     & KINEMATIC HARDENINGS
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
      CALL CCOEFT(CCOEF,CCOEB,CCOEQ,
     .            TEMPG,PROPS,
     .            CCOEFF,CCOEFA,CCOEBF,CCOEBA,CCOEQF,CCOEQA)
C
      CALL CKOEFT(CKOEF,CKOEB,CKOEQ,
     .            TEMPG,PROPS,
     .            CKOEFF,CKOEFA,CKOEBF,CKOEBA,CKOEQF,CKOEQA)
C
      RETURN
      END
