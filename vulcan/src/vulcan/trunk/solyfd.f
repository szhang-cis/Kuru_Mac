      SUBROUTINE SOLYFD(PROPS,TEMPG,
     .                  CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,CCOB7,CCOB8,
     .                  CCOB9,CCOB10,CCOB11,CCOB12,CLANK,
     .                  CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES SOME VARIABLES OF THE YIELD FUNCTION, FLOW
C     POTENTIAL & DAMAGE MODEL
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
      DIMENSION CLAN1(2),       CLAN2(2),       CLAN3(2),
     .          CLAN4(2),       CLAN5(2),       CLAN6(2)
C
C**** MATERIAL PROPERTIES (EXPLICIT TEMPERATURE FUNCTIONS)
C
      CALL CCOBTT(CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,CCOB7,CCOB8,
     .            CCOB9,CCOB10,CCOB11,CCOB12,CLANK,TEMPG,PROPS,
     .            CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6)
C
      RETURN
      END
