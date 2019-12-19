      SUBROUTINE SOL132(PROPS,TEMPG,
     .                  VISCO,EXPON,
     .                  VISCOF,VISCOA,EXPONF,EXPONA,
     .                  AZABA,AZABB,AZABC,
     .                  AZABAF,AZABAA,AZABBF,AZABBA,AZABCF,AZABCA)
#ifndef restricted
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES OTHER VARIABLES OF THE VISCOPLASTIC MODEL
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
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
      IF(IVIFL.EQ.1.OR.IVIFL.EQ.2.OR.IVIFL.EQ.3.OR.IVIFL.EQ.6)
     . CALL VISCOT(VISCO,TEMPG,PROPS,
     .             VISCOF,VISCOA)
      IF(IVIFL.EQ.4.OR.IVIFL.EQ.6.OR.IVIFL.EQ.7)
     . CALL ZABART(AZABA,AZABB,AZABC,TEMPG,PROPS,
     .             AZABAF,AZABAA,AZABBF,AZABBA,AZABCF,AZABCA)
C
      IF((IVIFL.GE.1.AND.IVIFL.LE.4).OR.IVIFL.EQ.6)
     . CALL EXPONT(EXPON,TEMPG,PROPS,
     .             EXPONF,EXPONA)
C
      IF(IVIFL.EQ.5) THEN
       IF(ITERME.LT.0) CENKEL=273.15D0   +  940.0D0   !!!!
       RECRY(21)=DEXP(RECRY(3)/(RECRY(4)*(TEMPG+CENKEL)))
      ENDIF
C
#endif
      RETURN
      END
