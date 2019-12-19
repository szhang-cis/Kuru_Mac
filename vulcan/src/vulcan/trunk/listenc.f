      SUBROUTINE LISTENC(SUBNA,NPRIN,ITAPE)
C***********************************************************************
C
C**** READS A STRING AND INTERPRETS IT AS WORDS AND PARAMETERS.
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'inpo_om.f'
C
      CHARACTER SUBNA*7
C
      GO TO (1,2,3,2,5,1,1,8), NMACHIM
C
    1 CALL LISTENC1(SUBNA,NPRIN,ITAPE)
      RETURN
C
    2 CALL LISTENC2(SUBNA,NPRIN,ITAPE)
      RETURN
C
    3 CALL LISTENC3(SUBNA,NPRIN,ITAPE)
      RETURN
C
    5 CALL LISTENC5(SUBNA,NPRIN,ITAPE)
      RETURN
C
    8 CALL LISTENC8(SUBNA,NPRIN,ITAPE)
      RETURN
C
      END
