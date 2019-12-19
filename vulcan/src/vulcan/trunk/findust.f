      SUBROUTINE FINDUST
C***********************************************************************
C
C**** THIS ROUTINE FINDS USEFUL CHARACTERS IN INPUT FILE NAME
C
C     Note: the number 80 in the DO loop is in relation with the
C           definition of input/output files in prob_omt.f
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
C
C**** FINDS USEFUL CHARACTERS
C
      DO INUS1T=1,80
       IF(CET(INUS1T:INUS1T).NE.' ') THEN
        ILUS1T=INUS1T
       END IF
      END DO
C
      RETURN
      END
