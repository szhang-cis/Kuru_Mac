      SUBROUTINE FINDUS
C***********************************************************************
C
C**** THIS ROUTINE FINDS USEFUL CHARACTERS IN INPUT FILE NAME
C
C     Note: the number 80 in the DO loop is in relation with the
C           definition of input/output files in prob_om.f
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** MCHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
C
C**** FINDS USEFUL CHARACTERS
C
      DO INUS1=1,80
       IF(CE(INUS1:INUS1).NE.' ') THEN
        ILUS1=INUS1
       END IF
      END DO
C
      RETURN
      END
