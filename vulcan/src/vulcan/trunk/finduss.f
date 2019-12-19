      SUBROUTINE FINDUSS
C***********************************************************************
C
C**** THIS ROUTINE FINDS USEFUL CHARACTERS IN INPUT FILE NAME
C
C     Note: the number 80 in the DO loop is in relation with the
C           definition of input/output files in prob_oms.f
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
C
C**** FINDS USEFUL CHARACTERS
C
      DO INUS1S=1,80
       IF(CES(INUS1S:INUS1S).NE.' ') THEN
        ILUS1S=INUS1S
       END IF
      END DO
C
      RETURN
      END
