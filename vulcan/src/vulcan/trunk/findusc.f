      SUBROUTINE FINDUSC
C***********************************************************************
C
C**** THIS ROUTINE FINDS USEFUL CHARACTERS IN INPUT FILE NAME
C
C     Note: the number 80 in the DO loop is in relation with the
C           definition of input/output files in nuec_om.f
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** FINDS USEFUL CHARACTERS
C
      DO INUS1C=1,80
       IF(CCOA(INUS1C:INUS1C).NE.' ') THEN
        ILUS1C=INUS1C
       END IF
      END DO
C
      RETURN
      END
