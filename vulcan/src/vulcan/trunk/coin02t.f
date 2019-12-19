      SUBROUTINE COIN02T(TEMP2,TEMP1,TEMPL,RCOOR,ISOLI)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE INTERFACE COORDINATES
C     IN A BIPHASE ELEMENT FOR THE CASES:
C
C     3-NODED BIDIMENSIONAL ELEMENT (IN ONE BOUNDARY)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      TETOL=1.0D-10
      TEMDI=TEMP2-TEMP1
      TEMDI=DABS(TEMDI)
      IF(TEMDI.GT.TETOL)THEN
       RCOOR=(TEMPL-TEMP1)/(TEMP2-TEMP1)      
       ISOLI=0
       IF(RCOOR.LE.1.0D0.AND.RCOOR.GE.0.0D0) ISOLI=1
      ELSE
       ISOLI=0
      ENDIF
C
      RETURN
      END
