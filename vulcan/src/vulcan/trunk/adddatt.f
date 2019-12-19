      SUBROUTINE ADDDATT(DISTOT,ELDATT,ELPRET,ELVART)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE MEMORY REQUIREMENTS & ADDRESS FOR 
C     DATABASE
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** GENERAL DIMENSIONS
C
      INCLUDE 'para_omt.f'
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION DISTOT(*), ELDATT(*), ELPRET(*), ELVART(*)
C
C**** RECOVER ELEMENT ARRAYS FROM DATA-FILE
C
      DO IELEMT=1,NELEMT
C
C**** READ ELDAT, ELPRE & ELVAR
C
       IF(NMEMO.EQ.1)
     .  CALL DATBAST(ELPRET,    2,    2)
       IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
     .  CALL DATBAST(ELVART,    3,    2)
C
C**** WRITE ELDAT, ELPRE & ELVAR
C
       IF(NMEMO.EQ.1)
     .  CALL DATBAST(ELPRET,    4,    1)
       IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
     .  CALL DATBAST(ELVART,    5,    1)
      ENDDO
C
C**** TRANSFER DISTO TO DATA-BASE
C
      CALL DATBAST(DISTOT,   12,    1)
C
      RETURN
      END
