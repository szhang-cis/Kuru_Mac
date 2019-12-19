      SUBROUTINE ADDDATS(DISTOS,ELDATS,ELPRES,ELVARS)
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
      INCLUDE 'para_oms.f'
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_oms.f'
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      DIMENSION DISTOS(*), ELDATS(*), ELPRES(*), ELVARS(*)
C
C**** RECOVER ELEMENT ARRAYS FROM DATA-FILE
C
      DO IELEMS=1,NELEMS
C
C**** READ ELDAT, ELPRE & ELVAR
C
       IF(NMEMOS.EQ.1)
     .  CALL DATBASS(ELPRES,    2,    2)
       IF(NMEMO3S.EQ.1.OR.NMEMO4S.EQ.1.OR.NMEMO5S.EQ.0)
     .  CALL DATBASS(ELVARS,    3,    2)
C
C**** WRITE ELDAT, ELPRE & ELVAR
C
       IF(NMEMOS.EQ.1)
     .  CALL DATBASS(ELPRES,    4,    1)
       IF(NMEMO3S.EQ.1.OR.NMEMO4S.EQ.1.OR.NMEMO5S.EQ.0)
     .  CALL DATBASS(ELVARS,    5,    1)
      ENDDO
C
C**** TRANSFER DISTO TO DATA-BASE
C
      CALL DATBASS(DISTOS,   12,    1)
C
      RETURN
      END
