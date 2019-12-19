      SUBROUTINE ADDDAT(DISTO,ELDAT,ELPRE,ELVAR)
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
      INCLUDE 'para_om.f'
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION DISTO(*), ELDAT(*), ELPRE(*), ELVAR(*)
C
C**** RECOVER ELEMENT ARRAYS FROM DATA-BASE
C
      DO IELEM=1,NELEM
C
C**** READ ELDAT, ELPRE & ELVAR
C
       IF(NMEMOM.EQ.1.OR.INITV.EQ.1)
     .  CALL DATBAS(ELPRE,    2,    2)
       CALL DATBAS(ELVAR,    3,    2)
C
C**** WRITE ELDAT, ELPRE & ELVAR
C
       IF(NMEMOM.EQ.1.OR.INITV.EQ.1)
     .  CALL DATBAS(ELPRE,    4,    1)
       CALL DATBAS(ELVAR,    5,    1)
      ENDDO
C
C**** TRANSFER DISTO TO DATABASE
C
      CALL DATBAS(DISTO,   12,    1)
C
c     WRITE(LURES,901) IBLOC
c     WRITE(LUPRI,901) IBLOC
C
      RETURN
c 901 FORMAT(//,5X,'NUMBER OF ELEMENT ARRAYS ALLOCATED IN VM:',I3,/)
      END
