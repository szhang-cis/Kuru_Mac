      SUBROUTINE TRASMEG(ELVAR,ELVART)
C***********************************************************************
C
C**** THIS ROUTINE TRANSFERS GAUSSIAN MECHANICAL COUPLING TERM FOR
C     THERMAL PROBLEM
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION ELVAR(NSTAT), ELVART(NSTATT)
C
      IF(NITERC.EQ.0.OR.NITERC.EQ.2.OR.NITERC.EQ.3) RETURN
C
C**** LOOP OVER ELEMENTS
C
      DO 1000 IELEMC=1,NELEMC
C
      IELEM=IELEMC
      IELEMT=IELEMC
C
C**** READ ELEM. DATA FROM DATA BASE
C
      CALL DATBAS(ELVAR,     3,    2)    ! current
      IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1)
     . CALL DATBAST(ELVART,    5,    2)  ! last converged
C
C**** TRANSFERS EHIST(IPLAS(6)) TO EHISTT(7+IX)
C
      CALL CAEHIST(ELVAR(ISTAT(2)),ELVART(ISTATT(2)))
C
C**** WRITE CURRENT STATE VARIABLES TO DATA BASE
C
      IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1)
     . CALL DATBAST(ELVART,    5,    1)  ! last converged
C
 1000 CONTINUE
C
      RETURN
      END
