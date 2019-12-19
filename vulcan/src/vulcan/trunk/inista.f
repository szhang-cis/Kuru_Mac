      SUBROUTINE INISTA(PWORK,PREAS,TGAPS)
C***********************************************************************
C
C**** THIS ROUTINE INITIALISES COUPLING & CONTACT MECHANICAL VARIABLES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
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
      INCLUDE 'inpo_om.f'
C
      DIMENSION PWORK(NPOIN,2), PREAS(NPOIN), TGAPS(NPOIN)
C
C**** INITIALISE COUPLING MECHANICAL TERM, NORMAL PRESSURE & NORMAL GAP
C
      IF(ITERME.GT.0) THEN
       DO IPOIN=1,NPOIN
        PWORK(IPOIN,1)=0.0D0
       ENDDO
      ENDIF
C
      IF(KGAPC.EQ.1) THEN
       DO IPOIN=1,NPOIN
        PREAS(IPOIN)=0.0D0
        TGAPS(IPOIN)=0.0D0
       ENDDO
      ENDIF
C
      RETURN
      END
