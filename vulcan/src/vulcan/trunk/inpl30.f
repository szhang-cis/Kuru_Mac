      SUBROUTINE INPL30(PROPS,EHIST)
C***********************************************************************
C
C**** THIS ROUTINE INITIALITES SOME INTERNAL VARIABLES
C     ( FOR ELEMENT NO. 30 )
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION PROPS(*), EHIST(NHIST,*)
C
C**** LOOP OVER INTEGRATION POINTS
C
      DO IGAUS=1,NGAUL
C
C**** INITIALIZATION OF SOME PARAMETERS
C
       CALL INPLAS(PROPS,EHIST(IPLAS(3),IGAUS))
C
      END DO ! IGAUS=1,NGAUL
C
      RETURN
      END
