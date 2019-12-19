      SUBROUTINE SETM04(VNORL)
C***********************************************************************
C
C**** THIS ROUTINE SETS UP SOME NEEDED CONSTANT MATRICES FOR
C     FUTURE USE ( ELEMENT NO. 4 )
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION VNORL(NDIME,NNODL)
C
C**** INITIALISES OUTWARD UNIT NORMAL (FOR CONTACT PROBLEM)
C
      DO IDIME=1,NDIME
       DO INODL=1,NNODL
        VNORL(IDIME,INODL)=0.0D+00
       ENDDO
      ENDDO
C
      RETURN
      END
