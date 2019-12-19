      SUBROUTINE ENERGYT(DISITT,REFORT)
C**********************************************************************
C
C**** THIS ROUTINE COMPUTES GZERO AND CORRECT THE LUMPED EXPWP
C
C**********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inte_omt.f'
C
      DIMENSION DISITT(*), REFORT(*)
C
      GZEROT=0.
      DO IDOFNT=1,NDOFNT
        DO IPOINT=1,NPOINT
          ITOTVT=(IPOINT-1)*NDOFCT+IDOFNT
          GZEROT=GZEROT+DISITT(ITOTVT)*REFORT(ITOTVT)
        ENDDO
      ENDDO
C
      RETURN
      END     
