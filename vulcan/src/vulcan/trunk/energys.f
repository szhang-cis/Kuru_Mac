      SUBROUTINE ENERGYS(DISITS,REFORS)
C**********************************************************************
C
C**** THIS ROUTINE COMPUTES GZERO AND CORRECT THE LUMPED EXPWP
C
C**********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'auxl_oms.f'
      INCLUDE 'inte_oms.f'
C
      DIMENSION DISITS(*), REFORS(*)
C
      GZEROS=0.
      DO IDOFNS=1,NDOFNS
        DO IPOINS=1,NPOINS
          ITOTVS=(IPOINS-1)*NDOFCS+IDOFNS
          GZEROS=GZEROS+DISITS(ITOTVS)*REFORS(ITOTVS)
        ENDDO
      ENDDO
C
      RETURN
      END
