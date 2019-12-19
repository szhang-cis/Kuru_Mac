      SUBROUTINE ENERGY(DISIT,REFOR)
C**********************************************************************
C
C**** THIS ROUTINE COMPUTES GZERO AND CORRECT THE LUMPED EXPWP
C
C**********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'auxl_om.f'
      INCLUDE 'inte_om.f'
C
      DIMENSION DISIT(*), REFOR(*)
C
      GZERO=0.
      DO IDOFN=1,NDOFN
        DO IPOIN=1,NPOIN
          ITOTV=(IPOIN-1)*NDOFC+IDOFN
          GZERO=GZERO+DISIT(ITOTV)*REFOR(ITOTV)
        ENDDO
      ENDDO
      IF(KPORE.EQ.2)THEN
        DO ITOTV=NDOFC,NTOTV,NDOFC
          DISIT(ITOTV)=DISIT(ITOTV)*WLUMP
        ENDDO
      ENDIF
C
      RETURN
      END     
