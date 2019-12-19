      SUBROUTINE DMGMTX(DMATX,DAMAG)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE D-MATRIX FOR DAMAGE MATERIALS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'auxl_om.f'
C
      DIMENSION DMATX(NSTRS,*)
C
      REDUC=1.0D00-DAMAG
      DO 15 ISTRE=1,NSTRS
      DO 15 JSTRE=1,NSTRS
   15 DMATX(ISTRE,JSTRE)=REDUC*DMATX(ISTRE,JSTRE)
C
      RETURN
      END
