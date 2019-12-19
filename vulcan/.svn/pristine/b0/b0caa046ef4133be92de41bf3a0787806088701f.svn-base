      SUBROUTINE SETI03(ELCOD,ELDIS,ELCO1,ELDI1)
C***********************************************************************
C
C**** THIS ROUTINE SETS UP SOME NEEDED CONSTANT MATRICES FOR
C     FUTURE USE ( ELEMENT NO. 4 )
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
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
      COMMON/TUNING4/RITEN,RITEF,TRUPL,TRUPM,TOLGA,TOLGAM
      COMMON/JACOBSA/IEROR,KEROR
C
      DIMENSION ELCOD(NDIME,*)
      DIMENSION ELCO1(NDIME,*)
      DIMENSION ELDIS(NDOFC,*),       ELDI1(NDOFC,*)
C
      DATA TWOPI/6.283185307179586/
C
      IF(NMEMO1M.EQ.0) THEN
       DO IDIME=1,NDIME
        DO INODL=1,NNODL
         ELCO1(IDIME,INODL)=ELCOD(IDIME,INODL)
        ENDDO
       ENDDO
      ENDIF
C
      IF(NMEMO5M.EQ.0) THEN
       DO IDOFC=1,NDOFC
        DO INODL=1,NNODL
         ELDI1(IDOFC,INODL)=ELDIS(IDOFC,INODL)
        ENDDO
       ENDDO
      ENDIF
C
      RETURN
      END
