      SUBROUTINE SETI04(VNORL,ELCOD,ELDIS,ELCO1,ELDI1,ITASK)
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
      DIMENSION VNORL(NDIME,NNODL)
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
      IF(NMEMO2M.EQ.0) THEN
       IF(ITASK.EQ.1) THEN
        CALL SETM04(VNORL)
       ELSE
        RETURN
       ENDIF               ! itask.eq.1
      ELSE
C
C**** NMEMO2M=1 >> CONTACT PROBLEM DEPENDENT OF ITASK
C
C     ITASK=1: initialises VNORL
C     ITASK>1: nothing
C
       IF(ITASK.EQ.1) THEN
        CALL SETM04(VNORL)
       ENDIF               ! itask.eq.1
      ENDIF                ! nmemo2m.eq.0
C
      RETURN
      END
