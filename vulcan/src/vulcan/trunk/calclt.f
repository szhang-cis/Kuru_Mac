      SUBROUTINE CALCLT(DMATX,DMAT1,ESTAB)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE LIQUID CONSTITUTIVE TENSOR
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
      INCLUDE 'auxl_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'prob_om.f'
C
      DIMENSION DMATX(NSTRS,*), DMAT1(6,6)
C
      IF(LARGE.EQ.3) RETURN
C
      ISTAN=2
      IF(LARGE.EQ.0) THEN
       ISTAN=1
      ELSE
       IF(IFREN.EQ.1.OR.IFREN.EQ.2.OR.IFREN.EQ.3.OR.
     .    IFREN.EQ.4.OR.IFREN.EQ.5.OR.IFREN.EQ.6.OR.IFREN.EQ.8.OR.
     .    IFREN.EQ.9.OR.IFREN.EQ.10) ISTAN=1
      ENDIF
C
C**** COMPUTES THE CONSTITUTIVE TENSOR IN THE LIQUID
C
      IF(ISTAN.EQ.1) THEN
       DMAT1(3,3)=DMATX(1,1)/ESTAB
       IF(NTYPE.EQ.4) THEN
        DMAT1(5,5)=DMATX(1,1)/ESTAB
        DMAT1(6,6)=DMATX(1,1)/ESTAB
       ENDIF
      ELSE                      ! (istan=2)
       IF(IFREN.EQ.7) THEN
        CALL RUNEND('ERROR: MODEL 7 NOT IMPLEMENTED IN CALCCL')
       ENDIF
      ENDIF                     ! istan.eq.1
C
      RETURN
      END
