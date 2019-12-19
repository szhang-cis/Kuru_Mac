      SUBROUTINE PLANES(STRAN,STRAP)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE "ZZ" COMPONENT OF THE DEFORMATION
C     TENSOR FOR THE PLANE STRESS CASE
C     (this component is not involved in the integration algorithm of 
C     the plastic evolution equations, see papers of Simo, etc.)
C
C     Notes:
C
C     This computation is intented to be LARGE index independent
C
C     STRAN(4)=-(PXZ/(1.0D00-PXZ))*STRAN(1)               ! orthotropic
C    .         -(PYZ/(1.0D00-PYZ))*STRAN(2)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'prob_om.f'
C
      COMMON/PROPSC/POISC1,POISC2,POISC3,STRRA,STRR4,DAMAP1,DAMAP2
C
      DIMENSION STRAN(*), STRAP(*), STRAT(4)
C
C**** THERMAL DEFORMATION
C
      STRAT(1)=STRRA
      STRAT(2)=STRRA
      STRAT(3)=0.0D0
      STRAT(4)=STRRA
C
      STRAN(4)=POISC1*(STRAN(1)-STRAP(1)-STRAT(1))+
     .         POISC2*(STRAN(2)-STRAP(2)-STRAT(2))+
     .         POISC3*(STRAN(3)-STRAP(3)-STRAT(3))+
     .         STRAP(4)+STRAT(4)+STRR4
C
      RETURN
      END
