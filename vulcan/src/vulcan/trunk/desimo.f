      SUBROUTINE DESIMO(XJACM,XJA3M,DETJM,STRAP,
     .                  RCGPI,FINET,TRABE)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES (IFREN=7):
C     1) INVERSE OF THE PLASTIC RIGHT CAUCHY-GREEN TENSOR (RCGPI)
C     2) ELASTIC LEFT CAUCHY-GREEN TENSOR (FINET)
C     3) TRACE OF FINET (TRABE)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION XJACM(NDIME,*), STRAP(*)
      DIMENSION RCGPI(*),       FINET(*)
      DIMENSION RCGPT(6),       UNOMA(6)               ! auxiliar arrays
C
      IF(LARGE.EQ.3) RETURN
C
      DO ISTR1=1,NSTR1
       UNOMA(ISTR1)=0.0D0
      END DO
      UNOMA(1)=1.0D0
      UNOMA(2)=1.0D0
      UNOMA(4)=1.0D0
C
C**** COMPUTES PLASTIC RIGHT CAUCHY-GREEN TENSOR
C
      DO ISTRS=1,NSTRS
       RCGPT(ISTRS)=2.0D0*STRAP(ISTRS)+UNOMA(ISTRS)
      ENDDO
C
C**** COMPUTES INVERSE OF PLASTIC RIGHT CAUCHY-GREEN TENSOR
C
C     Note: determinant of RCGPT is 1 (plastic incompressibility)
C
      RCGPI(1)=1.0D0/RCGPT(1)                ! 1D case
      IF(NTYPE.NE.5) THEN                    ! 2D cases
       RCGPI(1)= RCGPT(2)*RCGPT(4)
       RCGPI(2)= RCGPT(1)*RCGPT(4)
       RCGPI(3)=-RCGPT(3)*RCGPT(4)
       RCGPI(4)=1.0D0/RCGPT(4)
       IF(NTYPE.EQ.4) THEN                   ! 3D case
        RCGPI(1)= (RCGPT(2)*RCGPT(4)-RCGPT(6)*RCGPT(6))
        RCGPI(2)= (RCGPT(1)*RCGPT(4)-RCGPT(5)*RCGPT(5))
        RCGPI(3)=-(RCGPT(3)*RCGPT(4)-RCGPT(5)*RCGPT(6))
        RCGPI(4)= (RCGPT(1)*RCGPT(2)-RCGPT(3)*RCGPT(3))
        RCGPI(5)= (RCGPT(3)*RCGPT(6)-RCGPT(2)*RCGPT(5))
        RCGPI(6)=-(RCGPT(1)*RCGPT(6)-RCGPT(3)*RCGPT(5))
       ENDIF           ! ntype.eq.4
      ENDIF            ! ntype.ne.5
C
C**** COMPUTES ELASTIC LEFT CAUCHY-GREEN TENSOR
C
C     Note: this tensor is computed as the push-forward (stress
C           transformation) of RCGPI
C
      DETJX=1.0D0
      CALL PUFOBA(RCGPI,FINET,XJACM,XJA3M,DETJX,    1)    ! push-forward
C
      TRABE=FINET(1)                         ! 1D case
      IF(NTYPE.NE.5)                         ! 2D & 3D cases
     . TRABE=TRABE+FINET(2)+FINET(4)
C
      RETURN
      END
