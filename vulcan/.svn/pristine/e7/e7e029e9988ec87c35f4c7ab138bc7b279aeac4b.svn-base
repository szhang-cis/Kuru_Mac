      SUBROUTINE CALCFJ(XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .                  SHAPE,GPCOD,RCGTT,RCGTI,TRACC,
     .                  TRACI,ENEYE,DETJC,SEINC,FACTJ,
     .                  DCOFA,UNORA,DCINV,ROMTX,UNOMA)
#ifndef restricted
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE RIGHT CAUCHY-GREEN TENSOR & ITS INVERSE
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION XJACM(NDIME,*), XJACI(NDIME,*)
      DIMENSION SHAPE(*),       GPCOD(*)
      DIMENSION RCGTT(*),       RCGTI(*),
     .          ENEYE(*)
      DIMENSION DCOFA(6,*),     UNORA(6,*),
     .          DCINV(6,*)
      DIMENSION ROMTX(NDIME,*), D(3), ENEYX(6), UTSTR(6), UISTR(6),
     .          UNOMA(*)
C
      IF(LARGE.EQ.3) RETURN
C
      ISTAN=2
      IF(IFREN.EQ.2.OR.IFREN.EQ.3.OR.
     .   IFREN.EQ.5.OR.IFREN.EQ.6.OR.IFREN.EQ.7.OR.IFREN.EQ.8.OR.
     .   IFREN.EQ.9.OR.IFREN.EQ.10.OR.
     .   IFREN.EQ.51.OR.IFREN.EQ.52.OR.IFREN.EQ.53.OR.IFREN.EQ.54.OR.
     .   IFREN.EQ.55.OR.IFREN.EQ.56.OR.IFREN.EQ.57.OR.IFREN.EQ.58)
     . ISTAN=1
      IF(ISTAN.EQ.2) RETURN
C
C**** COMPUTES RIGHT CAUCHY-GREEN TENSOR, ITS TRACE & ITS INVERSE
C
C     Notes:
C
C     TRACI=SXJAI (trace of C^-1 = second invariant of F^-1)
C
C     An equivalent form to compute the inverse of C is: F^-1 F^-T
C     For example, for the axisymmetric case C^-1 is:
C     RCGTI(1)=XJACI(1,1)*XJACI(1,1)+XJACI(1,2)*XJACI(1,2)
C     RCGTI(2)=XJACI(2,1)*XJACI(2,1)+XJACI(2,2)*XJACI(2,2)
C     RCGTI(3)=XJACI(1,1)*XJACI(2,1)+XJACI(1,2)*XJACI(2,2)
C     RCGTI(4)=XJA3I*XJA3I ( = 1.0/(XJA3M*XJA3M) )
C
C     Invariants of C:
C     TRACC=trace of C
C     SEINC=0.5(TRACC**2-C:C)
C     DETJC=determinant of C
C
C     FACTJ=determinant of F to the power -2/3
C     DCOFA=derivative of cofactor(C_ij) with respect to C_kl
C           where cofactor(C_ij)=J*J*C^-1_ij
C     UNORA=delta_ik*delta_jl. Useful to compute the second derivative
C           of SEINC with respect to C = UNOMA*UNOMA-UNORA.
C     DCINV=derivative of C^-1_i with respect to C_kl multiplied by -2
C
C     ROMTX=rotation matrix
C
      FACTJ=DETJM**(-2.0D0/3.0D0)            ! J^(-2/3)
      DETJC=DETJM*DETJM                      ! determinant of C
C
      RCGTT(1)=XJACM(1,1)*XJACM(1,1)         ! right Cauchy-Green tensor
      TRACC=RCGTT(1)                         ! its trace
      SEINC=0.5D0*(TRACC*TRACC-RCGTT(1)*RCGTT(1))
      RCGTI(1)=1.0D0/RCGTT(1)                ! its inverse
      TRACI=RCGTI(1)                         ! trace of the inverse
      ENEYE(1)=RCGTI(1)*RCGTI(1)             ! C^-1 * C^-1
      ENEYX(1)=RCGTT(1)*RCGTT(1)             ! C * C
      SXJAI=XJACI(1,1)*XJACI(1,1)            ! second invariant of F^-1
      DCOFA(1,1)=0.0D0                       ! C deriv. of cofactor(C)
      UNORA(1,1)=0.0D0
      DCINV(1,1)=2.0D0*RCGTI(1)*RCGTI(1)
      IF(NTYPE.NE.5) THEN
       RCGTT(1)=XJACM(1,1)*XJACM(1,1)+XJACM(2,1)*XJACM(2,1)
       RCGTT(2)=XJACM(1,2)*XJACM(1,2)+XJACM(2,2)*XJACM(2,2)
       RCGTT(3)=XJACM(1,1)*XJACM(1,2)+XJACM(2,1)*XJACM(2,2)
       RCGTT(4)=1.0D0
       TRACC=RCGTT(1)+RCGTT(2)
       SEINC=0.5D0*(TRACC*TRACC-RCGTT(1)*RCGTT(1)-RCGTT(2)*RCGTT(2)-
     .                    2.0D0*RCGTT(3)*RCGTT(3))
       DCOFA(1,2)=RCGTT(4)
       DCOFA(1,3)=0.0D0
       DCOFA(1,4)=RCGTT(2)
       DCOFA(2,2)=0.0D0
       DCOFA(2,3)=0.0D0
       DCOFA(2,4)=RCGTT(1)
       DCOFA(3,3)=-RCGTT(4)
       DCOFA(3,4)=-RCGTT(3)
       DCOFA(4,4)=0.0D0
       UNORA(1,2)=0.0D0
       UNORA(1,3)=0.0D0
       UNORA(1,4)=0.0D0
       UNORA(2,2)=0.0D0
       UNORA(2,3)=0.0D0
       UNORA(2,4)=0.0D0
       UNORA(3,3)=1.0D0
       UNORA(3,4)=0.0D0
       UNORA(4,4)=0.0D0
       IF(NTYPE.EQ.2) THEN
        TRACC=TRACC+RCGTT(4)
        SEINC=0.5D0*(TRACC*TRACC-RCGTT(1)*RCGTT(1)-RCGTT(2)*RCGTT(2)-
     .                        2.0D0*RCGTT(3)*RCGTT(3)-RCGTT(4)*RCGTT(4))
       ENDIF
       RCGTI(1)= RCGTT(2)/DETJC
       RCGTI(2)= RCGTT(1)/DETJC
       RCGTI(3)=-RCGTT(3)/DETJC
       RCGTI(4)=1.0D0
       TRACI=RCGTI(1)+RCGTI(2)
       IF(NTYPE.EQ.2) TRACI=TRACI+RCGTI(4)
       ENEYE(1)=RCGTI(1)*RCGTI(1)+RCGTI(3)*RCGTI(3)
       ENEYE(2)=RCGTI(2)*RCGTI(2)+RCGTI(3)*RCGTI(3)
       ENEYE(3)=RCGTI(1)*RCGTI(3)+RCGTI(2)*RCGTI(3)
       ENEYE(4)=RCGTI(4)*RCGTI(4)
       ENEYX(1)=RCGTT(1)*RCGTT(1)+RCGTT(3)*RCGTT(3)
       ENEYX(2)=RCGTT(2)*RCGTT(2)+RCGTT(3)*RCGTT(3)
       ENEYX(3)=RCGTT(1)*RCGTT(3)+RCGTT(2)*RCGTT(3)
       ENEYX(4)=RCGTT(4)*RCGTT(4)
       SXJAI=XJACI(1,1)*XJACI(1,1)+XJACI(1,2)*XJACI(1,2)+
     .       XJACI(2,1)*XJACI(2,1)+XJACI(2,2)*XJACI(2,2)
       IF(NTYPE.EQ.2) SXJAI=SXJAI+XJA3I*XJA3I
       DCINV(1,1)=2.0D0*RCGTI(1)*RCGTI(1)
       DCINV(1,2)=2.0D0*RCGTI(3)*RCGTI(3)
       DCINV(1,3)=2.0D0*RCGTI(1)*RCGTI(3)
       DCINV(1,4)=0.0D0
       DCINV(2,2)=2.0D0*RCGTI(2)*RCGTI(2)
       DCINV(2,3)=2.0D0*RCGTI(2)*RCGTI(3)
       DCINV(2,4)=0.0D0
       DCINV(3,3)=RCGTI(1)*RCGTI(2)+RCGTI(3)*RCGTI(3)
       DCINV(3,4)=0.0D0
       DCINV(4,4)=2.0D0*RCGTI(4)*RCGTI(4)
       IF(NTYPE.EQ.3) THEN
        RCGTT(4)=XJA3M*XJA3M
        TRACC=TRACC+RCGTT(4)
        SEINC=0.5D0*(TRACC*TRACC-RCGTT(1)*RCGTT(1)-RCGTT(2)*RCGTT(2)-
     .                        2.0D0*RCGTT(3)*RCGTT(3)-RCGTT(4)*RCGTT(4))
        RCGTI(1)= RCGTT(2)*RCGTT(4)/DETJC
        RCGTI(2)= RCGTT(1)*RCGTT(4)/DETJC
        RCGTI(3)=-RCGTT(3)*RCGTT(4)/DETJC
        RCGTI(4)=1.0D0/RCGTT(4)
        TRACI=RCGTI(1)+RCGTI(2)+RCGTI(4)
        ENEYE(1)=RCGTI(1)*RCGTI(1)+RCGTI(3)*RCGTI(3)
        ENEYE(2)=RCGTI(2)*RCGTI(2)+RCGTI(3)*RCGTI(3)
        ENEYE(3)=RCGTI(1)*RCGTI(3)+RCGTI(2)*RCGTI(3)
        ENEYE(4)=RCGTI(4)*RCGTI(4)
        ENEYX(1)=RCGTT(1)*RCGTT(1)+RCGTT(3)*RCGTT(3)
        ENEYX(2)=RCGTT(2)*RCGTT(2)+RCGTT(3)*RCGTT(3)
        ENEYX(3)=RCGTT(1)*RCGTT(3)+RCGTT(2)*RCGTT(3)
        ENEYX(4)=RCGTT(4)*RCGTT(4)
        SXJAI=XJACI(1,1)*XJACI(1,1)+XJACI(1,2)*XJACI(1,2)+
     .        XJACI(2,1)*XJACI(2,1)+XJACI(2,2)*XJACI(2,2)+
     .        XJA3I*XJA3I
        DCOFA(1,2)=RCGTT(4)
        DCOFA(3,3)=-RCGTT(4)
        DCINV(1,1)=2.0D0*RCGTI(1)*RCGTI(1)
        DCINV(1,2)=2.0D0*RCGTI(3)*RCGTI(3)
        DCINV(1,3)=2.0D0*RCGTI(1)*RCGTI(3)
        DCINV(1,4)=0.0D0
        DCINV(2,2)=2.0D0*RCGTI(2)*RCGTI(2)
        DCINV(2,3)=2.0D0*RCGTI(2)*RCGTI(3)
        DCINV(2,4)=0.0D0
        DCINV(3,3)=RCGTI(1)*RCGTI(2)+RCGTI(3)*RCGTI(3)
        DCINV(3,4)=0.0D0
        DCINV(4,4)=2.0D0*RCGTI(4)*RCGTI(4)
       ENDIF
       IF(NTYPE.EQ.4) THEN
        RCGTT(1)=XJACM(1,1)*XJACM(1,1)+XJACM(2,1)*XJACM(2,1)+
     .           XJACM(3,1)*XJACM(3,1)
        RCGTT(2)=XJACM(1,2)*XJACM(1,2)+XJACM(2,2)*XJACM(2,2)+
     .           XJACM(3,2)*XJACM(3,2)
        RCGTT(3)=XJACM(1,1)*XJACM(1,2)+XJACM(2,1)*XJACM(2,2)+
     .           XJACM(3,1)*XJACM(3,2)
        RCGTT(4)=XJACM(1,3)*XJACM(1,3)+XJACM(2,3)*XJACM(2,3)+
     .           XJACM(3,3)*XJACM(3,3)
        RCGTT(5)=XJACM(1,1)*XJACM(1,3)+XJACM(2,1)*XJACM(2,3)+
     .           XJACM(3,1)*XJACM(3,3)
        RCGTT(6)=XJACM(1,2)*XJACM(1,3)+XJACM(2,2)*XJACM(2,3)+
     .           XJACM(3,2)*XJACM(3,3)
        TRACC=RCGTT(1)+RCGTT(2)+RCGTT(4)
        SEINC=0.5D0*(TRACC*TRACC-RCGTT(1)*RCGTT(1)-RCGTT(2)*RCGTT(2)-
     .                       2.0D0*RCGTT(3)*RCGTT(3)-RCGTT(4)*RCGTT(4)-
     .                  2.0D0*RCGTT(5)*RCGTT(5)-2.0D0*RCGTT(6)*RCGTT(6))
        RCGTI(1)= (RCGTT(2)*RCGTT(4)-RCGTT(6)*RCGTT(6))/DETJC
        RCGTI(2)= (RCGTT(1)*RCGTT(4)-RCGTT(5)*RCGTT(5))/DETJC
        RCGTI(3)=-(RCGTT(3)*RCGTT(4)-RCGTT(5)*RCGTT(6))/DETJC
        RCGTI(4)= (RCGTT(1)*RCGTT(2)-RCGTT(3)*RCGTT(3))/DETJC
        RCGTI(5)= (RCGTT(3)*RCGTT(6)-RCGTT(2)*RCGTT(5))/DETJC
        RCGTI(6)=-(RCGTT(1)*RCGTT(6)-RCGTT(3)*RCGTT(5))/DETJC
        TRACI=RCGTI(1)+RCGTI(2)+RCGTI(4)
        ENEYE(1)=RCGTI(1)*RCGTI(1)+RCGTI(3)*RCGTI(3)+
     .           RCGTI(5)*RCGTI(5)
        ENEYE(2)=RCGTI(2)*RCGTI(2)+RCGTI(3)*RCGTI(3)+
     .           RCGTI(6)*RCGTI(6)
        ENEYE(3)=RCGTI(1)*RCGTI(3)+RCGTI(2)*RCGTI(3)+
     .           RCGTI(5)*RCGTI(6)
        ENEYE(4)=RCGTI(4)*RCGTI(4)+RCGTI(5)*RCGTI(5)+
     .           RCGTI(6)*RCGTI(6)
        ENEYE(5)=RCGTI(1)*RCGTI(5)+RCGTI(3)*RCGTI(6)+
     .           RCGTI(4)*RCGTI(5)
        ENEYE(6)=RCGTI(2)*RCGTI(6)+RCGTI(3)*RCGTI(5)+
     .           RCGTI(4)*RCGTI(6)
        ENEYX(1)=RCGTT(1)*RCGTT(1)+RCGTT(3)*RCGTT(3)+
     .           RCGTT(5)*RCGTT(5)
        ENEYX(2)=RCGTT(2)*RCGTT(2)+RCGTT(3)*RCGTT(3)+
     .           RCGTT(6)*RCGTT(6)
        ENEYX(3)=RCGTT(1)*RCGTT(3)+RCGTT(2)*RCGTT(3)+
     .           RCGTT(5)*RCGTT(6)
        ENEYX(4)=RCGTT(4)*RCGTT(4)+RCGTT(5)*RCGTT(5)+
     .           RCGTT(6)*RCGTT(6)
        ENEYX(5)=RCGTT(1)*RCGTT(5)+RCGTT(3)*RCGTT(6)+
     .           RCGTT(4)*RCGTT(5)
        ENEYX(6)=RCGTT(2)*RCGTT(6)+RCGTT(3)*RCGTT(5)+
     .           RCGTT(4)*RCGTT(6)
        SXJAI=XJACI(1,1)*XJACI(1,1)+XJACI(1,2)*XJACI(1,2)+
     .        XJACI(1,3)*XJACI(1,3)+
     .        XJACI(2,1)*XJACI(2,1)+XJACI(2,2)*XJACI(2,2)+
     .        XJACI(2,3)*XJACI(2,3)+
     .        XJACI(3,1)*XJACI(3,1)+XJACI(3,2)*XJACI(3,2)+
     .        XJACI(3,3)*XJACI(3,3)
        DCOFA(1,2)=RCGTT(4)
        DCOFA(1,3)=0.0D0
        DCOFA(1,4)=RCGTT(2)
        DCOFA(2,2)=0.0D0
        DCOFA(2,3)=0.0D0
        DCOFA(2,4)=RCGTT(1)
        DCOFA(3,3)=-RCGTT(4)
        DCOFA(3,4)=-RCGTT(3)
        DCOFA(4,4)=0.0D0
        DCOFA(1,5)=0.0D0
        DCOFA(1,6)=-RCGTT(6)
        DCOFA(2,5)=-RCGTT(5)
        DCOFA(2,6)=0.0D0
        DCOFA(3,5)=RCGTT(6)
        DCOFA(3,6)=RCGTT(5)
        DCOFA(4,5)=0.0D0
        DCOFA(4,6)=0.0D0
        DCOFA(5,5)=-RCGTT(2)
        DCOFA(5,6)=RCGTT(3)
        DCOFA(6,6)=-RCGTT(1)
        UNORA(1,5)=0.0D0
        UNORA(1,6)=0.0D0
        UNORA(2,5)=0.0D0
        UNORA(2,6)=0.0D0
        UNORA(3,5)=0.0D0
        UNORA(3,6)=0.0D0
        UNORA(4,5)=0.0D0
        UNORA(4,6)=0.0D0
        UNORA(5,5)=1.0D0
        UNORA(5,6)=0.0D0
        UNORA(6,6)=1.0D0
        DCINV(1,1)=2.0D0*RCGTI(1)*RCGTI(1)
        DCINV(1,2)=2.0D0*RCGTI(3)*RCGTI(3)
        DCINV(1,3)=2.0D0*RCGTI(1)*RCGTI(3)
        DCINV(1,4)=2.0D0*RCGTI(5)*RCGTI(5)
        DCINV(1,5)=2.0D0*RCGTI(1)*RCGTI(5)
        DCINV(1,6)=2.0D0*RCGTI(3)*RCGTI(5)
        DCINV(2,2)=2.0D0*RCGTI(2)*RCGTI(2)
        DCINV(2,3)=2.0D0*RCGTI(2)*RCGTI(3)
        DCINV(2,4)=2.0D0*RCGTI(6)*RCGTI(6)
        DCINV(2,5)=2.0D0*RCGTI(3)*RCGTI(6)
        DCINV(2,6)=2.0D0*RCGTI(2)*RCGTI(6)
        DCINV(3,3)=RCGTI(1)*RCGTI(2)+RCGTI(3)*RCGTI(3)
        DCINV(3,4)=2.0D0*RCGTI(5)*RCGTI(6)
        DCINV(3,5)=RCGTI(1)*RCGTI(6)+RCGTI(3)*RCGTI(5)
        DCINV(3,6)=RCGTI(3)*RCGTI(6)+RCGTI(2)*RCGTI(5)
        DCINV(4,4)=2.0D0*RCGTI(4)*RCGTI(4)
        DCINV(4,5)=2.0D0*RCGTI(4)*RCGTI(5)
        DCINV(4,6)=2.0D0*RCGTI(4)*RCGTI(6)
        DCINV(5,5)=RCGTI(1)*RCGTI(4)+RCGTI(5)*RCGTI(5)
        DCINV(5,6)=RCGTI(3)*RCGTI(4)+RCGTI(5)*RCGTI(6)
        DCINV(6,6)=RCGTI(2)*RCGTI(4)+RCGTI(6)*RCGTI(6)
       ENDIF           ! ntype.eq.4
      ENDIF            ! ntype.ne.5
C
      DO ISTRS=1,NSTRS           ! not strictly necessary (see fourth.f)
       DO JSTRS=ISTRS,NSTRS
        DCOFA(JSTRS,ISTRS)=DCOFA(ISTRS,JSTRS)
        UNORA(JSTRS,ISTRS)=UNORA(ISTRS,JSTRS)
        DCINV(JSTRS,ISTRS)=DCINV(ISTRS,JSTRS)
       ENDDO
      ENDDO
C
C**** COMPUTES ROTATION TENSOR (ROMTX)
C
C     Note: closed-form evaluation of the rotation tensor by Simo
C           (not used now!)
C
      CALL PRINSI(RCGTT,D)                ! eigenvalues of C
      AL1=DSQRT(D(1))                     ! eigenvlaues of U
      AL2=DSQRT(D(2))
      AL3=DSQRT(D(3))
      AIN1U=AL1+AL2+AL3                   ! invariants of U
      AIN2U=AL1*AL2+AL1*AL3+AL2*AL3
      AIN3U=AL1*AL2*AL3
      AD=AIN1U*AIN2U-AIN3U                ! denominator
      DO ISTRS=1,NSTRS                    ! U & U^-1
       UTSTR(ISTRS)=1.0D0/AD*(-ENEYX(ISTRS)+
     .              (AIN1U*AIN1U-AIN2U)*RCGTT(ISTRS)+
     .              AIN1U*AIN3U*UNOMA(ISTRS))
       UISTR(ISTRS)=1.0D0/AIN3U*(RCGTT(ISTRS)-AIN1U*UTSTR(ISTRS)+
     .              AIN2U*UNOMA(ISTRS))
      ENDDO
      ROMTX(1,1)=XJACM(1,1)*UISTR(1)+XJACM(1,2)*UISTR(3)
      ROMTX(1,2)=XJACM(1,1)*UISTR(3)+XJACM(1,2)*UISTR(2)
      ROMTX(2,1)=XJACM(2,1)*UISTR(1)+XJACM(2,2)*UISTR(3)
      ROMTX(2,2)=XJACM(2,1)*UISTR(3)+XJACM(2,2)*UISTR(2)
      IF(NDIME.EQ.3) THEN
       ROMTX(1,1)=ROMTX(1,1)+XJACM(1,3)*UISTR(5)
       ROMTX(1,2)=ROMTX(1,2)+XJACM(1,3)*UISTR(6)
       ROMTX(1,3)=XJACM(1,1)*UISTR(5)+XJACM(1,2)*UISTR(6)+
     .            XJACM(1,3)*UISTR(4)
       ROMTX(2,1)=ROMTX(2,1)+XJACM(2,3)*UISTR(5)
       ROMTX(2,2)=ROMTX(2,2)+XJACM(2,3)*UISTR(6)
       ROMTX(2,3)=XJACM(2,1)*UISTR(5)+XJACM(2,2)*UISTR(6)+
     .            XJACM(2,3)*UISTR(4)
       ROMTX(3,1)=XJACM(3,1)*UISTR(1)+XJACM(3,2)*UISTR(3)+
     .            XJACM(3,3)*UISTR(5)
       ROMTX(3,2)=XJACM(3,1)*UISTR(3)+XJACM(3,2)*UISTR(2)+
     .            XJACM(3,3)*UISTR(6)
       ROMTX(3,3)=XJACM(3,1)*UISTR(5)+XJACM(3,2)*UISTR(6)+
     .            XJACM(3,3)*UISTR(4)
      ENDIF
C
      RETURN
#endif
      END
