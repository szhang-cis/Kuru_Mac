      SUBROUTINE SOURCB(QHEATT,HTLODT,ELCO1T,ELCOIT,GPCO1T,GPCODT,
     .                  SHAPET)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE TIME-SPACE-DEPENDENT SURFACE HEAT FLUX
C     ON BOUNDARY ELEMENTS 1D/2D ( ELEMENT NO. 101 )
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'   ! thermal-mechanical
      INCLUDE 'nued_om.f'   ! thermal-microstructural
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION HTLODT(NHLODT,NSUBFT,NFUNCT)
      DIMENSION ELCO1T(NDIMET,*), ELCOIT(NDIMET,*),
     .          GPCO1T(NDIMET),   GPCODT(NDIMET),   SHAPET(NNODLT)
      DIMENSION COORL(3),         VELOL(3),         RADIA(3)
      DIMENSION COORC(3), X1(3), X2(3), X3(3), XB(3), XBI(3), TB(3,3)
C
C**** COORDINATES AT GAUSS POINT AT MATERIAL OR SPATIAL CONFIGURATIONS
C
C     Notes: for ITYPL=101-104 (circular beam) and ITYPL=121-124
C            (ellipsoidal beam) the material configuration is used
C            (this simplificative assumption is considered even when
C             the energy equation is computed in the spatial
C             configuration, i.e., LARGET=3)
C
C            ELCOIT: material coordinates (=> GPCODT)
C            ELCO1T: spatial coordinates  (=> GPCO1T from set101t.f)
C
      DO IDIMET=1,NDIMET
       GPCODT(IDIMET)=0.0D0
       DO INODLT=1,NNODLT
        GPCODT(IDIMET)=GPCODT(IDIMET)+
     .                              ELCOIT(IDIMET,INODLT)*SHAPET(INODLT)
       END DO
      END DO
C
C**** IDENTIFIES HEAT FUNCTION
C
      QHEATT=0.0D0                     ! initialization
      DO IFUNCT=1,NFUNCT
       DO ISUBFT=1,NSUBFT
        ITYPL=DINT(HTLODT(1,ISUBFT,IFUNCT))
C
        IF(ITYPL.EQ.101.OR.ITYPL.EQ.102) THEN         ! linear
         LGRUPX=DINT(HTLODT(2,ISUBFT,IFUNCT))
         IF(LGRUPX.EQ.LGRUPT) THEN
          TIME0=HTLODT(3,ISUBFT,IFUNCT)
          TIMEF=HTLODT(4,ISUBFT,IFUNCT)
          RADIL=HTLODT(5,ISUBFT,IFUNCT)*HTLODT(5,ISUBFT,IFUNCT)
          DO IDIMET=1,NDIMET
           COORL(IDIMET)=HTLODT(5+IDIMET,ISUBFT,IFUNCT)
           VELOL(IDIMET)=HTLODT(5+NDIMET+IDIMET,ISUBFT,IFUNCT)
          END DO
          IF(TTIMET.GT.TIME0.AND.TTIMET.LT.TIMEF) THEN
           RADIU=0.0D0
           DO IDIMET=1,NDIMET
            RADIU=RADIU+
     .      (GPCODT(IDIMET)-VELOL(IDIMET)*(TTIMET-TIME0)-COORL(IDIMET))*
     .      (GPCODT(IDIMET)-VELOL(IDIMET)*(TTIMET-TIME0)-COORL(IDIMET))
           END DO
           RADIX=RADIL
c          IF(ITYPL.EQ.102) RADIX=1.498D0*RADIL  ! input heat of 95%
           IF(RADIU.LE.RADIX) THEN
            FF=1.0D0                             ! constant distribution
            IF(ITYPL.EQ.102)                     ! Gaussian distribution
     .       FF=2.0D0*DEXP(-2.0D0*RADIU/RADIL)
            QHEATT=QHEATT+HTLODT(6+2*NDIMET,ISUBFT,IFUNCT)*FF
           END IF
          END IF
         END IF     ! lgrupx=lgrupt
        ELSE IF(ITYPL.EQ.103.OR.ITYPL.EQ.104) THEN    ! circular
         LGRUPX=DINT(HTLODT(2,ISUBFT,IFUNCT))
         IF(LGRUPX.EQ.LGRUPT) THEN
          TIME0=HTLODT(3,ISUBFT,IFUNCT)
          TIMEF=HTLODT(4,ISUBFT,IFUNCT)
          RADIL=HTLODT(5,ISUBFT,IFUNCT)*HTLODT(5,ISUBFT,IFUNCT)
          DO IDIMET=1,NDIMET
           COORC(IDIMET)=HTLODT(5+IDIMET,ISUBFT,IFUNCT)
           COORL(IDIMET)=HTLODT(5+NDIMET+IDIMET,ISUBFT,IFUNCT)
           VELOL(IDIMET)=HTLODT(5+2*NDIMET+IDIMET,ISUBFT,IFUNCT)
          END DO
          IF(TTIMET.GT.TIME0.AND.TTIMET.LT.TIMEF) THEN
           X1MOD=0.0D0                 ! local system
           DO IDIMET=1,NDIMET
            X1MOD=X1MOD+
     .       (COORL(IDIMET)-COORC(IDIMET))*(COORL(IDIMET)-COORC(IDIMET))
           END DO
           X1MOD=DSQRT(X1MOD)
           DO IDIMET=1,NDIMET
            X1(IDIMET)=(COORL(IDIMET)-COORC(IDIMET))/X1MOD
           END DO
           X3MOD=0.0D0
           DO IDIMET=1,NDIMET
            X3MOD=X3MOD+VELOL(IDIMET)*VELOL(IDIMET)
           END DO
           X3MOD=DSQRT(X3MOD)
           DO IDIMET=1,NDIMET
            X3(IDIMET)=VELOL(IDIMET)/X3MOD
           END DO
           X2MOD=0.0D0
           DO IDIMA=1,NDIMET           ! compute normal as cross product
            IDIMB=IDIMA+1-IDIMA/3*3
            IDIMC=IDIMB+1-IDIMB/3*3
            X2(IDIMA)=X3(IDIMB)*X1(IDIMC)-X1(IDIMB)*X3(IDIMC)
            X2MOD=X2MOD+X2(IDIMA)*X2(IDIMA)
           END DO
           DO IDIMET=1,NDIMET
            X2(IDIMET)=X2(IDIMET)/X2MOD
           END DO
           XBI(1)=X1MOD*DCOS(X3MOD*(TTIMET-TIME0))
           XBI(2)=X1MOD*DSIN(X3MOD*(TTIMET-TIME0))
           XBI(3)=0.0D0
           DO IDIMET=1,NDIMET
            TB(1,IDIMET)=X1(IDIMET)
            TB(2,IDIMET)=X2(IDIMET)
            TB(3,IDIMET)=X3(IDIMET)
           END DO
           DO IDIMET=1,NDIMET
            XB(IDIMET)=0.0D0
            DO JDIMET=1,NDIMET
             XB(IDIMET)=XB(IDIMET)+TB(JDIMET,IDIMET)*XBI(JDIMET)
            END DO
            XB(IDIMET)=XB(IDIMET)+COORC(IDIMET)
           END DO
           RADIU=0.0D0
           DO IDIMET=1,NDIMET
            RADIU=RADIU+
     .           (GPCODT(IDIMET)-XB(IDIMET))*(GPCODT(IDIMET)-XB(IDIMET))
           END DO
           RADIX=RADIL
c          IF(ITYPL.EQ.104) RADIX=1.498D0*RADIL  ! input heat of 95%
           IF(RADIU.LE.RADIX) THEN
            FF=1.0D0                             ! constant distribution
            IF(ITYPL.EQ.104)                     ! Gaussian distribution
     .       FF=2.0D0*DEXP(-2.0D0*RADIU/RADIL)
            QHEATT=QHEATT+HTLODT(6+3*NDIMET,ISUBFT,IFUNCT)*FF
           END IF
          END IF
         END IF     ! lgrupx=lgrupt
        ELSE IF(ITYPL.GE.105.AND.ITYPL.LE.120) THEN
         CALL RUNENDT('ERROR: ITYPL GE 105 & LE 120 NOT IMPLEM. YET')
        ELSE IF(ITYPL.EQ.121.OR.ITYPL.EQ.122) THEN    ! linear-ellipse
         LGRUPX=DINT(HTLODT(2,ISUBFT,IFUNCT))
         IF(LGRUPX.EQ.LGRUPT) THEN
          TIME0=HTLODT(3,ISUBFT,IFUNCT)
          TIMEF=HTLODT(4,ISUBFT,IFUNCT)
          DO IDIMET=1,NDIMET
           RADIA(IDIMET)=HTLODT(4+IDIMET,ISUBFT,IFUNCT)*
     .                   HTLODT(4+IDIMET,ISUBFT,IFUNCT)
           COORL(IDIMET)=HTLODT(4+NDIMET+IDIMET,ISUBFT,IFUNCT)
           VELOL(IDIMET)=HTLODT(4+2*NDIMET+IDIMET,ISUBFT,IFUNCT)
          END DO
          IF(TTIMET.GT.TIME0.AND.TTIMET.LT.TIMEF) THEN
           RADIU=0.0D0
           DO IDIMET=1,NDIMET
            RADIU=RADIU+
     .      (GPCODT(IDIMET)-VELOL(IDIMET)*(TTIMET-TIME0)-COORL(IDIMET))*
     .      (GPCODT(IDIMET)-VELOL(IDIMET)*(TTIMET-TIME0)-COORL(IDIMET))/
     .      RADIA(IDIMET)
           END DO
           RADIX=1.0D0
c                                                ! to be revised !!!
           IF(ITYPL.EQ.122) RADIX=1.498D0        ! input heat of 95%
           IF(RADIU.LE.RADIX) THEN
            FF=1.0D0                             ! constant distribution
            IF(ITYPL.EQ.122)                     ! Gaussian distribution
     .       FF=2.0D0*DEXP(-2.0D0*RADIU)
            QHEATT=QHEATT+HTLODT(5+3*NDIMET,ISUBFT,IFUNCT)*FF
           END IF
          END IF
         END IF     ! lgrupx=lgrupt
        ELSE IF(ITYPL.EQ.123.OR.ITYPL.EQ.124) THEN    ! circular-ellipse
         LGRUPX=DINT(HTLODT(2,ISUBFT,IFUNCT))
         IF(LGRUPX.EQ.LGRUPT) THEN
          TIME0=HTLODT(3,ISUBFT,IFUNCT)
          TIMEF=HTLODT(4,ISUBFT,IFUNCT)
          DO IDIMET=1,NDIMET
           RADIA(IDIMET)=HTLODT(4+IDIMET,ISUBFT,IFUNCT)*
     .                   HTLODT(4+IDIMET,ISUBFT,IFUNCT)
           COORC(IDIMET)=HTLODT(4+IDIMET,ISUBFT,IFUNCT)
           COORL(IDIMET)=HTLODT(4+2*NDIMET+IDIMET,ISUBFT,IFUNCT)
           VELOL(IDIMET)=HTLODT(4+3*NDIMET+IDIMET,ISUBFT,IFUNCT)
          END DO
          IF(TTIMET.GT.TIME0.AND.TTIMET.LT.TIMEF) THEN
           X1MOD=0.0D0                 ! local system
           DO IDIMET=1,NDIMET
            X1MOD=X1MOD+
     .       (COORL(IDIMET)-COORC(IDIMET))*(COORL(IDIMET)-COORC(IDIMET))
           END DO
           X1MOD=DSQRT(X1MOD)
           DO IDIMET=1,NDIMET
            X1(IDIMET)=(COORL(IDIMET)-COORC(IDIMET))/X1MOD
           END DO
           X3MOD=0.0D0
           DO IDIMET=1,NDIMET
            X3MOD=X3MOD+VELOL(IDIMET)*VELOL(IDIMET)
           END DO
           X3MOD=DSQRT(X3MOD)
           DO IDIMET=1,NDIMET
            X3(IDIMET)=VELOL(IDIMET)/X3MOD
           END DO
           X2MOD=0.0D0
           DO IDIMA=1,NDIMET           ! compute normal as cross product
            IDIMB=IDIMA+1-IDIMA/3*3
            IDIMC=IDIMB+1-IDIMB/3*3
            X2(IDIMA)=X3(IDIMB)*X1(IDIMC)-X1(IDIMB)*X3(IDIMC)
            X2MOD=X2MOD+X2(IDIMA)*X2(IDIMA)
           END DO
           DO IDIMET=1,NDIMET
            X2(IDIMET)=X2(IDIMET)/X2MOD
           END DO
           XBI(1)=X1MOD*DCOS(X3MOD*(TTIMET-TIME0))
           XBI(2)=X1MOD*DSIN(X3MOD*(TTIMET-TIME0))
           XBI(3)=0.0D0
           DO IDIMET=1,NDIMET
            TB(1,IDIMET)=X1(IDIMET)
            TB(2,IDIMET)=X2(IDIMET)
            TB(3,IDIMET)=X3(IDIMET)
           END DO
           DO IDIMET=1,NDIMET
            XB(IDIMET)=0.0D0
            DO JDIMET=1,NDIMET
             XB(IDIMET)=XB(IDIMET)+TB(JDIMET,IDIMET)*XBI(JDIMET)
            END DO
            XB(IDIMET)=XB(IDIMET)+COORC(IDIMET)
           END DO
           RADIU=0.0D0
           DO IDIMET=1,NDIMET
            RADIU=RADIU+
     .          (GPCODT(IDIMET)-XB(IDIMET))*(GPCODT(IDIMET)-XB(IDIMET))/
     .          RADIA(IDIMET)
           END DO
           RADIX=1.0D0
c                                                ! to be revised !!!
           IF(ITYPL.EQ.124) RADIX=1.498D0        ! input heat of 95%
           IF(RADIU.LE.RADIX) THEN
            FF=1.0D0                             ! constant distribution
            IF(ITYPL.EQ.124)                     ! Gaussian distribution
     .       FF=2.0D0*DEXP(-2.0D0*RADIU)
            QHEATT=QHEATT+HTLODT(5+4*NDIMET,ISUBFT,IFUNCT)*FF
           END IF
          END IF
         END IF     ! lgrupx=lgrupt
        ELSE IF(ITYPL.GE.125.AND.ITYPL.LE.140) THEN
         CALL RUNENDT('ERROR: ITYPL GE 125 & LE 140 NOT IMPLEM. YET')
        END IF
C
       END DO       ! isubft=1,nsubft
      END DO        ! ifunct=1,nfunct
C
      RETURN
      END
