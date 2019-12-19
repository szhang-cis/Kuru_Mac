      SUBROUTINE MICROST(PROPST,TGAUAT,TGAUSX,TGAUST,TGAUIT,
     .                   DTEMPX,DTEMPT,
     .                   BASMM ,BASCC ,BASKK ,SOUR2T,DSOURT,SOUC2T,
     .                   SOUX2T,
     .                   ALPHAM)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE PHASE-CHANGE FUNCTION, DENSITY,
C     CAPACITY & CONDUCTIVITY ACCORDING WITH A MICROSTRUCTURAL MODEL
C
C     Notes:
C
C     This routine should be modified to deal with filling materials.
C
C     For advective problems (ICONVT=1), the temperature derivative of
C     f_pc must be evaluated as better as possible because it takes
C     place not only on the Jacobian but also on the residual (this is
C     also true for macroscopic problems). Therefore, the option
C     IFPCDT=3 should be used for computing such derivative in routines
C     micros*.f. Other possibility (not implemented yet) is to compute
C     this derivative in two forms: one for the Jacobian and other for
C     the residual (like in an advectionless problem).
C
C***********************************************************************
C
C             IPCFO=0 TOTAL PHASE-CHANGE FORMULATION
C                     f_pc(T) is bijective
C
C                     IPCMO=1 LINEAR
C                     IPCMO=2 IDEM 1 WITH f_s
C                     IPCMO=3 PARABOLIC
C                     IPCMO=4 CUBIC
C                     IPCMO=5 SCHEIL'S EQUATION
C                     IPCMO=6 LEVER'S EQUATION
C                     IPCMO=7 ....
C
C             IPCFO=1 RATE PHASE-CHANGE FORMULATION
C                     f_pc(\alpha_m_k,T) where:
C                     \alpha_m_k=microscopical variables
C                     => the computation of f^._pc is necessary
C
C                     IPCMO=1 MODEL ...
C                     IPCMO=2 MODEL ...
C
C***********************************************************************
C
C     Index of variables
C
C     BASMM = Density
C     BASCC = Capacity coefficient
C     BASKK = Conductivity
C     SOUR2T= L*Phase-change function (total)
C     DSOURT= L*Phase-change function rate (total)
C     SOUC2T= L*d(Phase-change function)/dT (total)
C
C     ALPHAM= array of microstructural (microscopical) variables
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nued_om.f'   ! thermal-microstructural
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      COMMON/CONVE3/ICONV2
C
      DIMENSION PROPST(*), ALPHAM(NHISTM)
      DIMENSION BASKK(NDIMETO)
C
      DIMENSION TSOE1(5), TSOE2(5)
      DIMENSION TSOC1(5), TSAUX(5)
C
C**** INITIALIZATION
C
c     SOUX2T=SOUC2T                ! see fruc05t.f & stuc05t.f
C
      SOUR1T=0.0D0                 ! L*f_pc (total) at time t
      SOUR2T=0.0D0                 ! L*f_pc (total) at time t+dt
      SOUC2T=0.0D0                 ! L*(df_pc/dT) (total) at time t+dt
      DSOURT=0.0D0                 ! L*(\delta f_pc)
      IF(NPLAT.EQ.0) RETURN
C
      VELTOT=1.D-10
      VELA1T=DABS(DTEMPT)
C
      IF(IMICO.EQ.0) THEN
       TGAUSX=TGAUAT                           ! at time t
       TGAUAT=TGAUSX-DTEMPX*DTIMET             ! at time t-dt
      ENDIF
      IF(LINPC.EQ.1) TGAUSX=TGAUAT             ! linearized phase-change
C
      INUPC=0                      ! ALPHAM index for NPLAT > 1
      INDEX3=0                     ! APLUOT index for NPLAT > 1
      INUPM=0                      ! number of microstructural ph-ch
      DO IPLAT=1,NPLAT
       TSOE1(IPLAT)=0.0D0          ! each f_pc at time t+dt
       TSOE2(IPLAT)=0.0D0          ! each f_pc at time t
       TSOC1(IPLAT)=0.0D0
       IF(ICONVT.EQ.1) TSAUX(IPLAT)=0.0D0
       IPCFO=INT(VPLAT(IPLAT,4))
       IPCMO=INT(VPLAT(IPLAT,5))
C
       IF(IPCFO.EQ.0) THEN
C
C**** TOTAL PHASE-CHANGE FORMULATION (IPCFO=0)
C
        TEINF=VPLAT(IPLAT,1)             ! each T_solidus
        TESUP=VPLAT(IPLAT,2)             ! each T_liquidus
        HENER=VPLAT(IPLAT,3)             ! each latent heat
        DELEE=TESUP-TEINF                ! each (T_l-T_s)
C
C**** CALCULATES THE LATENT HEAT RELEASED
C
        IF(IPCMO.EQ.-1)                  ! no real phase-change
     .   CALL RUNENDT('ERROR: NO R_P_C WITH MICRO. MODELS NOT POSSIBLE')
C
        IF(IPCMO.EQ.0) THEN              ! f_pc linear
C
         IF(LINPC.EQ.0) THEN
          IF(TGAUSX.LE.TEINF) TSOE1(IPLAT)=0.0D0
          IF(TGAUSX.GT.TEINF.AND.TGAUSX.LE.TESUP)
     .     TSOE1(IPLAT)=1.0D0/DELEE*(TGAUSX-TEINF)
          IF(TGAUSX.GT.TESUP) TSOE1(IPLAT)=1.0D0
C
          IF(TGAUAT.LE.TEINF) TSOE2(IPLAT)=0.0D0
          IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
     .     TSOE2(IPLAT)=1.0D0/DELEE*(TGAUAT-TEINF)
          IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=1.0D0
C
          TSOE1(IPLAT)=TSOE1(IPLAT)*HENER
          TSOE2(IPLAT)=TSOE2(IPLAT)*HENER
C
          ICACOT=1
          IF(VELA1T.GT.VELTOT) THEN
           IF(IITERT.GT.0)
     .      TSOC1(IPLAT)=(TSOE1(IPLAT)-TSOE2(IPLAT))/(DTEMPT*DTIMET)
          ENDIF
         ELSE                                  ! linearized phase-change
          IF(TGAUSX.LE.TEINF.OR.TGAUSX.GT.TESUP) TSOC1(IPLAT)=0.0D0
          IF(TGAUSX.GT.TEINF.AND.TGAUSX.LE.TESUP)
     .     TSOC1(IPLAT)=1.0D0/DELEE*HENER
         ENDIF        ! linpc.eq.0
C
        ELSE
         IF(LINPC.EQ.0) THEN
          ILAHET=1
          CALL MICPHAS(TGAUSX,TGAUAT,ILAHET,
     .                 IPLAT,HENER,TSOE1,TSOE2,
     .                 TEINF,TESUP,DELEE,IPCMO,    0)
C
          ICACOT=1
          IF(VELA1T.GT.VELTOT) THEN
           IF(IITERT.GT.0)
     .      TSOC1(IPLAT)=(TSOE1(IPLAT)-TSOE2(IPLAT))/(DTEMPT*DTIMET)
          ENDIF
C
          IF(ICONVT.EQ.1) THEN               ! see stuc05t.f
           ILAHET=2
           ISINRT=1
           IF(ICONV2.NE.0) THEN
            CALL MICPHAS(TGAUSX,TGAUAT,ILAHET,
     .                   IPLAT,HENER,TSAUX,TSOE2,
     .                   TEINF,TESUP,DELEE,IPCMO,    0)
            IF(IITERT.GT.0)
     .       TSOC1(IPLAT)=TSAUX(IPLAT)
           ELSE
            if(iitert.gt.1) ILAHET=1
            CALL MICPHAS(TGAUSX,TGAUAT,ILAHET,
     .                   IPLAT,HENER,TSAUX,TSOE2,
     .                   TEINF,TESUP,DELEE,IPCMO,    0)
            if(iitert.gt.1)
     .       TSOC1(IPLAT)=(TSOE1(IPLAT)-TSOE2(IPLAT))/(DTEMPT*DTIMET)
           ENDIF
          ENDIF       ! iconvt.eq.1
         ELSE                                  ! linearized phase-change
          ILAHET=2
          ISINRT=1
          CALL MICPHAS(TGAUSX,TGAUAT,ILAHET,
     .                 IPLAT,HENER,TSAUX,TSOE2,
     .                 TEINF,TESUP,DELEE,IPCMO,    0)
          TSOC1(IPLAT)=TSAUX(IPLAT)
         ENDIF        ! linpc.eq.0
        ENDIF         ! ipcmo.eq.0
C
       ELSE           ! ipcfo.eq.0
C
C**** RATE PHASE-CHANGE FORMULATION (IPCFO=1) > MICROSCOPICAL MODELS
C
        INUPM=INUPM+1
        GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14), IPCMO
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 1
C
    1   CALL MICROS1(TGAUAT,TGAUSX,TGAUIT,DTEMPX,
     .                BASMM, BASCC, BASKK,
     .                TSOE2, TSOE1, TSOC1,
     .                IPLAT,
     .               ALPHAM, INUPC)
        GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 2
C
    2   CALL MICROS2(TGAUAT,TGAUSX,TGAUIT,DTEMPX,
     .                BASMM, BASCC, BASKK,
     .                TSOE2, TSOE1, TSOC1,
     .                IPLAT,
     .               ALPHAM, INUPC,INDEX3)
        GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 3
C
    3   CALL MICROS3(TGAUAT,TGAUSX,TGAUIT,DTEMPX,
     .                BASMM, BASCC, BASKK,
     .                TSOE2, TSOE1, TSOC1,
     .                IPLAT,
     .               ALPHAM, INUPC)
        GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 4
C
    4   CALL MICROS4(TGAUAT,TGAUSX,TGAUIT,DTEMPX,
     .                BASMM, BASCC, BASKK,
     .                TSOE2, TSOE1, TSOC1,
     .                IPLAT,
     .               ALPHAM, INUPC)
        GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 5
C
    5   CALL MICROS5(TGAUAT,TGAUSX,TGAUIT,DTEMPX,
     .                BASMM, BASCC, BASKK,
     .                TSOE2, TSOE1, TSOC1,
     .                IPLAT,
     .               ALPHAM, INUPC)
        GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 6
C
    6   CALL MICROS6(TGAUAT,TGAUSX,TGAUIT,DTEMPX,
     .                BASMM, BASCC, BASKK,
     .                TSOE2, TSOE1, TSOC1,
     .                IPLAT, INUPM,
     .               ALPHAM, INUPC)
        GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 7
C
    7   CALL MICROS7(TGAUAT,TGAUSX,TGAUIT,DTEMPX,
     .                BASMM, BASCC, BASKK,
     .                TSOE2, TSOE1, TSOC1,
     .                IPLAT,
     .               ALPHAM, INUPC)
        GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 8
C
    8   CALL MICROS8(TGAUAT,TGAUSX,TGAUIT,DTEMPX,
     .                BASMM, BASCC, BASKK,
     .                TSOE2, TSOE1, TSOC1,
     .                IPLAT,
     .               ALPHAM, INUPC)
        GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 9
C
    9   CALL MICROS9(TGAUAT,TGAUSX,TGAUIT,DTEMPX,
     .                BASMM, BASCC, BASKK,
     .                TSOE2, TSOE1, TSOC1,
     .                IPLAT,
     .               ALPHAM, INUPC)
        GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 10
C
   10   CALL MICROS10(TGAUAT,TGAUSX,TGAUIT,DTEMPX,
     .                 BASMM, BASCC, BASKK,
     .                 TSOE2, TSOE1, TSOC1,
     .                 IPLAT,
     .                ALPHAM, INUPC)
        GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 11
C
   11   CALL MICROS11(TGAUAT,TGAUSX,TGAUIT,DTEMPX,
     .                 BASMM, BASCC, BASKK,
     .                 TSOE2, TSOE1, TSOC1,
     .                 IPLAT, INUPM,
     .                ALPHAM, INUPC)
        GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 12
C
   12   CALL MICROS12(TGAUAT,TGAUSX,TGAUIT,DTEMPX,
     .                 BASMM, BASCC, BASKK,
     .                 TSOE2, TSOE1, TSOC1,
     .                 IPLAT, INUPM,
     .                ALPHAM, INUPC)
        GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 13
C
   13   CALL RUNENDT('ERROR: MODEL 13 NOT IMPLEMENTED YET')
        GO TO 100
C
C**** CONSIDERS MICROSTRUCTURAL MODEL 14
C
   14   CALL MICROS14(TGAUAT,TGAUSX,TGAUIT,DTEMPX,
     .                 BASMM, BASCC, BASKK,
     .                 TSOE2, TSOE1, TSOC1,
     .                 IPLAT, INUPM,
     .                ALPHAM, INUPC)
        GO TO 100
C
  100   CONTINUE
C
       ENDIF        ! ipcfo.eq.0
C
      ENDDO         ! iplat=1,nplat
C
      DO IPLAT=1,NPLAT
       SOUR1T=SOUR1T+TSOE2(IPLAT)
       SOUR2T=SOUR2T+TSOE1(IPLAT)
       SOUC2T=SOUC2T+TSOC1(IPLAT)
      ENDDO
C
C**** PHASE-CHANGE FUNCTION RATE (*DTIMET)
C
      DSOURT=SOUR2T-SOUR1T
C
C**** DEALS WITH LINEARIZED PHASE-CHANGE
C
      IF(LINPC.EQ.1) THEN
       SOUR2T=SOUC2T
       SOUR1T=0.0D0
       DSOURT=SOUR2T-SOUR1T
      ENDIF
C
      RETURN
      END
